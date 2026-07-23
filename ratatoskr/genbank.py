from datetime import datetime
from importlib import resources
from io import BytesIO
import os
import re
import subprocess
import sys
import time

import asyncio
import aiohttp

from Bio import Entrez, SeqIO
from loguru import logger
import tqdm
import polars as pl


from ratatoskr.misc import get_haves_and_have_nots, fetch_data, tidy_genome_dir, tidy_16S_dir
from ratatoskr.utils import get_credentials, unzip_file, make_dir, move_and_rename, delete_thing
from ratatoskr.outputs import output_metadata


class RateLimiter:
    def __init__(self, max_rate: float):
        self.max_rate = max_rate          # requests per second
        self.interval = 1 / max_rate
        self._lock = asyncio.Lock()
        self._last_call = 0

    async def acquire(self):
        async with self._lock:
            now = time.monotonic()
            wait = self.interval - (now - self._last_call)
            if wait > 0:
                await asyncio.sleep(wait)
            self._last_call = time.monotonic()


async def esearch(term, db="nuccore", retmax=10**4, email=None, api_key=None, rate_limiter=None):
    await rate_limiter.acquire()
    params = {"db": db, "term": term, "retmax": str(retmax), "retmode": "xml", "email": email, "api_key": api_key, "usehistory": "y"}
    async with aiohttp.ClientSession() as session:
        async with session.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params=params) as resp:
            text = await resp.read()
            return BytesIO(text)
        

async def esummary(query_key, webenv, db="nuccore", email=None, api_key=None, rate_limiter=None):
    await rate_limiter.acquire()
    params = {"db": db, "query_key":str(query_key), "webenv": str(webenv), "retmode": "xml", "email": email, "api_key": api_key}
    async with aiohttp.ClientSession() as session:
        async with session.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params=params) as resp:
            text = await resp.read()
            return BytesIO(text)
    

def get_genbank_api_info(dev_mode=False):
    """
    Set up the LPSN client.
    """
    email, _, api_key = get_credentials(email=True, password=False, api_key=True, api_being_accessed="genbank", dev_mode=dev_mode)
    Entrez.email = email
    Entrez.api_key = api_key
    return email, api_key


async def request_ncbi_taxon_ids(query_terms, api_key, sem, pbar=None):
    headers = {
        "accept": "application/json",
        "api-key": api_key
    }
    session = aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60000),connector=aiohttp.TCPConnector(limit=9))
    sem = asyncio.Semaphore(9)
    urls = [f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{"%2C".join(subset).replace(" ", "%20")}?returned_content=TAXIDS&page_size=1000" for subset in [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]]
    query_lengths = [len(subset) for subset in [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]]
    parrallel_tasks = [fetch_data(url, session, headers, sem, 'taxonomy_nodes', pbar, query_length) for url, query_length in zip(urls, query_lengths)]
    data_list = await asyncio.gather(*parrallel_tasks)
    await session.close()
    return [item for sublist in data_list for item in sublist if sublist is not None]
    

def retrieve_ncbi_taxon_ids(lpsn_types, api_key):

    logger.info("Retrieving NCBI Taxon IDs for LPSN type strains.")
    
    query_terms = []
    has_ncbi_taxid, missing_ncbi_taxid = get_haves_and_have_nots(lpsn_types, "species_ncbi_tax_id")
    missing_ncbi_taxid.extend([x for x in has_ncbi_taxid if x.parent_subspecies is not None])
    has_ncbi_taxid = [x for x in has_ncbi_taxid if x.parent_subspecies is None]
    for type_strain in missing_ncbi_taxid:
        if type(type_strain.binomial_synonyms) != list:
            type_strain.binomial_synonyms = type_strain.binomial_synonyms.split(",")
        if type_strain.parent_subspecies is not None:
            type_strain.binomial_synonyms.append(type_strain.parent_species)
        query_terms.extend(type_strain.binomial_synonyms)

    pbar = tqdm.tqdm(total=len(query_terms), desc="Requesting NCBI Taxon IDs", unit="query", ncols=100, colour="magenta")
    data_list = asyncio.run(request_ncbi_taxon_ids(query_terms,  api_key, sem=None, pbar=pbar))
    if len (data_list) == 0 or all([x.get("taxonomy") is None for x in data_list]):
        logger.info("No NCBI Taxon IDs retrieved from GenBank. Continuing")
        return has_ncbi_taxid + missing_ncbi_taxid

    pbar.close()
    data_dict = pl.DataFrame(data_list).explode('query').select(pl.col('query'), pl.col('taxonomy').struct.field("*")).drop_nulls().rows_by_key('query', unique=True)
    
    for type_strain in tqdm.tqdm(missing_ncbi_taxid, desc="Processing NCBI Taxon IDs", unit="type strain", ncols=100, colour="magenta"):
        if type(type_strain.binomial_synonyms) != list:
            type_strain.binomial_synonyms = type_strain.binomial_synonyms.split(",")
        hits = []
        for synonym in type_strain.binomial_synonyms:
            if data_dict.get(synonym) is not None:
                hits.append(*data_dict.get(synonym))
        type_strain.species_ncbi_tax_id = hits
        if len(hits) > 1:
            logger.debug(f"Multiple Taxon IDs found for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species}. Will carry forward using all matches.")
        if len(hits) == 0:
            logger.debug(f"No Taxon ID found for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species}")
            type_strain.species_ncbi_tax_id = []

    return has_ncbi_taxid + missing_ncbi_taxid


async def search_16s(batch, email, api_key, pbar, rate_limiter, db_to_search):
    max_retries = 5
    batch_df = pl.DataFrame([{
        "strain_names": [type_strain.type_names] if isinstance(type_strain.type_names, str) else type_strain.type_names,
        "ncbi_tax_id": type_strain.strain_ncbi_tax_id + type_strain.species_ncbi_tax_id, 
        "parent_species_id": type_strain.parent_species_id, 
        "parent_subspecies_id": type_strain.parent_subspecies_id} for type_strain in batch], 
                schema={"strain_names": pl.List(pl.String), 
                        "ncbi_tax_id": pl.List(pl.String), 
                        "parent_species_id": pl.Int64, 
                        "parent_subspecies_id": pl.Int64}
    ).lazy().explode(["ncbi_tax_id"])
    ids = batch_df.select(pl.col("ncbi_tax_id")).unique().collect().to_series().drop_nulls().to_list()
    have_rRNA = []
    missing_rRNA = []
    if db_to_search == "refseq":
        term = "(" + " OR ".join([f"txid{i}[Organism:exp]" for i in ids if i is not None]) + ")" + " AND (33175[BioProject] OR 33317[BioProject])"
    elif db_to_search == "genbank":
        term = "(" + " OR ".join([f"txid{i}[Organism:exp]" for i in ids if i is not None]) + ")" + " AND 16S[Title]"
    for attempt in range(1, max_retries + 1):
        try:
            search_result_xml = await esearch(db="nucleotide", term=term, email = email, retmax=10**6, api_key=api_key, rate_limiter=rate_limiter)
            search_result = Entrez.read(search_result_xml)
            webenv = search_result["WebEnv"]
            query_key = search_result["QueryKey"]
            break
        except Exception as e:
            if attempt == max_retries:
                logger.error(f"Failed to retrieve search results for batch after {max_retries} attempts. Error: {e}")
                pbar.update(len(batch))
                return {"hits": have_rRNA, "misses" : batch}
            await asyncio.sleep(2 ** attempt)
    if search_result.get("Count") != "0":
        for attempt in range(1, max_retries + 1):
            try:
                records = await esummary(query_key, webenv, db="nuccore", api_key=api_key, email=email, rate_limiter=rate_limiter)
                records = Entrez.read(records)
                break
            except Exception as e:
                if attempt == max_retries:
                    logger.error(f"Failed to retrieve esummary results for batch after {max_retries} attempts. Error: {e}")
                    pbar.update(len(batch))
                    return {"hits": have_rRNA, "misses" : batch}
                await asyncio.sleep(2 ** attempt)

        records_df = (
            pl.LazyFrame(
                {'Length': int(x.get("Length")), 
                'Status': x.get("Status"), 
                'AccessionVersion': x.get("AccessionVersion"), 
                'CreateDate': x.get("CreateDate"), 
                'TaxId': str(int(x.get("TaxId"))),
                'Title': x.get("Title")}
                for x in records).filter(
                    (pl.col("Length") >= 750) &
                    (pl.col("Length") <= 2500) & 
                    (pl.col("Status") == "live") & 
                    (pl.col("AccessionVersion").str.len_chars() < 12)
                )
            )

        batch_df = (
            batch_df.with_columns(
                strain_names = pl.col("strain_names").list.concat(
                    pl.col("strain_names").list.eval(
                        pl.element().str.replace_all(" ", "-", literal=True)
                    )
                ).list.concat(
                    pl.col("strain_names").list.eval(
                        pl.element().str.replace_all(" ", "_", literal=True)
                    )
                ).list.concat(
                        pl.col("strain_names").list.eval(
                            pl.element().str.replace_all(r" ", "", literal=True)
                    )
                ).list.concat(
                        pl.col("strain_names").list.eval(
                            pl.element().str.replace_all(r"[^a-zA-Z0-9]", "")
                    )
            )).with_columns(
               strain_names = pl.col("strain_names").list.eval(pl.element().map_elements(lambda x: f" {x} ", return_dtype=pl.String)),
            ).join(
                records_df.select(pl.col("TaxId").alias("ncbi_tax_id"), 
                                  pl.col("AccessionVersion"), pl.col("CreateDate"), pl.col("Title")), 
                                  on="ncbi_tax_id", 
                                  how="left"
                ).with_columns(
                    Title = pl.col("Title").str.replace_all(r"[^a-zA-Z0-9\s]", "")
                ).filter(
                    pl.col("Title").str.extract_many(pl.col("strain_names")).list.len() > 0
                ).sort(
                    "CreateDate", 
                    descending=False
                ).group_by(
                    "parent_species_id", 
                    "parent_subspecies_id"
                ).agg(
                    pl.col("AccessionVersion").first().alias("rRNA_acc")
                )
        ).collect()
        subspecies_lookup = (batch_df.filter(pl.col("parent_subspecies_id").is_not_null())
            .select("parent_subspecies_id", "rRNA_acc")
            .to_dict(as_series=False))
        subspecies_map = dict(zip(subspecies_lookup["parent_subspecies_id"],
                            subspecies_lookup["rRNA_acc"]))

        species_lookup = (
            batch_df.filter(
                pl.col("parent_subspecies_id").is_null(), 
                pl.col("parent_species_id").is_not_null()
            ).select(
                "parent_species_id", 
                "rRNA_acc"
            ).to_dict(as_series=False)
        )
        species_map = dict(zip(species_lookup["parent_species_id"], species_lookup["rRNA_acc"]))

        for type_strain in batch:
            if type_strain.parent_subspecies_id is not None:
                hit = subspecies_map.get(type_strain.parent_subspecies_id)
            else:
                hit = species_map.get(type_strain.parent_species_id)

            if hit is not None:
                type_strain.rRNA_acc = hit.split(".")[0]
                if db_to_search == "refseq": 
                    type_strain.rRNA_info = {"source": "RefSeq", "note": "RefSeq reference 16S rRNA", "98.5%_match": None, "95%_match": None}
                if db_to_search == "genbank":
                    type_strain.rRNA_info = {"source": "GenBank", "note": "Retrieved from GenBank", "98.5%_match": None, "95%_match": None}
                have_rRNA.append(type_strain)
                pbar.update(1)
            else:
                missing_rRNA.append(type_strain)
    else:
        missing_rRNA.extend(batch)
    return {"hits": have_rRNA, "misses" : missing_rRNA}


async def search_16s_binomial(type_strain, email, api_key, pbar, rate_limiter):
    max_retries = 5
    binomial_search_term = " OR ".join([f'("{x}"[Organism])' for x in type_strain.binomial_synonyms])
    try:
        strain_search_term = " OR ".join(set([f'("{x.replace(" ", y)}"[Strain])' for x in type_strain.type_names for y in ["_", "-", "."]]))
    except Exception as e:
        logger.error(f"Error creating strain search term for type strain {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species}: {e}")
        strain_search_term = ""
        pbar.update(1)
        return type_strain
    for attempt in range(1, max_retries + 1):
        try:

            search_result_xml = await esearch(db="nucleotide", term=f"({binomial_search_term}) AND ({strain_search_term}) AND (16S[Title])", email = email, retmax=10**6, api_key=api_key, rate_limiter=rate_limiter)
            search_result = Entrez.read(search_result_xml)
            webenv = search_result["WebEnv"]
            query_key = search_result["QueryKey"]
            break
        except Exception as e:
            if attempt == max_retries:
                logger.error(f"Failed to retrieve search results for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species} after {max_retries} attempts. Error: {e}")
                logger.debug(f"Search result XML: {search_result_xml.getvalue().decode()}")
                pbar.update(1)
                return type_strain
            await asyncio.sleep(2 ** attempt)
    if search_result.get("Count") != "0":
        for attempt in range(1, max_retries + 1):
            try:
                records = await esummary(query_key, webenv, db="nuccore", api_key=api_key, email=email, rate_limiter=rate_limiter)
                records = Entrez.read(records)
                record = [x for x in records if 750 <= int(x.get("Length")) <= 2500 and x.get("Status") == "live" and len(x.get("AccessionVersion")) < 12]
                if len(record) > 0:
                    record = sorted(record, key=lambda x: datetime.strptime(x.get("CreateDate"), "%Y/%m/%d"))[0].get("AccessionVersion")
                    type_strain.rRNA_acc = record.split(".")[0]
                    type_strain.rRNA_info = {"source": "GenBank", "note": "Retrieved from GenBank", "98.5%_match": None, "95%_match": None}
                    pbar.update(1)
                    return type_strain
                break
            except Exception as e:
                if attempt == max_retries:
                    logger.error(f"Failed to retrieve esummary results for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species} after {max_retries} attempts. Error: {e}")
                    logger.debug(f"Esummary result XML: {records}")
                    pbar.update(1)
                    return type_strain
                await asyncio.sleep(2 ** attempt)

    pbar.update(1)
    return type_strain


def batch_missing_rrna(missing_rRNA, max_ids=80):
    batches = []
    current_batch = []
    current_count = 0

    for item in missing_rRNA:
        if item.strain_ncbi_tax_id is None:
            item.strain_ncbi_tax_id = []
        if item.species_ncbi_tax_id is None:
            item.species_ncbi_tax_id = [] 
        item_count = (
            len(item.strain_ncbi_tax_id) +
            len(item.species_ncbi_tax_id)
        )

        # If a single item exceeds the limit, force it into its own batch
        if item_count > max_ids:
            if current_batch:
                batches.append(current_batch)
                current_batch = []
                current_count = 0
            batches.append([item])
            continue

        # Start a new batch if needed
        if current_count + item_count > max_ids:
            batches.append(current_batch)
            current_batch = []
            current_count = 0

        current_batch.append(item)
        current_count += item_count

    if current_batch:
        batches.append(current_batch)

    return batches


def run_mmseqs_against_SILVA(fasta_file, output_path, threads):
    type_sequences_db_path = output_path / "sequences" / "16S" / "tmp_16S_mmseqs_db"
    create_mmseqs_db(fasta_file, type_sequences_db_path, threads)
    SILVA_db_path = resources.files("ratatoskr.data").joinpath("SILVA_138.2_SSURef_NR99_tax_silva.mmseqs")
    mmseqs_out_db = output_path / "sequences" / "16S" / "16S_mmseqs_hits.db"
    mmseqs_out_tmp = output_path / "sequences" / "16S" / "16S_mmseqs_hits.tmp"
    mmseqs_out_result = output_path / "sequences" / "16S" / "16S_mmseqs_hits.tsv"
    mmseqs_cmd = [
        "mmseqs",
        "search",
        str(type_sequences_db_path),
        str(SILVA_db_path),
        str(mmseqs_out_db),
        str(mmseqs_out_tmp),
        "--search-type", "3",
        "--min-seq-id", "0.729",
        "--remove-tmp-files", "1",
        "--threads", str(threads),
        "-s", "7.5",
        "--cov-mode", "1", 
        "-c", "0.5"
    ]
    subprocess.run(mmseqs_cmd, check=True)

    convert_cmd = [
        "mmseqs",
        "convertalis",
        "-v", "1",
        str(type_sequences_db_path),
        str(SILVA_db_path),
        str(mmseqs_out_db),
        str(mmseqs_out_result),
        "--format-output",
        "query,target,pident,alnlen,evalue,bits,theader,qlen,tlen"
    ]
    subprocess.run(convert_cmd, check=True)


def create_mmseqs_db(input_fasta, db_out_path, threads):
    db_cmd = [
        "mmseqs",
        "createdb",
        "--threads", str(threads),
        "-v", "1",
        str(input_fasta),
        str(db_out_path), 
    ]
    subprocess.run(db_cmd, check=True)  
    


def process_mmseqs_results(lpsn_types, output_path, mmseqs_16S_check):
    mmseqs_output_file = output_path / "sequences" / "16S" / "16S_mmseqs_hits.tsv"

    if mmseqs_16S_check == "all":
        search_sources = ["GenBank", "RefSeq", "LPSN", "BacDive"]
    elif mmseqs_16S_check == "genbank":
        search_sources = ["GenBank"]

    mmseqs_results = (
        pl.read_csv(mmseqs_output_file, 
                    separator="\t", 
                    has_header=False, 
                    schema = {"qseqid": pl.Utf8, 
                                 "sseqid": pl.Utf8, 
                                 "pident": pl.Float64, 
                                 "length": pl.Float64, 
                                 "evalue": pl.Float64, 
                                 "bitscore": pl.Float64, 
                                 "stitle": pl.Utf8,
                                 "qlen": pl.Int64,
                                 "slen": pl.Int64}
            ).lazy().filter(
                pl.col("length") >= pl.min_horizontal("qlen", "slen") * 0.9
            ).with_columns(
                stitle = pl.col("stitle").str.replace_all("[", "", literal=True).str.replace_all("]", "", literal=True)
            ).with_columns(
                qseqid = pl.col("qseqid").str.split(".").list.get(0),
                species_97point2 = pl.when((pl.col("pident") >= 97.2) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) & 
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(1).list.get(0).str.split(" ").list.head(2).list.join(" ")
                                ).otherwise(None),
                genus_97point2 = pl.when(pl.col("pident") >= 97.2).then(
                                    pl.col("stitle").str.split(";").list.tail(2).list.get(0)
                                ).otherwise(None),
                family_97point2 = pl.when(pl.col("pident") >= 97.2).then(
                                    pl.col("stitle").str.split(";").list.tail(3).list.get(0)
                                ).otherwise(None),
                order_97point2 = pl.when(pl.col("pident") >= 97.2).then(
                                    pl.col("stitle").str.split(";").list.tail(4).list.get(0)
                                ).otherwise(None),
                species_90point1 = pl.when((pl.col("pident") >= 90.1) & (pl.col("pident") < 97.2) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(1).list.get(0).str.split(" ").list.head(2).list.join(" ")
                                ).otherwise(None),
                genus_90point1 = pl.when((pl.col("pident") >= 90.1) & (pl.col("pident") < 97.2)).then(
                                    pl.col("stitle").str.split(";").list.tail(2).list.get(0)
                                ).otherwise(None),
                family_90point1 = pl.when((pl.col("pident") >= 90.1) & (pl.col("pident") < 97.2)).then(
                                    pl.col("stitle").str.split(";").list.tail(3).list.get(0)
                                ).otherwise(None),
                order_90point1 = pl.when((pl.col("pident") >= 90.1) & (pl.col("pident") < 97.2)).then(
                                    pl.col("stitle").str.split(";").list.tail(4).list.get(0)
                                ).otherwise(None),
                species_80point1 = pl.when((pl.col("pident") >= 80.1) & (pl.col("pident") < 90.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(1).list.get(0).str.split(" ").list.head(2).list.join(" ")
                                ).otherwise(None),
                genus_80point1 = pl.when((pl.col("pident") >= 80.1) & (pl.col("pident") < 90.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(2).list.get(0)
                                ).otherwise(None),
                family_80point1 = pl.when((pl.col("pident") >= 80.1) & (pl.col("pident") < 90.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(3).list.get(0)
                                ).otherwise(None),
                order_80point1 = pl.when((pl.col("pident") >= 80.1) & (pl.col("pident") < 90.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(4).list.get(0)
                                ).otherwise(None),
                species_72point9 = pl.when((pl.col("pident") >= 72.9) & (pl.col("pident") < 80.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(1).list.get(0).str.split(" ").list.head(2).list.join(" ")
                                ).otherwise(None),
                genus_72point9 = pl.when((pl.col("pident") >= 72.9) & (pl.col("pident") < 80.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(2).list.get(0)
                                ).otherwise(None),
                family_72point9 = pl.when((pl.col("pident") >= 72.9) & (pl.col("pident") < 80.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(3).list.get(0)
                                ).otherwise(None),
                order_72point9 = pl.when((pl.col("pident") >= 72.9) & (pl.col("pident") < 80.1) & 
                                           (pl.col("stitle").str.contains("(?i)uncultured") == False) & 
                                           (pl.col("stitle").str.contains(" sp.", literal=True) == False) &
                                           (pl.col("stitle").str.contains("(?i)metagenome") == False) 
                                ).then(
                                    pl.col("stitle").str.split(";").list.tail(4).list.get(0)
                                ).otherwise(None),
            ).select("qseqid", 
                     "species_97point2", "genus_97point2", "family_97point2", "order_97point2",
                     "species_90point1", "genus_90point1", "family_90point1", "order_90point1",
                     "species_80point1", "genus_80point1", "family_80point1", "order_80point1",
                     "species_72point9", "genus_72point9", "family_72point9", "order_72point9"
            ).group_by("qseqid").agg(
                pl.col("species_97point2").unique().alias("species_97point2"),
                pl.col("genus_97point2").unique().alias("genus_97point2"), 
                pl.col("family_97point2").unique().alias("family_97point2"), 
                pl.col("order_97point2").unique().alias("order_97point2"),
                pl.col("species_90point1").unique().alias("species_90point1"),
                pl.col("genus_90point1").unique().alias("genus_90point1"),
                pl.col("family_90point1").unique().alias("family_90point1"),
                pl.col("order_90point1").unique().alias("order_90point1"),
                pl.col("species_80point1").unique().alias("species_80point1"),
                pl.col("genus_80point1").unique().alias("genus_80point1"),
                pl.col("family_80point1").unique().alias("family_80point1"),
                pl.col("order_80point1").unique().alias("order_80point1"),
                pl.col("species_72point9").unique().alias("species_72point9"),
                pl.col("genus_72point9").unique().alias("genus_72point9"),
                pl.col("family_72point9").unique().alias("family_72point9"),
                pl.col("order_72point9").unique().alias("order_72point9")
            ).with_columns(
                species_97point2 = pl.col("species_97point2").list.drop_nulls(),
                genus_97point2 = pl.col("genus_97point2").list.drop_nulls(),
                family_97point2 = pl.col("family_97point2").list.drop_nulls(),
                order_97point2 = pl.col("order_97point2").list.drop_nulls(),
                species_90point1 = pl.col("species_90point1").list.drop_nulls(),
                genus_90point1 = pl.col("genus_90point1").list.drop_nulls(),
                family_90point1 = pl.col("family_90point1").list.drop_nulls(),
                order_90point1 = pl.col("order_90point1").list.drop_nulls(),
                species_80point1 = pl.col("species_80point1").list.drop_nulls(),
                genus_80point1 = pl.col("genus_80point1").list.drop_nulls(),
                family_80point1 = pl.col("family_80point1").list.drop_nulls(),
                order_80point1 = pl.col("order_80point1").list.drop_nulls(),
                species_72point9 = pl.col("species_72point9").list.drop_nulls(),
                genus_72point9 = pl.col("genus_72point9").list.drop_nulls(),
                family_72point9 = pl.col("family_72point9").list.drop_nulls(),
                order_72point9 = pl.col("order_72point9").list.drop_nulls()
            )
    ).lazy()

    types_df = (
            pl.DataFrame([{
            "order": type_strain.parent_order,
            "family": type_strain.parent_family,
            "genus": list(set([x.split(" ")[0] for x in type_strain.binomial_synonyms]).union(set([type_strain.parent_genus]))),
            "species": list(set(type_strain.binomial_synonyms).union(set([type_strain.parent_species]))),
            "rRNA_acc": type_strain.rRNA_acc.split(".")[0],
            "source": type_strain.rRNA_info.get("source", None) if type_strain.rRNA_info is not None else None
                } for type_strain in lpsn_types if type_strain.rRNA_acc is not None],
        ).lazy().join(
            mmseqs_results, 
            left_on="rRNA_acc", 
            right_on="qseqid", 
            how="left"
        ).with_columns(
            result_97point2 = pl.when(
                            (pl.col("order_97point2").is_null())
                        ).then(
                            pl.lit("No hits")
                        ).when(
                          (pl.col("species").list.set_intersection("species_97point2").list.len() > 0)
                        ).then(
                            pl.lit("Species match")
                        ).when(
                            (pl.col("genus").list.set_intersection("genus_97point2").list.len() > 0)
                        ).then(
                            pl.lit("Genus match")
                        ).when(
                            (pl.col("family").cast(pl.List(pl.String)).list.set_intersection("family_97point2").list.len() > 0)
                        ).then(
                            pl.lit("Family match")
                        ).when(
                            (pl.col("order").cast(pl.List(pl.String)).list.set_intersection("order_97point2").list.len() > 0)
                        ).then(
                            pl.lit("Order match")
                        ).otherwise(
                            pl.lit('Hits but no matches')
                        )
        ).with_columns(
            result_90point1 = pl.when(
                            (pl.col("order_90point1").is_null())
                        ).then(
                            pl.lit("No hits")
                        ).when(
                          (pl.col("species").list.set_intersection("species_90point1").list.len() > 0)
                        ).then(
                            pl.lit("Species match")
                        ).when(
                            (pl.col("genus").list.set_intersection("genus_90point1").list.len() > 0)
                        ).then(
                            pl.lit("Genus match")
                        ).when(
                            (pl.col("family").cast(pl.List(pl.String)).list.set_intersection("family_90point1").list.len() > 0)
                        ).then(
                            pl.lit("Family match")
                        ).when(
                            (pl.col("order").cast(pl.List(pl.String)).list.set_intersection("order_90point1").list.len() > 0)
                        ).then(
                            pl.lit("Order match")
                        ).otherwise(
                            pl.lit('Hits but no matches')
                        )
        ).with_columns(
            result_80point1 = pl.when(
                            (pl.col("order_80point1").is_null())
                        ).then(
                            pl.lit("No hits")
                        ).when(
                          (pl.col("species").list.set_intersection("species_80point1").list.len() > 0)
                        ).then(
                            pl.lit("Species match")
                        ).when(
                            (pl.col("genus").list.set_intersection("genus_80point1").list.len() > 0)
                        ).then(
                            pl.lit("Genus match")
                        ).when(
                            (pl.col("family").cast(pl.List(pl.String)).list.set_intersection("family_80point1").list.len() > 0)
                        ).then(
                            pl.lit("Family match")
                        ).when(
                            (pl.col("order").cast(pl.List(pl.String)).list.set_intersection("order_80point1").list.len() > 0)
                        ).then(
                            pl.lit("Order match")
                        ).otherwise(
                            pl.lit('Hits but no matches')
                        )
        ).with_columns(
            result_72point9 = pl.when(
                            (pl.col("order_72point9").is_null())
                        ).then(
                            pl.lit("No hits")
                        ).when(
                          (pl.col("species").list.set_intersection("species_72point9").list.len() > 0)
                        ).then(
                            pl.lit("Species match")
                        ).when(
                            (pl.col("genus").list.set_intersection("genus_72point9").list.len() > 0)
                        ).then(
                            pl.lit("Genus match")
                        ).when(
                            (pl.col("family").cast(pl.List(pl.String)).list.set_intersection("family_72point9").list.len() > 0)
                        ).then(
                            pl.lit("Family match")
                        ).when(
                            (pl.col("order").cast(pl.List(pl.String)).list.set_intersection("order_72point9").list.len() > 0)
                        ).then(
                            pl.lit("Order match")
                        ).otherwise(
                            pl.lit('Hits but no matches')
                        )
        ).with_columns(
            rRNA_note = pl.when(
                (pl.col("result_97point2") == "No hits") & 
                (pl.col("result_90point1") == "No hits") & 
                (pl.col("result_80point1") == "No hits") & 
                (pl.col("result_72point9") == "No hits")
            ).then(
                pl.lit("Warning: No hits vs SILVA")
            ).when(
                (pl.col("result_97point2").is_in(["Hits but no matches", "No hits"])) & 
                (pl.col("result_90point1").is_in(["Hits but no matches", "No hits"])) & 
                (pl.col("result_80point1").is_in(["Hits but no matches", "No hits"])) & 
                (pl.col("result_72point9").is_in(["Hits but no matches", "No hits"]))  
            ).then(
                pl.lit("Warning: SILVA hits found but none match the expected taxonomy")
            )
            .otherwise(
                pl.lit(None)
            )
        ).with_columns(
            result_97point2 = pl.when(~pl.col("source").is_in(search_sources)).then(None).otherwise(pl.col("result_97point2")),
            result_90point1 = pl.when(~pl.col("source").is_in(search_sources)).then(None).otherwise(pl.col("result_90point1")),
            result_80point1 = pl.when(~pl.col("source").is_in(search_sources)).then(None).otherwise(pl.col("result_80point1")),
            result_72point9 = pl.when(~pl.col("source").is_in(search_sources)).then(None).otherwise(pl.col("result_72point9")),
            rRNA_note = pl.when(~pl.col("source").is_in(search_sources)).then(None).otherwise(pl.col("rRNA_note"))
        ).select(
            "rRNA_acc", 
            "result_97point2", 
            "result_90point1", 
            "result_80point1", 
            "result_72point9", 
            "rRNA_note"
        )
    ).collect()        

    updated_types = []
    
    for i in lpsn_types:
        if i.rRNA_acc is not None:
            match = types_df.filter(pl.col("rRNA_acc") == i.rRNA_acc.split('.')[0]).to_dict(as_series=False)
            if len(match["rRNA_acc"]) > 0:
                i.rRNA_info = {
                    "97.2%_match": match["result_97point2"][0],
                    "90.1%_match": match["result_90point1"][0],
                    "80.1%_match": match["result_80point1"][0],
                    "72.9%_match": match["result_72point9"][0],
                    "note": match["rRNA_note"][0] if match["rRNA_note"][0] is not None else i.rRNA_info.get("note", None),
                    "source": i.rRNA_info.get("source", None)
                }
        updated_types.append(i)
          
    return updated_types

def fasta_to_dict(fasta_file):
    
    return {record.id.split('.')[0]: record for record in SeqIO.parse(fasta_file, "fasta")}


def get_tmp_genbank_16S_fasta(fasta, tmp_fasta, type_hit_list):

    genbank_accs = [x.rRNA_acc for x in type_hit_list if x.rRNA_info and (x.rRNA_info.get("source", None) == "GenBank")]
    genbank_records = fasta_to_dict(fasta)
    genbank_seqs = [genbank_records[acc] for acc in genbank_accs if acc in genbank_records]

    with open(tmp_fasta, "w") as f:
        SeqIO.write(genbank_seqs, f, "fasta")


    

def confirm_genbank_16S_hits_with_mmseqs(type_hit_list, output_path, threads, mmseqs_16S_check):

    if mmseqs_16S_check == "none":
        logger.info("Skipping MMSeqs2 vs SILVA check for 16S rRNAs.")
        return type_hit_list
    
    logger.info(f"Running MMSeqs2 vs SILVA on {mmseqs_16S_check} 16S rRNAs.")

    default_fasta = output_path / "sequences" / "16S" / "16S.fasta"

    if mmseqs_16S_check == "all":
        if len([x for x in type_hit_list if x.rRNA_acc is not None]) > 0:
            run_mmseqs_against_SILVA(default_fasta, output_path, threads)
            type_hit_list = process_mmseqs_results(type_hit_list, output_path, mmseqs_16S_check)
    if mmseqs_16S_check == "genbank":
        if len([x for x in type_hit_list if x.rRNA_info and (x.rRNA_info.get("source", None) == "GenBank")]) > 0:
            tmp_fasta = output_path / "sequences" / "16S" / "tmp_genbank_16S.fasta"
            get_tmp_genbank_16S_fasta(default_fasta, tmp_fasta, type_hit_list)
            run_mmseqs_against_SILVA(tmp_fasta, output_path, threads)
            type_hit_list = process_mmseqs_results(type_hit_list, output_path, mmseqs_16S_check)
    return type_hit_list


async def retrieve_missing_16S_from_refseq(missing_rRNA, max_terms, email, api_key, max_concurrency=10, pbar=None):
    
    semaphore = asyncio.Semaphore(max_concurrency)
    rate_limiter = RateLimiter(max_rate=8)
    batches = batch_missing_rrna(missing_rRNA, max_ids=max_terms)
    async def limited_search(batch, pbar, rate_limiter):
        async with semaphore:
            return await search_16s(batch, email, api_key, pbar=pbar, rate_limiter=rate_limiter, db_to_search="refseq")
    tasks = [limited_search(batch, pbar=pbar, rate_limiter=rate_limiter) for batch in batches]
    batch_results = await asyncio.gather(*tasks)
    hits = [x for r in batch_results for x in r["hits"]]
    misses = [x for r in batch_results for x in r["misses"]]

    return hits, misses
    

async def retrieve_missing_16S_from_genbank(missing_rRNA, max_terms, email, api_key, max_concurrency=10, pbar=None):
    
    semaphore = asyncio.Semaphore(max_concurrency)
    rate_limiter = RateLimiter(max_rate=8)
    batches = batch_missing_rrna(missing_rRNA, max_ids=max_terms)
    async def limited_search(batch, pbar, rate_limiter):
        async with semaphore:
            return await search_16s(batch, email, api_key, pbar=pbar, rate_limiter=rate_limiter, db_to_search="genbank")
    tasks = [limited_search(batch, pbar=pbar, rate_limiter=rate_limiter) for batch in batches]
    batch_results = await asyncio.gather(*tasks)
    hits = [x for r in batch_results for x in r["hits"]]
    misses = [x for r in batch_results for x in r["misses"]]

    if len(misses) > 0:
        async def limited_binomial_search(type_strain, pbar, rate_limiter):
            async with semaphore:
                return await search_16s_binomial(type_strain, email, api_key, pbar=pbar, rate_limiter=rate_limiter)
        binomial_tasks = [limited_binomial_search(type_strain, pbar=pbar, rate_limiter=rate_limiter) for type_strain in misses]
        binomial_results = await asyncio.gather(*binomial_tasks)

        hits += [x for x in binomial_results if x.rRNA_acc is not None]
        misses = [x for x in binomial_results if x.rRNA_acc is None]

    return hits, misses
    

def retrieve_missing_16S_info(lpsn_types, email, api_key, output_path):
    """
    Retrieve missing 16S rRNA gene information from NCBI (first RefSeq then GenBank) for the given LPSN type strains.
    """
    logger.info("Retrieving missing 16S rRNA gene information from NCBI.")
    Entrez.api_key = api_key
    has_rRNA, missing_rRNA = get_haves_and_have_nots(lpsn_types, "rRNA_acc")
    
    for type_strain in missing_rRNA:
        if not hasattr(type_strain, "parent_species_id"):
            logger.error(f"Type strain {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species} is missing parent_species_id attribute. This is required for 16S retrieval. Skipping this type strain for 16S retrieval.")
            logger.debug(type_strain)
    
    pbar = tqdm.tqdm(total=len(missing_rRNA), desc="Getting missing 16S rRNAs from RefSeq", unit="type strain", ncols=100, colour="magenta")
    refSeq_rRNA, still_missing = asyncio.run(retrieve_missing_16S_from_refseq(missing_rRNA, max_terms=99, email=email, api_key=api_key, max_concurrency=8, pbar=pbar))
    pbar.close()

    pbar = tqdm.tqdm(total=len(still_missing), desc="Retrieving missing 16S rRNA gene info from GenBank", unit="type strain", ncols=100, colour="magenta")
    genbank_rRNA, still_missing = asyncio.run(retrieve_missing_16S_from_genbank(still_missing, max_terms=99, email=email, api_key=api_key, max_concurrency=8, pbar=pbar))
    pbar.close()
 
    return has_rRNA + refSeq_rRNA + genbank_rRNA + still_missing


async def request_ncbi_genomes(query_terms, api_key, pbar, strain_names):
    headers = {
        "accept": "application/json",
        "api-key": api_key
    }
    logger.debug(f"Requesting genome information for {len(query_terms)} taxon IDs from GenBank.")
    session = aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60000), connector=aiohttp.TCPConnector(limit=9))
    sem = asyncio.Semaphore(9)

    chunks = [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]
    urls = [
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{'%2C'.join(subset).replace(' ', '%20')}/dataset_report?returned_content=COMPLETE&page_size=1000"
        for subset in chunks
    ]
    query_lengths = [len(subset) for subset in chunks]

    parrallel_tasks = [
        fetch_data(url, session, headers, sem, 'reports', pbar, query_length)
        for url, query_length in zip(urls, query_lengths)
    ]
    async with sem:
        data_list = await asyncio.gather(*parrallel_tasks)
    await session.close()
    
    data_list = [item for sublist in data_list for item in sublist if sublist is not None]

    data_list = [item for item in data_list 
                    if re.sub(r'[^a-zA-Z0-9]', '', item.get('organism', {}).get('infraspecific_names', {}).get('strain', '')) in strain_names or 
                    re.sub(r'[^a-zA-Z0-9]', '', item.get('organism', {}).get('infraspecific_names', {}).get('isolate', '')) in strain_names or
                    re.sub(r'[^a-zA-Z0-9]', '', item.get('assembly_info', {}).get('biosample', {}).get('strain', '')) in strain_names]

    return data_list


async def request_ncbi_checkm(query_terms, api_key, pbar):
    headers = {
        "accept": "application/json",
        "api-key": api_key
    }
    logger.debug(f"Requesting Checkm information for {len(query_terms)} genomes from GenBank.")
    session = aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60000),connector=aiohttp.TCPConnector(limit=9))
    sem = asyncio.Semaphore(9)
    urls = [f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{"%2C".join(subset).replace(" ", "%20")}/dataset_report?returned_content=COMPLETE&page_size=1000" for subset in [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]]
    query_lengths = [len(subset) for subset in [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]]
    parrallel_tasks = [fetch_data(url, session, headers, sem, 'reports', pbar, query_length) for url, query_length in zip(urls, query_lengths)]
    async with sem:
        data_list = await asyncio.gather(*parrallel_tasks)
    await session.close()
    return [item for sublist in data_list for item in sublist if sublist is not None]


def fill_missing_report_fields(ncbi_accession_report):

    for item in ncbi_accession_report:
        if 'atypical' not in item.get('assembly_info', {}):
            item['assembly_info']['atypical'] = {'is_atypical': None, 'warnings': []}
        if 'assembly_status' not in item.get('assembly_info', {}):
            item['assembly_info']['assembly_status'] = ''
        if 'paired_assembly' not in item.get('assembly_info', {}):
            item['assembly_info']['paired_assembly'] = {"accession": '',
                                                        "status": '',
                                                    }
        if 'genome_notes' not in item.get('assembly_info', {}):
            item['assembly_info']['genome_notes'] = []
        if 'checkm_info' not in item:
            item['checkm_info'] = {'checkm_species_tax_id': '', "completeness": '', "contamination": ''}
        if 'infraspecific_names' not in item.get('organism', {}):
            item['organism']['infraspecific_names'] = {'strain': '', 'isolate': ''}
        if 'strain' not in item.get('organism', {}).get('infraspecific_names', {}):
            item['organism']['infraspecific_names']['strain'] = ''
        if 'isolate' not in item.get('organism', {}).get('infraspecific_names', {}):
            item['organism']['infraspecific_names']['isolate'] = ''
        if 'type_material' not in item:
            item['type_material'] = {'type_display_text': '', 'type_label': ''}
        if 'type_display_text' not in item.get('type_material', {}):
            item['type_material']['type_display_text'] = ''
        if 'type_label' not in item.get('type_material', {}):
            item['type_material']['type_label'] = ''
        if 'average_nucleotide_identity' not in item:
            item['average_nucleotide_identity'] = {"submitted_ani_match": {"ani": None, 
                                                                           "organism_name": None, 
                                                                           "assembly": None,
                                                                           "assembly_coverage": None
                                                                           },
                                                    "best_ani_match": {"best_ani_match": None, 
                                                                        "organism_name": None, 
                                                                        "assembly": None,
                                                                        "assembly_coverage": None
                                                                        }
                                                    }
            
        for field in ['submitted_ani_match', 'best_ani_match']:
            if field not in item.get('average_nucleotide_identity', {}):
                item['average_nucleotide_identity'][field] = {"ani": None, "organism_name": None, "assembly": None, "assembly_coverage": None}
            for subfield in ['ani', 'organism_name', 'assembly', 'assembly_coverage']:
                if subfield not in item.get('average_nucleotide_identity', {}).get(field, {}):
                    item['average_nucleotide_identity'][field][subfield] = None

    return ncbi_accession_report


def convert_ncbi_accession_report_to_polars_df(ncbi_accession_report):

    needed_structs = [
        "accession", "organism", "assembly_info", "checkm_info",
        "type_material", "average_nucleotide_identity",
    ]

    ncbi_accession_report = fill_missing_report_fields(ncbi_accession_report)

    df = (
        pl.DataFrame(ncbi_accession_report).lazy(
            ).select(
                *needed_structs
            )
            .with_columns(
                pl.col('accession').cast(pl.String).str.replace('\\.[0-9]+$', '')
            ).unique(
                subset=["accession"], keep="first"
            ).unnest(
                'organism'
            ).with_columns(
                pl.col('assembly_info').struct.field('atypical')
            ).unnest(
                'atypical'
            ).with_columns(
                genome_notes = pl.col("assembly_info").struct.field("genome_notes"),
                checkm_completeness = pl.col('checkm_info').struct.field('completeness'),
                checkm_contamination = pl.col('checkm_info').struct.field('contamination'),
                biosample = pl.col('assembly_info').struct.field('biosample').struct.field('strain').cast(pl.String),
                assembly_level = pl.col('assembly_info').struct.field('assembly_level').cast(pl.String),
                isolate = pl.col('infraspecific_names').struct.field('isolate').cast(pl.String),
                strain = pl.col('infraspecific_names').struct.field('strain').cast(pl.String),
                type_display_text = pl.col('type_material').struct.field('type_display_text').cast(pl.String),
                type_label = pl.col('type_material').struct.field('type_label').cast(pl.String),
                ncbis_species_match_ani = pl.col('average_nucleotide_identity').struct.field('submitted_ani_match').struct.field('ani'),
                ncbis_species_match_name = pl.col('average_nucleotide_identity').struct.field('submitted_ani_match').struct.field('organism_name'),
                ncbis_species_match_assembly = pl.col('average_nucleotide_identity').struct.field('submitted_ani_match').struct.field('assembly'),
                ncbis_species_match_coverage = pl.col('average_nucleotide_identity').struct.field('submitted_ani_match').struct.field('assembly_coverage'),
                ncbis_best_match_ani = pl.col('average_nucleotide_identity').struct.field('best_ani_match').struct.field('ani'),
                ncbis_best_match_name = pl.col('average_nucleotide_identity').struct.field('best_ani_match').struct.field('organism_name'),
                ncbis_best_match_assembly = pl.col('average_nucleotide_identity').struct.field('best_ani_match').struct.field('assembly'),
                ncbis_best_match_coverage = pl.col('average_nucleotide_identity').struct.field('best_ani_match').struct.field('assembly_coverage'),
                is_atypical = pl.col('is_atypical').cast(pl.Boolean),
                paired_assembly_status = pl.col('assembly_info').struct.field('paired_assembly').struct.field('status').cast(pl.String),
                assembly_status = pl.col('assembly_info').struct.field('assembly_status').cast(pl.String)
            ).with_columns(
                strain_name = pl.concat_list([pl.col('strain'), pl.col('isolate'), pl.col('biosample')]).list.drop_nulls().list.unique()
            ).drop(
                ['biosample', 'isolate', 'strain']
            ).sort(
                [pl.col('checkm_completeness'), pl.col('accession')], 
                descending=[True, False]
            ).select("accession", 
                    "assembly_level",
                    "tax_id", 
                    "strain_name", 
                    "checkm_completeness", 
                    "checkm_contamination", 
                    "ncbis_species_match_ani", 
                    "ncbis_species_match_name", 
                    "ncbis_species_match_assembly", 
                    "ncbis_species_match_coverage",
                    "ncbis_best_match_ani", 
                    "ncbis_best_match_name", 
                    "ncbis_best_match_assembly", 
                    "ncbis_best_match_coverage", 
                    "is_atypical", 
                    "paired_assembly_status", 
                    "assembly_status", 
                    "genome_notes", 
                    "type_display_text", 
                    "type_label")
    )
   
    return df


def screen_accessions_for_atypical(ncbi_accession_report_df):
    atypical_accs = ncbi_accession_report_df.filter(
        pl.col("is_atypical") | 
        (pl.col("paired_assembly_status") == "suppressed") |
        (pl.col("assembly_status") == "suppressed")
    )
    typical_accs = ncbi_accession_report_df.filter(
        (pl.col("is_atypical").not_().fill_null(True)) & 
        (pl.col("paired_assembly_status") != "suppressed") &
        (pl.col("assembly_status") != "suppressed")
        )
    return atypical_accs, typical_accs


def screen_accessions_for_bad_types(ncbi_accession_report_df):
    bad_types = ncbi_accession_report_df.filter(
        pl.col("genome_notes").list.contains("not used as type") | 
        (pl.col("type_display_text") == "not used as type") |
        (pl.col("type_label") == "EXCLUDE_TYPE_MATERIAL")
        )
    good_types = ncbi_accession_report_df.filter(
        pl.col("genome_notes").list.contains("not used as type").not_().fill_null(True) &
        (pl.col("type_display_text") != "not used as type") &
        (pl.col("type_label") != "EXCLUDE_TYPE_MATERIAL")
        )
    return bad_types, good_types


def QC_screen_genome_accessions(types_with_genomes, api_key):

    failed_types = []
    passed_types = []
    
    genome_accessions = [str(type_strain.genome_acc.get("accession")) for type_strain in types_with_genomes]
    
    if len(genome_accessions) > 0:
        
        logger.info(f"Screening {len(genome_accessions)} genomes.")
        pbar = tqdm.tqdm(total=len(genome_accessions), desc="Confirming known genome info", unit="type strain", ncols=100, colour="magenta")
        ncbi_accession_report = asyncio.run(request_ncbi_checkm(genome_accessions, api_key, pbar))
        pbar.close()

        ncbi_accession_df = convert_ncbi_accession_report_to_polars_df(ncbi_accession_report).collect(engine="streaming")

        atypical_accs, good_accs = screen_accessions_for_atypical(ncbi_accession_df)
        bad_type_accs, good_accs = screen_accessions_for_bad_types(ncbi_accession_df)

        both_bad_count = len(set(atypical_accs["accession"].to_list()).intersection(set(bad_type_accs["accession"].to_list())))
        atypical_count, bad_type_count = atypical_accs.height - both_bad_count, bad_type_accs.height - both_bad_count
        all_bad_accs = set(atypical_accs["accession"].to_list() + bad_type_accs["accession"].to_list())

        logger.debug(f"Found {atypical_count} atypical, {bad_type_count} bad type, and {both_bad_count} accession(s) that were both.")
        if len(all_bad_accs) > 0:
            logger.debug(f"Removing these accessions and updating type strains accordingly.")

        for type_strain in types_with_genomes:
            if type_strain.genome_acc.get("accession") in all_bad_accs:
                logger.debug(f"Removing genome accession {type_strain.genome_acc.get('accession')} from {type_strain.parent_species} {sorted(type_strain.type_names)[0]}.")
                type_strain.genome_acc = None
                failed_types.append(type_strain)
            else:
                passed_types.append(type_strain)

    else:
        logger.info("No genome accessions to screen. Continuing.")

    return passed_types, failed_types


def sort_by_prefix(items, custom_order):
    custom_order = [c.lower() for c in custom_order]

    def sort_key(item):
        item_lower = item.lower()
        for i, prefix in enumerate(custom_order):
            if item_lower.startswith(prefix):
                return (i, item_lower)
        return (len(custom_order), item_lower)

    return sorted(items, key=sort_key)


def fetch_missing_genomes(missing_genome, api_key):

    updated_genomes = []

    query_terms = []
    strain_names = []
    for type_strain in missing_genome:
        if type_strain.species_ncbi_tax_id is not None:
            if isinstance(type_strain.species_ncbi_tax_id, list):
                query_terms.extend([str(tax_id) for tax_id in type_strain.species_ncbi_tax_id])
            else:
                query_terms.append(str(type_strain.species_ncbi_tax_id))
        if type_strain.type_names is not None:
            if isinstance(type_strain.type_names, list):
                strain_names.extend([variant for name in type_strain.type_names for variant in (name, re.sub(f'{type_strain.parent_genus} sp. ', '', name, flags=re.IGNORECASE), re.sub("strain ", "", name, flags=re.IGNORECASE), name.replace(" ", ""), name.replace(" ", "-"), name.replace(" ", "_"), re.sub(r'[^a-zA-Z0-9\s]', '', name), re.sub(r'[^a-zA-Z0-9]', '', name))])
            else:
                strain_names.extend([type_strain.type_names, re.sub(f'{type_strain.parent_genus} sp. ', '', type_strain.type_names, flags=re.IGNORECASE), re.sub("strain ", "", type_strain.type_names, flags=re.IGNORECASE), type_strain.type_names.replace(" ", ""), type_strain.type_names.replace(" ", "-"), type_strain.type_names.replace(" ", "_"), re.sub(r'[^a-zA-Z0-9\s]', '', type_strain.type_names), re.sub(r'[^a-zA-Z0-9]', '', type_strain.type_names)])
    strain_names = set(strain_names)

    pbar = tqdm.tqdm(total=len(query_terms), desc="Retrieving genome info", unit="type strain", ncols=100, colour="magenta")
    data_list = asyncio.run(request_ncbi_genomes(query_terms, api_key, pbar, strain_names))
    pbar.close()
    logger.debug(f"Retrieved {len(data_list)} genome records from GenBank.")

    if len(data_list) == 0:
        logger.info("No genome records retrieved from GenBank. Continuing")
        return missing_genome

    logger.info("Filtering genome records for atypical, suppressed, and bad type material.")

    ncbi_data_df = convert_ncbi_accession_report_to_polars_df(data_list).filter(
        (pl.col("is_atypical").not_().fill_null(True)) &
        (pl.col("genome_notes").list.contains("not used as type").not_().fill_null(True)) &
        (pl.col("type_display_text") != "not used as type") &
        (pl.col("type_label") != "EXCLUDE_TYPE_MATERIAL") &
        (pl.col("paired_assembly_status") != "suppressed") &
        (pl.col("assembly_status") != "suppressed")
    ).with_columns(
        corrected_strain_name=pl.col('strain_name').list.eval(pl.element().str.replace_all(r"[^a-zA-Z0-9]", ""))
    ).collect(engine="streaming")

    for type_strain in tqdm.tqdm(missing_genome, desc="Processing missing genome info", unit="type strain", ncols=100, colour="magenta"):
        filtered = ncbi_data_df.filter(pl.col('tax_id').is_in(type_strain.species_ncbi_tax_id))
        if filtered.height > 0:
            try:
                if type(type_strain.type_names) == str:
                    type_strain.type_names = [ts.strip() for ts in type_strain.type_names.split(",")]
                candidate_names = set(
                    type_strain.type_names
                    + [re.sub(r'[^a-zA-Z0-9]', "", f'{type_strain.parent_genus} sp. {x}') for x in type_strain.type_names]
                    + [re.sub("strain ", "", x, flags=re.IGNORECASE) for x in type_strain.type_names]
                    + [x.replace(" ", "") for x in type_strain.type_names]
                    + [x.replace(" ", "-") for x in type_strain.type_names]
                    + [x.replace(" ", "_") for x in type_strain.type_names]
                    + [re.sub(r'[^a-zA-Z0-9\s]', '', x) for x in type_strain.type_names]
                    + [re.sub(r'[^a-zA-Z0-9]', '', x) for x in type_strain.type_names]
                )
                
                strain_filtered = filtered.filter(
                    pl.col('corrected_strain_name').list.eval(pl.element().is_in(candidate_names, nulls_equal=True)).list.any()
                )
            except Exception as e:
                logger.error(f"Error filtering strains for {type_strain}: {e}")
                strain_filtered = filtered.clear()
            if strain_filtered.height > 0:
                best_hit = strain_filtered.row(0, named=True)
                best_hit_strain = best_hit.get("strain_name", [None])[0]

                # sort type_names so dsm, atcc, nctc come first, and the best hit strain name is first if it is in the type names, then alphabetically
                if type_strain.type_names is not None:
                    custom_order = [best_hit_strain.lower(), "dsm", "atcc", "nctc"]
                    type_strain.type_names = sort_by_prefix(type_strain.type_names, custom_order)

                type_strain.genome_acc = {"accession": best_hit['accession'].split('.')[0],
                                          "assembly level": best_hit['assembly_level'],
                                          "checkm_completeness": best_hit['checkm_completeness'],
                                          "checkm_contamination": best_hit['checkm_contamination'] if best_hit['checkm_contamination'] is not None else 0 if best_hit['checkm_completeness'] is not None else '',
                                          "ncbis_species_match_ani": best_hit['ncbis_species_match_ani'],
                                          "ncbis_species_match_name": best_hit['ncbis_species_match_name'],
                                          "ncbis_species_match_assembly": best_hit['ncbis_species_match_assembly'],
                                          "ncbis_species_match_coverage": best_hit['ncbis_species_match_coverage'],
                                          "ncbis_best_match_ani": best_hit['ncbis_best_match_ani'],
                                          "ncbis_best_match_name": best_hit['ncbis_best_match_name'],
                                          "ncbis_best_match_assembly": best_hit['ncbis_best_match_assembly'],
                                          "ncbis_best_match_coverage": best_hit['ncbis_best_match_coverage']
                                          }
            else:
                logger.debug(f"No genome data found for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species} matching type strain names. Just FYI")
        else:
            logger.debug(f"No genome data found for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species}. Just FYI")
        updated_genomes.append(type_strain)

    return updated_genomes


def retrieve_missing_genome_info(lpsn_types, api_key):
    
    logger.info("Retrieving genome information from GenBank.")

    query_terms = []
    has_genome, missing_genome = get_haves_and_have_nots(lpsn_types, "genome_acc")      
    
    passed_qc, failed_qc = QC_screen_genome_accessions(has_genome, api_key)

    missing_genome = missing_genome + failed_qc
    if len(missing_genome) == 0:
        logger.info("No missing genome information to retrieve from GenBank. Continuing")
        return passed_qc
    
    updated_genomes = fetch_missing_genomes(missing_genome, api_key)

    return passed_qc + updated_genomes


def retrieve_genome_sequences(lpsn_types, output_path, threads, api_key):

    logger.info("Retrieving genome sequences from GenBank.")

    has_genome_seq, missing_genome_seq = get_haves_and_have_nots(lpsn_types, "genome_acc")
    logger.info(f"{len(missing_genome_seq)} type strains have no known genome sequence. Downloading sequences for remaining {len(has_genome_seq)} type strains.")
    with open(output_path / "genome_accessions.txt", "w") as f:
        accessions = {
            x.genome_acc["accession"]
            for x in has_genome_seq
            if x.genome_acc is not None
        }
        if len(accessions) == 0:
            logger.info("No genome accessions found for any type strains. Skipping genome sequence retrieval.")
            return
        f.write("\n".join(accessions))

    make_dir(output_path / "sequences" / "genomes")
    command = f'datasets download genome accession --inputfile {output_path}/genome_accessions.txt --api-key {api_key} --filename {output_path}/sequences/genomes/genomes.zip --dehydrated --include genome,gbff'
    logger.info("Downloading dehydrated genome sequences from GenBank")    
    try:
        info = subprocess.Popen(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            output = info.stderr.readline()
            if not output:
                break
            logger.info(output.strip())
        logger.info("Finished downloading genome sequences.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error downloading genome sequences: {e.stderr}")
        sys.exit(1)


def retrieve_checkm_info(lpsn_types, api_key):
    
    logger.info("Retrieving checkm information from GenBank.")

    query_terms = []
    
    has_genome, missing_genome = get_haves_and_have_nots(lpsn_types, "genome_acc")        
    query_terms = [str(type_strain.genome_acc.get("accession")) for type_strain in has_genome if type_strain.genome_acc is not None and type_strain.genome_acc.get("checkm_completeness") is None]
    if len(has_genome) == 0:
        logger.info("No checkm information to retrieve from GenBank. Continuing")
        return has_genome + missing_genome
    pbar = tqdm.tqdm(total=len(query_terms), desc="Retrieving checkm info", unit="type strain", ncols=100, colour="magenta")
    data_list = asyncio.run(request_ncbi_checkm(query_terms, api_key, pbar))
    pbar.close()
    
    logger.debug(f"Retrieved {len(data_list)} checkm records from GenBank.")
    if len(data_list) == 0:
        logger.info("No checkm information retrieved from GenBank. Continuing")
        return has_genome + missing_genome
        
    logger.info("Processing checkm information.")

    df = convert_ncbi_accession_report_to_polars_df(data_list).collect(engine="streaming")
       
    for type_strain in tqdm.tqdm(has_genome, desc="Processing checkm info", unit="type strain", ncols=100, colour="magenta"):
        filtered = df.filter(pl.col('accession').str.replace('\\.[0-9]$', '') == type_strain.genome_acc.get("accession")).fill_null('')
        if len(filtered) > 0:
            best_hit = filtered.rows(named=True)[0]
            best_hit_strain = best_hit.get("strain_name", [None])[0]
            best_hit_completeness = best_hit.get("checkm_completeness", '') if best_hit.get("checkm_completeness") is not None else ''
            best_hit_contamination = best_hit.get("checkm_contamination", '') if best_hit.get("checkm_contamination") is not None else 0 if best_hit.get("checkm_completeness") is not None else ''
            best_hit_ncbi_species_match_ani = best_hit.get("ncbis_species_match_ani", '') if best_hit.get("ncbis_species_match_ani") is not None else ''
            best_hit_ncbi_species_match_name = best_hit.get("ncbis_species_match_name", '') if best_hit.get("ncbis_species_match_name") is not None else ''
            best_hit_ncbi_species_match_assembly = best_hit.get("ncbis_species_match_assembly", '') if best_hit.get("ncbis_species_match_assembly") is not None else ''
            best_hit_ncbi_species_match_assembly_coverage = best_hit.get("ncbis_species_match_coverage", '') if best_hit.get("ncbis_species_match_coverage") is not None else ''
            best_hit_ncbi_best_match_ani = best_hit.get("ncbis_best_match_ani", '') if best_hit.get("ncbis_best_match_ani") is not None else ''
            best_hit_ncbi_best_match_name = best_hit.get("ncbis_best_match_name", '') if best_hit.get("ncbis_best_match_name") is not None else ''
            best_hit_ncbi_best_match_assembly = best_hit.get("ncbis_best_match_assembly", '') if best_hit.get("ncbis_best_match_assembly") is not None else ''
            best_hit_ncbi_best_match_assembly_coverage = best_hit.get("ncbis_best_match_coverage", '') if best_hit.get("ncbis_best_match_coverage") is not None else ''

            if type_strain.type_names is not None:
                    custom_order = [best_hit_strain.lower(), "dsm", "atcc", "nctc"]
                    type_strain.type_names = sort_by_prefix(type_strain.type_names, custom_order)
            type_strain.genome_acc.update({"checkm_completeness": best_hit_completeness, 
                                           "checkm_contamination": best_hit_contamination, 
                                           "ncbis_species_match_ani": best_hit_ncbi_species_match_ani, 
                                           "ncbis_species_match_name": best_hit_ncbi_species_match_name, 
                                           "ncbis_species_match_assembly": best_hit_ncbi_species_match_assembly, 
                                           "ncbis_species_match_coverage": best_hit_ncbi_species_match_assembly_coverage,
                                           "ncbis_best_match_ani": best_hit_ncbi_best_match_ani, 
                                           "ncbis_best_match_name": best_hit_ncbi_best_match_name, 
                                           "ncbis_best_match_assembly": best_hit_ncbi_best_match_assembly,
                                           "ncbis_best_match_coverage": best_hit_ncbi_best_match_assembly_coverage})
    return has_genome + missing_genome   


def rehydrate_genome_sequences(output_path, threads, input_taxon):

    logger.info("Rehydrating downloaded genome sequences.")
    try:
        unzip_file(output_path / "sequences" / "genomes" / "genomes.zip", output_path / "sequences" / "genomes")
    except FileNotFoundError:
        logger.warning("No genomes.zip file found to unzip. Skipping genome rehydration.")
        return

    command = f'datasets rehydrate --directory {output_path / "sequences" / "genomes"} --max-workers {threads}'

    try:
        info = subprocess.Popen(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            output = info.stdout.readline()
            if not output:
                break
            logger.info(output.strip())
    except subprocess.CalledProcessError as e:
        logger.error(f"Error rehydrating genome sequences: {e.stderr}")
        sys.exit(1)


def check_16S_retrieval(output_path, input_taxon, accessions):

    retrieved_seqs = []    
    for record in SeqIO.parse( output_path / "sequences" / "16S" / "16S.fasta", "fasta"):
        retrieved_seqs.append(record.id.split(".")[0])
    
    missing = set(accessions) - set(retrieved_seqs)
    if len(missing) > 0:
        logger.warning(f"Could not retrieve {len(missing)} 16S rRNA gene sequences.")
        logger.debug(f"Missing accessions: {', '.join(missing)}")


async def validate_accessions(accessions, email, api_key, pbar=None, rate_limiter=None):
    term = " OR ".join(f"{x}[Accession]" for x in accessions)
    valid_subset = set()
    max_retries = 5
    for attempt in range(max_retries):
        try:
            # esearch async call
            search_result_xml = await esearch(term, db="nuccore", retmax=10**4, api_key=api_key, email=email, rate_limiter=rate_limiter)
            search_result = Entrez.read(search_result_xml)
            webenv = search_result["WebEnv"]
            query_key = search_result["QueryKey"]
            break
        except Exception as e:
            if attempt < max_retries - 1:
                await asyncio.sleep(2 ** attempt)  # Exponential backoff
                logger.warning(f"Attempt {attempt + 1} failed to retrieve search results. Retrying...")
                logger.debug(f"Error details: {e}")
                logger.debug(f"Search result XML: {search_result_xml.getvalue().decode() if 'search_result_xml' in locals() else 'No search result XML available'}")
            else:
                logger.error(f"Failed to retrieve search results after {max_retries} attempts.")
                exit(1)
                

    if search_result.get("Count")!=0 : 
        for attempt in range(max_retries):
            try:
                # esummary async call
                records = await esummary(query_key, webenv, db="nuccore", api_key=api_key, email=email, rate_limiter=rate_limiter)
                records = Entrez.read(records)
                break
            except Exception as e:
                if attempt < max_retries - 1:
                    await asyncio.sleep(2 ** attempt)  # Exponential backoff
                else:
                    logger.error(f"Failed to retrieve esummary results after {max_retries} attempts.")
                    exit(1)
        # Filter valid records
        for x in records:
            length = x.get("Length")
            status = x.get("Status")
            acc_ver = x.get("AccessionVersion")
            if (
                length is not None
                and 1000 <= int(length) <= 2500
                and status == "live"
                and acc_ver
                and len(acc_ver) < 12
            ):
                valid_subset.add(acc_ver.split(".")[0])
        if pbar:
            pbar.update(len(accessions))
    return valid_subset


async def grouped_validate_accessions(accessions, max_terms, email, api_key, max_concurrency=9, pbar=None):
    semaphore = asyncio.Semaphore(max_concurrency)
    batches = [accessions[i:i + max_terms] for i in range(0, len(accessions), max_terms)]
    rate_limiter = RateLimiter(max_rate=8)

    async def limited_validate(batch, pbar=pbar, rate_limiter=None):
        async with semaphore:
            return await validate_accessions(batch, email, api_key, pbar=pbar, rate_limiter=rate_limiter)
    tasks = [limited_validate(batch, pbar=pbar, rate_limiter=rate_limiter) for batch in batches]
    batch_results = await asyncio.gather(*tasks)
    
    # Merge all sets
    valid_accessions = set().union(*batch_results)
    return valid_accessions


def get_acc_seq_lengths(acc_list, api_key=None, email=None):
    accessions = [x.split(".")[0] for x in acc_list]
    if not accessions:
        return []
    valid_accessions = set()
    pbar = tqdm.tqdm(total=len(accessions), desc="Validating 16S rRNA gene sequences", unit="type strain", ncols=100, colour="magenta")
    valid_accessions = asyncio.run(grouped_validate_accessions(accessions, max_terms=100, email=email, api_key=api_key, pbar=pbar))
    pbar.close()

    return [x for x in accessions if x in valid_accessions]


async def retrieve_accessions(accessions, email, api_key, rate_limiter=None):
    term = " OR ".join(f"{x}[Accession]" for x in accessions)
    max_retries = 5
    for attempt in range(max_retries):
        try:
            # esearch async call
            search_result_xml = await esearch(term, db="nuccore", retmax=10**4, api_key=api_key, email=email, rate_limiter=rate_limiter)
            search_result = Entrez.read(search_result_xml)
            break
        except Exception as e:
            if attempt < max_retries - 1:
                await asyncio.sleep(2 ** attempt)  # Exponential backoff
            else:
                logger.error(f"Failed to retrieve search results after {max_retries} attempts.")
                raise e

    return set(search_result.get("IdList", []))


async def grouped_retrieve_16S(accessions, max_terms, email, api_key, max_concurrency=9):
    semaphore = asyncio.Semaphore(max_concurrency)
    batches = [accessions[i:i + max_terms] for i in range(0, len(accessions), max_terms)]
    rate_limiter = RateLimiter(max_rate=9)

    async def limited_retrieve(batch, rate_limiter=None):
        async with semaphore:
            return await retrieve_accessions(batch, email, api_key, rate_limiter=rate_limiter)
    tasks = [limited_retrieve(batch, rate_limiter=rate_limiter) for batch in batches]
    batch_results = await asyncio.gather(*tasks)
    gi_list = set().union(*batch_results)
    return gi_list

 
def retrieve_16S_sequences(lpsn_types, output_path, email, input_taxon, api_key):
    
    logger.info("Retrieving 16S rRNA gene sequences from GenBank.")

    retmax = 10**6

    accessions = [x.rRNA_acc.split(".")[0] for x in lpsn_types if x.rRNA_acc is not None]
    giList = asyncio.run(grouped_retrieve_16S(accessions, max_terms=100, email=email, api_key=api_key))
    if len(giList) == 0:
        logger.warning(f"No 16S rRNA gene sequences found. Continuing without 16S sequences.")
        search_results = None
    
    if len(giList) > 0:
        logger.info(f"Found {len(giList)} 16S rRNA gene sequences to retrieve from GenBank.")
        try:
            search_handle = Entrez.epost(db="nucleotide", id=",".join(giList), email=email, api_key=api_key)
            search_results = Entrez.read(search_handle)
            webenv, query_key = search_results["WebEnv"], search_results["QueryKey"] 
        except Exception as e:
            logger.error(f"Error posting search results: {e}")
            sys.exit(1)

    make_dir( output_path / "sequences" / "16S" )

    if search_results is not None:
        with open( output_path / "sequences" / "16S" / "16S.fasta", "w" ) as f:
            for start in tqdm.tqdm(range(0, len(giList), 100), desc="Downloading 16S sequences", unit="batch", ncols=100, colour="magenta"):
                handle = Entrez.efetch(db='nucleotide', rettype="fasta", retmode='text', retstart = start, retmax=100, webenv= webenv, query_key= query_key, email=email, api_key=api_key)
                data = handle.read().replace("\n\n", "\n")
                f.write(data)
        check_16S_retrieval(output_path, input_taxon, accessions)
    else:
        open(output_path / "sequences" / "16S" / "16S.fasta", "a").close()



def retrieve_sequences_workflow(lpsn_types, output_path, threads, dev_mode, input_taxon, sequence_download, email=None, api_key=None, mmseqs_16s_check="all"):
    
    logger.info("Step 4 of 4: Generating ouputs")
    
    if sequence_download == "all":
        retrieve_genome_sequences(lpsn_types, output_path, threads, api_key)
        rehydrate_genome_sequences(output_path, threads, input_taxon)
        retrieve_16S_sequences(lpsn_types, output_path, email, input_taxon, api_key)
        lpsn_types = confirm_genbank_16S_hits_with_mmseqs(lpsn_types, output_path, threads, mmseqs_16s_check)
        tidy_16S_dir(output_path)
    elif sequence_download == "16S":
        logger.info("Skipping genome sequence retrieval as per user request.\n")
        retrieve_16S_sequences(lpsn_types, output_path, email, input_taxon, api_key)
        lpsn_types = confirm_genbank_16S_hits_with_mmseqs(lpsn_types, output_path, threads, mmseqs_16s_check)
        tidy_16S_dir(output_path)
    elif sequence_download == "genomes":
        logger.info("Skipping 16S sequence retrieval as per user request.\n")
        retrieve_genome_sequences(lpsn_types, output_path, threads, api_key)
        rehydrate_genome_sequences(output_path, threads, input_taxon)
    elif sequence_download == "none":
        logger.info("Skipping sequence download steps as per user request.\n")
    else:
        logger.error(f"Invalid sequence_download option: {sequence_download}. Must be one of 'all', '16S', 'genomes', or 'none'.")
        sys.exit(1)

    output_metadata(lpsn_types, output_path)

    tidy_genome_dir(output_path)

    logger.info("Sequence retrieval workflow complete.\n")

    logger.info("###################################")
    logger.info("###     Ratatoskr finished!     ###")
    logger.info("###################################\n")


def retrieve_info_from_genbank(lpsn_types, output_path, threads, dev_mode, input_taxon, sequence_download, email=None, api_key=None, mmseqs_16s_check="all"):
    """
    Retrieve sequences from GenBank for the given LPSN type strains.
    """
    logger.info("Step 3 of 4: Retrieving information from GenBank.")
    lpsn_types = retrieve_ncbi_taxon_ids(lpsn_types, api_key)

    lpsn_types = retrieve_missing_16S_info(lpsn_types, email, api_key, output_path)
    lpsn_types = retrieve_missing_genome_info(lpsn_types, api_key)
    lpsn_types = retrieve_checkm_info(lpsn_types, api_key)
    for type_strain in lpsn_types:
        if type_strain.genome_acc is not None:
            type_strain.genome_acc['accession'] = type_strain.genome_acc['accession'].split(".")[0]

    logger.success("Metadata retrieval from GenBank complete.\n")

    retrieve_sequences_workflow(lpsn_types, output_path, threads, dev_mode, input_taxon, sequence_download, email=email, api_key=api_key, mmseqs_16s_check=mmseqs_16s_check)


