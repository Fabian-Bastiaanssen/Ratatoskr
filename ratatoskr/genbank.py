from datetime import datetime
import os
import subprocess
import sys
import asyncio
import aiohttp
import time
from io import BytesIO

from Bio import Entrez, SeqIO
from loguru import logger
import tqdm
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl

from ratatoskr.misc import get_haves_and_have_nots, fetch_data, tidy_genome_dir
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

    return has_ncbi_taxid +missing_ncbi_taxid
async def search_16s(batch, email, api_key, pbar, rate_limiter):
    max_retries = 5
    batch_df = pl.LazyFrame([{"ncbi_tax_id": type_strain.strain_ncbi_tax_id +  type_strain.species_ncbi_tax_id, "parent_species_id": type_strain.parent_species_id, "parent_subspecies_id": type_strain.parent_subspecies_id} for type_strain in batch], schema={"ncbi_tax_id": pl.List(pl.String), "parent_species_id": pl.Int64, "parent_subspecies_id": pl.Int64})
    batch_df = batch_df.explode(["ncbi_tax_id"])
    ids = batch_df.select(pl.col("ncbi_tax_id")).unique().collect().to_series().drop_nulls().to_list()
    have_rRNA = []
    missing_rRNA = []
    # ids =[item for v in [type_strain.strain_ncbi_tax_id,type_strain.species_ncbi_tax_id] for item in (v if isinstance(v, list) else [v])]
    term = " OR ".join([f"txid{i}[Organism:exp]" for i in ids if i is not None]) + " AND 16S[Title]"
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
        records_df = pl.LazyFrame({'Length': int(x.get("Length")), 'Status': x.get("Status"), 'AccessionVersion': x.get("AccessionVersion"), 'CreateDate': x.get("CreateDate"), 'TaxId': str(int(x.get("TaxId")))} for x in records)
        records_df = records_df.filter((pl.col("Length") >= 750) & (pl.col("Length") <= 2500) & (pl.col("Status") == "live") & (pl.col("AccessionVersion").str.len_chars() < 12))
        batch_df = batch_df.join(records_df.select(pl.col("TaxId").alias("ncbi_tax_id"), pl.col("AccessionVersion"), pl.col("CreateDate")), on="ncbi_tax_id", how="left")
        batch_df = batch_df.sort("CreateDate", descending=False).group_by("parent_species_id", "parent_subspecies_id").agg(pl.col("AccessionVersion").first().alias("rRNA_acc")).collect()
        subspecies_lookup = (batch_df.filter(pl.col("parent_subspecies_id").is_not_null())
            .select("parent_subspecies_id", "rRNA_acc")
            .to_dict(as_series=False))
        subspecies_map = dict(zip(subspecies_lookup["parent_subspecies_id"],
                            subspecies_lookup["rRNA_acc"]))

        species_lookup = (batch_df.filter(pl.col("parent_subspecies_id").is_null(), pl.col("parent_species_id").is_not_null())
            .select("parent_species_id", "rRNA_acc")
            .to_dict(as_series=False))
        species_map = dict(zip(species_lookup["parent_species_id"], species_lookup["rRNA_acc"]))

        for type_strain in batch:
            if type_strain.parent_subspecies_id is not None:
                hit = subspecies_map.get(type_strain.parent_subspecies_id)
            else:
                hit = species_map.get(type_strain.parent_species_id)

            if hit is not None:
                type_strain.rRNA_acc = hit
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
                    type_strain.rRNA_acc = record
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

async def retrieve_missing_16S(missing_rRNA, max_terms, email, api_key, max_concurrency=10, pbar=None):
    semaphore = asyncio.Semaphore(max_concurrency)
    rate_limiter = RateLimiter(max_rate=8)
    batches = batch_missing_rrna(missing_rRNA, max_ids=max_terms)
    async def limited_search(batch, pbar, rate_limiter):
        async with semaphore:
            return await search_16s(batch, email, api_key, pbar=pbar, rate_limiter=rate_limiter)
    tasks = [limited_search(batch, pbar=pbar, rate_limiter=rate_limiter) for batch in batches]
    batch_results = await asyncio.gather(*tasks)
    hits = [x for r in batch_results for x in r["hits"]]
    misses = [x for r in batch_results for x in r["misses"]]
    if len(misses) > 0:
        logger.debug(f"{len(misses)} type strains missing 16S rRNA gene information after initial search. Attempting binomial and strain name search for these strains.")
        async def limited_binomial_search(type_strain, pbar, rate_limiter):
            async with semaphore:
                return await search_16s_binomial(type_strain, email, api_key, pbar=pbar, rate_limiter=rate_limiter)
        binomial_tasks = [limited_binomial_search(type_strain, pbar=pbar, rate_limiter=rate_limiter) for type_strain in misses]
        binomial_results = await asyncio.gather(*binomial_tasks)
        return hits + binomial_results
    return hits            
def retrieve_missing_16S_info(lpsn_types, email, api_key):
    """
    Retrieve missing 16S rRNA gene information from GenBank for the given LPSN type strains.
    """
    logger.info("Retrieving missing 16S rRNA gene information from GenBank.")
    Entrez.api_key = api_key
    incorrect_types = [ts for ts in lpsn_types if type(ts) == str]
    has_rRNA, missing_rRNA = get_haves_and_have_nots(lpsn_types, "rRNA_acc")
    for type_strain in missing_rRNA:
        # check if type_strain is missing attribute species_id
        if not hasattr(type_strain, "parent_species_id"):
            logger.error(f"Type strain {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species} is missing parent_species_id attribute. This is required for 16S retrieval. Skipping this type strain for 16S retrieval.")
            logger.debug(type_strain)
    pbar = tqdm.tqdm(total=len(missing_rRNA), desc="Retrieving missing 16S rRNA gene info", unit="type strain", ncols=100, colour="magenta")
    missing_rRNA = asyncio.run(retrieve_missing_16S(missing_rRNA, max_terms=99, email=email, api_key=api_key, max_concurrency=8, pbar=pbar))
    pbar.close()

    return has_rRNA + missing_rRNA


async def request_ncbi_genomes(query_terms, api_key, pbar):
    headers = {
        "accept": "application/json",
        "api-key": api_key
    }
    logger.debug(f"Requesting genome information for {len(query_terms)} taxon IDs from GenBank.")
    session = aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60000),connector=aiohttp.TCPConnector(limit=9))
    sem = asyncio.Semaphore(9)
    urls = [f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{"%2C".join(subset).replace(" ", "%20")}/dataset_report?returned_content=COMPLETE&page_size=1000" for subset in [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]]
    query_lengths = [len(subset) for subset in [query_terms[i:i + 100] for i in range(0, len(query_terms), 100)]]
    parrallel_tasks = [fetch_data(url, session, headers, sem, 'reports', pbar, query_length) for url, query_length in zip(urls, query_lengths)]
    async with sem:
        data_list = await asyncio.gather(*parrallel_tasks)
    await session.close()
    return [item for sublist in data_list for item in sublist if sublist is not None]
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

def retrieve_missing_genome_info(lpsn_types, api_key):
    
    logger.info("Retrieving genome information from GenBank.")

    query_terms = []
    has_genome_first, missing_genome = get_haves_and_have_nots(lpsn_types, "genome_acc")      
    query_terms_has_genome = [str(type_strain.genome_acc.get("accession")) for type_strain in has_genome_first if type_strain.genome_acc is not None]
    has_genome =[]
    if len(query_terms_has_genome) !=0:
        pbar = tqdm.tqdm(total=len(query_terms_has_genome), desc="Confirming known genome info", unit="type strain", ncols=100, colour="magenta")
        data_list_has = asyncio.run(request_ncbi_checkm(query_terms_has_genome, api_key, pbar))
        pbar.close()
        df = pl.DataFrame(data_list_has).with_columns(pl.col('accession').str.replace('\\.[0-9]+$', ''))
        for type_strain in has_genome_first:
            if 'atypical' not in df.schema.get('assembly_info', pl.Struct({'none': pl.String})).to_schema():
                df = df.with_columns(assembly_info={'atypical': {'is_atypical': None, 'warnings': [None]}})
            filtered_df = df.filter(pl.col('accession') == str(type_strain.genome_acc.get("accession"))).with_columns(pl.col('assembly_info').struct.field('atypical')).unnest('atypical')
            if len(filtered_df) == 0:
                logger.info(f'{type_strain.parent_species} was assigned accesion {type_strain.genome_acc["accession"]}, but this record is supressed or missing. Redoing assignment')
                type_strain.genome_acc['accession'] =''
                missing_genome.append(type_strain)
            elif filtered_df.select(pl.col('is_atypical').any()).item():
                if filtered_df.select(pl.col('warnings')).item() is not [None]:
                    reasons = "reasons: " + ", ".join(filtered_df.select(pl.col('warnings')).item())
                else:
                    reasons = "unknown reasons"
                logger.info(f'{type_strain.parent_species} was assigned accesion {type_strain.genome_acc["accession"]}, but this record is flagged as atypical in GenBank for {reasons}. Redoing assignment')
                type_strain.genome_acc['accession'] =''
                missing_genome.append(type_strain)
            else:
                has_genome.append(type_strain)
                
        
    query_terms = [str(ncbi_id) for type_strain in missing_genome for ncbi_id in type_strain.species_ncbi_tax_id if type_strain.species_ncbi_tax_id is not None]
    if len(missing_genome) == 0:
        logger.info("No missing genome information to retrieve from GenBank. Continuing")
        return has_genome + missing_genome
    pbar = tqdm.tqdm(total=len(query_terms), desc="Retrieving genome info", unit="type strain", ncols=100, colour="magenta")
    data_list = asyncio.run(request_ncbi_genomes(query_terms, api_key, pbar))
    pbar.close()
    
    logger.debug(f"Retrieved {len(data_list)} genome records from GenBank.")
    if len(data_list) == 0:
        logger.info("No genome information retrieved from GenBank. Continuing")
        return has_genome + missing_genome
    
    logger.info("Processing genome information.")

    df = pl.LazyFrame(data_list)
    if 'checkm_info' not in df.collect_schema().names():
        df = df.with_columns(checkm_info={'checkm_species_tax_id': None, "completeness": None, "contamination": None})

    df= df.select([pl.col('accession'), pl.col('organism').struct.field("*"), pl.col('assembly_info').struct.field(['assembly_level','biosample']), pl.col('checkm_info').struct.field('checkm_species_tax_id').alias('checkm_tax_id'), pl.col('checkm_info').struct.field('completeness').alias('checkm_completeness'), pl.col('checkm_info').struct.field('contamination').alias('checkm_contamination')])
    if 'strain' not in df.collect_schema().get('infraspecific_names').to_schema():
        df = df.with_columns(pl.col('infraspecific_names').struct.with_fields(strain=pl.lit(None)))
    if 'strain' not in df.collect_schema().get('biosample').to_schema():
        df = df.with_columns(pl.col('biosample').struct.with_fields(strain=pl.lit(None)))

    df = df.with_columns(pl.col('infraspecific_names').struct.field("strain").alias('infraspecific_names'), pl.col('biosample').struct.field('strain').alias('biosample'))
    enum_type = pl.Enum(("Complete Genome", "Chromosome", "Scaffold", "Contig"))
    df = df.with_columns(strain=pl.concat_list(pl.col('biosample'), pl.col('infraspecific_names')).list.unique()).unique().sort(pl.col('assembly_level').cast(enum_type), pl.col('checkm_completeness'), descending=[False, True], nulls_last=True).collect().group_by('tax_id').all()
    for type_strain in tqdm.tqdm(missing_genome, desc="Processing missing genome info", unit="type strain", ncols=100, colour="magenta"):
        filtered = df.filter(pl.col('tax_id').is_in(type_strain.species_ncbi_tax_id)| pl.col("checkm_tax_id").list.set_intersection(type_strain.species_ncbi_tax_id).list.len() > 0)
        if len(filtered) > 0:
            try:
                if type(type_strain.type_names) == str:
                    type_strain.type_names = [ts.strip() for ts in type_strain.type_names.split(",")]
                strain_filtered = filtered.explode(['accession', 'assembly_level', 'strain', 'checkm_completeness', 'checkm_contamination']).filter(pl.col('strain').list.eval(pl.element().is_in(type_strain.type_names, nulls_equal=True)).list.any())
            except Exception as e:
                logger.error(f"Error filtering strains for {type_strain}: {e}")
            if len(strain_filtered) > 0:
                best_hit = strain_filtered.rows(named=True)[0]
                best_hit_strain = best_hit.get("strain", [None])[0]
                # sort type_names so dsm, atcc, ntcc come first, and the best hit strain name is first if it is in the type names, then alphabetically
                if type_strain.type_names is not None:
                    custom_order = [best_hit_strain, "dsm", "atcc", "ntcc"]
                    type_strain.type_names =  sorted(type_strain.type_names, key=lambda x: (x not in custom_order, x))
                type_strain.genome_acc = {"accession": best_hit['accession'].split('.')[0], "assembly level": best_hit['assembly_level'], "checkm_completeness": best_hit['checkm_completeness'], "checkm_contamination": best_hit['checkm_contamination']}
        #     else:
        #         logger.debug(f"No genome data found for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species} matching type strain names. Just FYI")
        # else:
        #     logger.debug(f"No genome data found for {type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species}. Just FYI")
        

    return has_genome + missing_genome


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

    df = pl.LazyFrame(data_list)
    if 'checkm_info' not in df.collect_schema().names():
        df = df.with_columns(checkm_info={"completeness": '', "contamination": ''})
    df= df.select([pl.col('accession'), pl.col('organism').struct.field("*"), pl.col('assembly_info').struct.field(['assembly_level','biosample']),  pl.col('checkm_info').struct.field('completeness').alias('checkm_completeness'), pl.col('checkm_info').struct.field('contamination').alias('checkm_contamination')])

    if 'infraspecific_names' not in df.collect_schema().names():
        df = df.with_columns(infraspecific_names={"strain": pl.lit(None)})
    if 'strain' not in df.collect_schema().get('infraspecific_names').to_schema():
        df = df.with_columns(pl.col('infraspecific_names').struct.with_fields(strain=pl.lit(None)))
    if 'biosample' not in df.collect_schema().names():
        df = df.with_columns(biosample={"strain": pl.lit(None)})
    if 'strain' not in df.collect_schema().get('biosample').to_schema():
        df = df.with_columns(pl.col('biosample').struct.with_fields(strain=pl.lit(None)))

    df = df.with_columns(pl.col('infraspecific_names').struct.field("strain").alias('infraspecific_names'), pl.col('biosample').struct.field('strain').alias('biosample'))
    df = df.with_columns(strain=pl.concat_list(pl.col('biosample'), pl.col('infraspecific_names')).list.unique()).unique()

    df= df.select([pl.col('accession'), pl.col('checkm_completeness'), pl.col('checkm_contamination'), pl.col('strain')])
    for type_strain in tqdm.tqdm(has_genome, desc="Processing checkm info", unit="type strain", ncols=100, colour="magenta"):
        filtered = df.filter(pl.col('accession').str.replace('\\.[0-9]$', '') == type_strain.genome_acc.get("accession")).fill_null('').collect()
        if len(filtered) > 0:
            best_hit = filtered.rows(named=True)[0]
            best_hit_strain = best_hit.get("strain", [None])[0]
            best_hit_completeness = best_hit.get("checkm_completeness", '') if best_hit.get("checkm_completeness") is not None else ''
            best_hit_contamination = best_hit.get("checkm_contamination", '') if best_hit.get("checkm_contamination") is not None else ''
            # sort type_names so dsm, atcc, ntcc come first, and the best hit strain name is first if it is in the type names, then alphabetically
            if type_strain.type_names is not None:
                custom_order = [best_hit_strain, "dsm", "atcc", "ntcc"]
                type_strain.type_names =  sorted(type_strain.type_names, key=lambda x: (x not in custom_order, x))
            
            type_strain.genome_acc.update({"checkm_completeness": best_hit_completeness, "checkm_contamination": best_hit_contamination})
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


def retrieve_info_from_genbank(lpsn_types, output_path, threads, dev_mode, input_taxon, skip_download, email=None, api_key=None):
    """
    Retrieve sequences from GenBank for the given LPSN type strains.
    """
    logger.info("Step 3 of 4: Retrieving information from GenBank.")
    lpsn_types = retrieve_ncbi_taxon_ids(lpsn_types, api_key)

    lpsn_types = retrieve_missing_16S_info(lpsn_types, email, api_key)
    lpsn_types = retrieve_missing_genome_info(lpsn_types, api_key)
    lpsn_types = retrieve_checkm_info(lpsn_types, api_key)
    for type_strain in lpsn_types:
        if type_strain.genome_acc is not None:
            type_strain.genome_acc['accession'] = type_strain.genome_acc['accession'].split(".")[0]

    if skip_download != "all" and skip_download != "genomes":
        # logger.info("Skipping sequence download steps as per user request.\n")
        retrieve_genome_sequences(lpsn_types, output_path, threads, api_key)
    if skip_download != "all" and skip_download != "16s":
        retrieve_16S_sequences(lpsn_types, output_path, email, input_taxon, api_key)

    logger.success("Metadata retrieval from GenBank complete.\n")

    return lpsn_types


def retrieve_sequences_workflow(lpsn_types, output_path, threads, dev_mode, input_taxon, skip_download):
    
    logger.info("Step 4 of 4: Generating ouputs")
    output_metadata(lpsn_types, output_path)
    if skip_download == "all":
        logger.info("Skipping sequence download steps as per user request.\n")
        logger.info("###################################")
        logger.info("###     Ratatoskr finished!     ###")
        logger.info("###################################\n")
        return
    rehydrate_genome_sequences(output_path, threads, input_taxon)
    tidy_genome_dir(output_path, input_taxon)

    logger.info("Sequence retrieval workflow complete.\n")

    logger.info("###################################")
    logger.info("###     Ratatoskr finished!     ###")
    logger.info("###################################\n")
