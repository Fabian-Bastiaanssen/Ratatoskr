from loguru import logger
import asyncio
import aiohttp

from ratatoskr.utils import make_dir, move_and_rename, delete_thing

def tidy_genome_dir(output_path, input_taxon):
    base = output_path / "sequences" / "genomes"
    fastas_dir = base / "fastas"
    genbanks_dir = base / "genbank"

    make_dir(fastas_dir)
    make_dir(genbanks_dir)

    data_root = base / "ncbi_dataset" / "data"

    move_and_rename(data_root / "*" / "*.fna", fastas_dir, "fasta")
    move_and_rename(data_root / "*" / "*.gbff", genbanks_dir, "genbank")
    

    for i in [base / "ncbi_dataset", base / "genomes.zip", base / "README.md", base / "md5sum.txt", output_path / "genome_accessions.txt"]:
        delete_thing(i)

async def fetch_data(url, session, headers, sem, mode, pbar=None, query_length=1):
    results = []
    next_url = url
    async with sem:
        while next_url:
            for attempt in range(1, 4):
                try:
                    async with session.get(next_url, headers=headers) as resp:
                        if resp.status in [429, 500, 502, 503, 504]:
                            await asyncio.sleep(2 ** attempt)
                            continue
                        data = await resp.json()
                        results.extend(data.get(mode, []))
                        next_token = data.get("next_page_token")
                        if next_token:
                            # logger.debug(f"Fetching next page with token {next_token}")
                            base_url = next_url.split("&page_token=")[0]
                            next_url = f"{base_url}&page_token={next_token}"
                        else:
                            next_url = None
                        break  # break retry loop
                except aiohttp.ClientError:
                    await asyncio.sleep(2 ** attempt)
            else:
                # all retry attempts failed
                if pbar:
                    pbar.update(query_length)
                return results  # or None
    if pbar:
        pbar.update(query_length)
    return results


def get_haves_and_have_nots(lpsn_types, attribute_name):
    has = []
    has_not = []
    for x in lpsn_types:
        if getattr(x, attribute_name) is not None:
            if getattr(x, attribute_name) != [] and getattr(x, attribute_name) != "":
                has.append(x)
            else:
                has_not.append(x)
        else:
            has_not.append(x)
    return has, has_not


def get_taxaomic_levels(input):

    if input.lower() in ["bacteria", "archaea", "miscellanea"]:
        taxonomic_level = "domain"
    elif input.lower().endswith(("ida", "ati")):
        taxonomic_level = "kingdom"
    elif input.lower().endswith("ota"):
        taxonomic_level = "phylum"
    elif input.lower().endswith("aceae"):
        taxonomic_level = "family"
    elif input.lower().endswith(("li", "iia", "lia", "cutes", 
                                 "eae", "bacteria", "cia", "lobi", 
                                 "enia", "adia", "ditrichia", "phagia",
                                 "flexia", "mycetes")):
        taxonomic_level = "class"
    elif input.lower().endswith("ales"):
        taxonomic_level = "order"
    else:
        taxonomic_level = "genus"

    return taxonomic_level
