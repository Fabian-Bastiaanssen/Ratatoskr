from importlib import resources
import sys

from loguru import logger
from async_dsmz import lpsn_async
import pickle as pkl
import tqdm
import os
import asyncio

os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl

from ratatoskr.utils import suppress_stdout
from ratatoskr.misc import get_taxaomic_levels
from ratatoskr.type_strain import TypeStrain
from ratatoskr.genbank import get_acc_seq_lengths
accepted_ranks = [
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies"
]

bad_ids = [37563, 37583, 65201, 37506, 6001, 
           29696, 6002, 6000, 29695, 6009, 6025, 
           6024, 6023, 62824, 6007, 6031, 6030, 
           6013, 29698, 29697, 29699, 6019, 6028, 
           6029, 6021, 6006, 29701, 6016, 6010, 
           46673, 29700, 6008]

def set_lpsn_client(email, password):
    """
    Set up the LPSN client.
    """
    
    logger.debug(f"Setting up LPSN client with email {email}.")
    
    try:
        with suppress_stdout():
            lpsn_client = lpsn_async(email, password)
        logger.success("LPSN client set up.")
    except Exception as e:
        logger.error(f"LPSN client could not be set up: {e}")
        sys.exit(1)
    
    return lpsn_client


def match_parent(main_df, child_df, level, previous_level):
    """
    Match parent taxon IDs from main_df to child_df and aggregate relevant info.
    """
    if "species" in level:
            child_df = (child_df.with_columns(pl.col('lpsn_taxonomic_status').str.starts_with('correct name').alias('correct'),
                                            pl.coalesce(pl.col('lpsn_correct_name_id'), pl.col('id')).alias(f'parent_{level}_id')).sort(by='correct', descending=True))
            child_df = (child_df.group_by(f'parent_{level}_id')
                    .agg(pl.col('lpsn_parent_id').first(),
                        pl.col('full_name').first().alias(f"parent_{level}"),
                        pl.col("full_name").alias(f"binomial_synonyms"),
                        pl.col('lpsn_taxonomic_status'),
                        pl.col('type_strain_names').first().alias(f"type_names"),
                        pl.col('correct'),
                        pl.col('molecules').list.eval(pl.element().struct.field('identifier')).list.first(),
                        pl.col('ijsem_list_doi').first().alias('list_ref'),
                        pl.col('publication_doi').first().alias('pub_ref'),
                        pl.col('authority').first().alias('authority')))
            
            child_df = (child_df.filter(pl.col('correct').list.any()).with_columns(pl.col('molecules').list.filter(pl.element().str.len_chars() < 12).list.first().alias('rRNA_acc'))
                    .drop(['correct', 'lpsn_taxonomic_status']))
    else:
            child_df = (child_df.with_columns(pl.col('lpsn_taxonomic_status').str.starts_with('correct name').alias('correct'),
                                            pl.coalesce(pl.col('lpsn_correct_name_id'), pl.col('id')).alias(f'parent_{level}_id')))

            child_df = (child_df.sort(by='correct', descending=True).group_by(f'parent_{level}_id')
                    .agg(pl.col('lpsn_parent_id').first(),
                        pl.col('full_name').first().alias(f"parent_{level}"),
                        pl.col('lpsn_taxonomic_status'),
                        pl.col('correct'))) 
            child_df = (child_df
                    .drop(['correct', 'lpsn_taxonomic_status'])
                )
    main_df = main_df.join(child_df, left_on=f'parent_{previous_level}_id', right_on='lpsn_parent_id', how='left')
    return main_df


def coalesce_all_lpsn_taxon_df(df):
    return (
        df.with_columns(pl.coalesce(pl.col('type_names_right'), pl.col('type_names')).alias('type_names'),
                        pl.coalesce(pl.col('rRNA_acc_right'), pl.col('rRNA_acc')).alias('rRNA_acc'),
                        pl.coalesce(pl.col('list_ref_right'), pl.col('list_ref')).alias('list_ref'),
                        pl.coalesce(pl.col('pub_ref_right'), pl.col('pub_ref')).alias('pub_ref'),
                        pl.coalesce(pl.col('binomial_synonyms_right'), pl.col('binomial_synonyms')).alias('binomial_synonyms'),
                        pl.coalesce(pl.col('authority_right'), pl.col('authority')).alias('authority'))
            .drop(['type_names_right', 'rRNA_acc_right', 'list_ref_right', 'molecules', 
                   'molecules_right', 'pub_ref_right', 'binomial_synonyms_right', 'authority_right'])
        )


def search_all_lpsn(lpsn_client):
    """
    Search LPSN for all taxonomic levels and compile a full DataFrame.
    """
    logger.info("Retrieving full LPSN taxonomy data.")
    
    with tqdm.tqdm(accepted_ranks, desc="LPSN taxonomic levels", unit="level", ncols=100, colour="magenta") as pbar:
        for taxon in pbar:
            pbar.set_postfix(level = taxon.capitalize())
            if lpsn_client.search(category=taxon) > 0:
                hits = lpsn_client.retrieve()
                if taxon == 'domain':
                    df = (
                        pl.LazyFrame(hits).rename({'id': 'parent_domain_id'})
                            .select([pl.col('parent_domain_id'), pl.col('full_name').alias('parent_domain')])
                    )
                else:
                    df = match_parent(df, pl.LazyFrame(hits), taxon, previous_level)
                previous_level = taxon
            else:
                logger.error(f"No hits found for taxonomic level {taxon} in LPSN. This sould not happen, please report this issue.")
                sys.exit(1)
    
    df = coalesce_all_lpsn_taxon_df(df)
    return df


def filter_dataframe(df, input, taxonomic_level):
    """
    Check if the input taxon is a synonym in LPSN DataFrame.
    """
    logger.info(f"Filtering LPSN data for {taxonomic_level} {input}.")
    lpsn_hits = (
        df.filter(pl.col(f"parent_{taxonomic_level}").str.to_lowercase() == input.lower())
          .filter(~((pl.col(f"parent_species").is_null()) & (pl.col(f"parent_subspecies").is_null())))
    ).collect()

    
    if len(lpsn_hits) == 0:
        
        logger.error(f"Could not uniquely identify {input} as {taxonomic_level} on LPSN.")
        # sys.exit(1)
        level_df = (df.lazy()
        # 1) build a per-column equality predicate
        .with_columns(
            [
                (pl.col(f"parent_{lvl}").str.to_lowercase() == input.lower())
                .alias(f"match_{lvl}")
                for lvl in accepted_ranks
            ]
        )
        # 2) filter rows where any of the taxonomic levels match
        .filter(pl.any_horizontal([pl.col(f"match_{lvl}") for lvl in accepted_ranks]))
        # 3) keep only rows where at least one of species/subspecies is non-null
        .filter(
            ~(
                (pl.col("parent_species").is_null())
                & (pl.col("parent_subspecies").is_null())
            )
        ).collect())
        # 4) add a column indicating which taxonomic level matched
        if len(level_df) == 0:
            logger.error(f"No matches found for {input} at any taxonomic level on LPSN.")
            print(df.filter(pl.col(f"parent_class").str.to_lowercase() == input.lower()).collect())
            sys.exit(1)
        else:
            taxonomic_level =level_df.with_columns(
                pl.when(pl.col("match_species")).then(pl.lit("species"))
                .when(pl.col("match_subspecies")).then(pl.lit("subspecies"))
                .when(pl.col("match_genus")).then(pl.lit("genus"))
                .when(pl.col("match_family")).then(pl.lit("family"))
                .when(pl.col("match_order")).then(pl.lit("order"))
                .when(pl.col("match_class")).then(pl.lit("class"))
                .when(pl.col("match_phylum")).then(pl.lit("phylum"))
                .when(pl.col("match_kingdom")).then(pl.lit("kingdom"))
                .when(pl.col("match_domain")).then(pl.lit("domain"))
                .otherwise(pl.lit(None))
                .alias("matched_level")
            )['matched_level'].first()
        if taxonomic_level:
            logger.info(f"Input taxon {input} matched at level {taxonomic_level}. Retrying filter.")
            lpsn_hits = (
                df.filter(pl.col(f"parent_{taxonomic_level}").str.to_lowercase() == input.lower())
                  .filter(~((pl.col(f"parent_species").is_null()) & (pl.col(f"parent_subspecies").is_null())))
            ).collect()
            if len(lpsn_hits) == 0:
                logger.error(f"Could not uniquely identify {taxonomic_level} as {input} on LPSN after second attempt.")
                sys.exit(1)
        else:
            logger.error(f"Could not uniquely identify {taxonomic_level} as {input} on LPSN after second attempt.")
            sys.exit(1)
   
    return lpsn_hits


def polars_to_type_strain_list(df):
    logger.debug("Converting LPSN DataFrame to list of TypeStrain objects.")
    return [TypeStrain(**row) for row in df.rows(named=True)]


def load_rRNA_cache(cache = None, no_cache=False):
    logger.info(f"Loading rRNA accession cache from {cache}.")
    if no_cache:
        logger.info("Skipping cache loading due to --no_cache flag.")
        return {"good": set(), "bad": set()}
    
    if cache is None:
        rRNA_cache_path = resources.files("ratatoskr.data").joinpath("rRNA_cache.pkl")
        if rRNA_cache_path.exists():
            logger.info("No cache path provided, loading default cache from package resources.")
            with rRNA_cache_path.open("rb") as f:
                rRNA_cache = pkl.load(f)
        else:
            logger.warning("No rRNA cache path provided and default rRNA cache file does not exist. Starting with empty cache.")
            rRNA_cache = {"good": set(), "bad": set()}
    else:
        if os.path.exists(cache):
            logger.debug(f"Loading rRNA accession cache from {cache}.")
            with open(cache, "rb") as f:
                rRNA_cache = pkl.load(f)
        else:
            logger.warning(f"Cache file {cache} does not exist. Starting with empty cache.")
            rRNA_cache = {"good": set(), "bad": set()}

    return rRNA_cache


def check_lpsn_rRNA_accs(lpsn_hits, dev_mode, no_cache, cache=None, api_key=None, email=None):
    
    logger.info("Checking LPSN rRNA accessions for validity.")

    rRNA_cache = load_rRNA_cache(cache, no_cache)

    all_accs = set([hit.rRNA_acc for hit in lpsn_hits if hit.rRNA_acc is not None])
    accs_needing_checked = all_accs - rRNA_cache["good"] - rRNA_cache["bad"]
   
    good_lengths = set()
    if len(accs_needing_checked) > 0:
        good_lengths = get_acc_seq_lengths(accs_needing_checked, api_key=api_key, email=email)
    good_lengths = set(good_lengths) | rRNA_cache["good"]
   
    updated_hits = []
    for hit in lpsn_hits:
        if hit.rRNA_acc is not None:
            if hit.rRNA_acc.split(".")[0] not in good_lengths:
                logger.debug(f"rRNA accession {hit.rRNA_acc} for {hit.parent_species} type strain {hit.type_names[0]} is shorter than 1000 bp or longer than 2000 bp, which may indicate an issue with the sequence. Removing.")
                hit.rRNA_acc = None
        updated_hits.append(hit)

    if dev_mode and cache is not None:
        logger.info("Updating rRNA accession cache with new information.")
        rRNA_cache["good"] = good_lengths
        rRNA_cache["bad"] = all_accs - good_lengths
        logger.debug(f"Saving rRNA accession cache to {cache}.")
        with open(cache, "wb") as f:
            pkl.dump(rRNA_cache, f)
    
    return updated_hits

   
def retrieve_LPSN_type_info(input, output_path, threads, level, lpsn_client, dev_mode, no_cache, cache, api_key=None, email=None):
 
    logger.info("Step 1 of 4: Retrieving taxonomical data from LPSN.")
    logger.info("No cached LPSN hits found, retrieving from LPSN API.")
    taxonomic_level = get_taxaomic_levels(input) if level == "auto" else level
    full_lpsn_df = search_all_lpsn(lpsn_client)
    lpsn_hits = filter_dataframe(full_lpsn_df, input, taxonomic_level)
    lpsn_hits = polars_to_type_strain_list(lpsn_hits)
    lpsn_hits = check_lpsn_rRNA_accs(lpsn_hits, dev_mode, no_cache, cache, api_key, email)
    logger.success(f"Retrieved metadata for {len(lpsn_hits)} type strains from LPSN.\n")
    return lpsn_hits
