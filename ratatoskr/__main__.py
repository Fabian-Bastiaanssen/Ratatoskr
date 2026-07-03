#!/usr/bin/env python3
"""ratatoskr"""

######### Import libraries ###########

from pathlib import Path

import rich_click as click

from ratatoskr.utils import get_version
from ratatoskr.initialisation import set_up_logger, initialise_clients, fetch_SILVA_seqs
from ratatoskr.lpsn import retrieve_LPSN_type_info
from ratatoskr.genbank import retrieve_info_from_genbank, retrieve_sequences_workflow, get_genbank_api_info
from ratatoskr.bacdive import retrieve_extra_info_from_bacdive

######################################

###### define global variables #######


############ define main cli ############

@click.group()
@click.version_option(get_version(), "--version", "-v")
@click.help_option("--help", "-h")
def main_cli():
    """
    \b      
    8888888b.           888             888                      888              
    888    888          888             888                      888              
    888   d88P  8888b.  888888  8888b.  888888  .d88b.  .d8888b  888  888 888d888 
    8888888P"      "88b 888        "88b 888    d88""88b 88K      888 .88P 888P"   
    888 T88b   .d888888 888    .d888888 888    888  888 "Y8888b. 888888K  888     
    888  T88b  888  888 Y88b.  888  888 Y88b.  Y88..88P      X88 888 "88b 888     
    888   T88b "Y888888  "Y888 "Y888888  "Y888  "Y88P"   88888P' 888  888 888
    \b
    Search for and retrieve type strains info, 16S rRNA and genome sequences for 
    a given genus, family, or order. 
    """
    pass

#########################################

############ define options #############

@click.command()

@click.option(
        "-i",
        "--input",
        help="Name of the taxa of interest.",
        type=str,
        required=True,
)
@click.option(
        "-o",
        "--output_path",
        "--output",
        help="Specify the desired path for creation of the output folder",
        required=True,
        type=click.Path(exists=False, dir_okay=True, readable=True, path_type=Path)
    )
@click.option(
        "-l",
        "--level",
        help="Specify the taxanomic level to check for your taxon.",
        type=click.Choice(['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'auto'], case_sensitive=False),
        callback=lambda ctx, param, value: value.lower(),
        default = "auto",
        show_default = True,
        required=False,
    )
@click.option(
        "-t",
        "--threads",
        type=click.IntRange(min=1),
        help="Number of threads to be used.",
        required=False,
        default = 1,
        show_default = True
    )
@click.option(
        "-f",
        "--force",
        help="Force overwrite of output directory.",
        required=False,
        is_flag=True
    )
@click.option(
        "-d",
        "--dev_mode",
        help="Run in development mode.",
        required=False,
        is_flag=True
    )
@click.option(
        "-s",
        "--sequence_download",
        help="Specify the sequence type to download. All choices still return accessions, selection is to download associated FASTA files.",
        type=click.Choice(['all', '16S', 'genomes', 'none'], case_sensitive=False),
        callback=lambda ctx, param, value: value.lower(),
        show_default = True,
        default = 'all',
        required=False
    )
@click.option(
        "-m",
        "--mmseqs_16S_check",
        help="Specify whether to check 16S rRNA sequences with mmseqs against SILVA for either all sequences or only those from GenBank ('all' can be time-consuming step for higher taxonomic levels)",
        type=click.Choice(['all', 'genbank'], case_sensitive=False),
        callback=lambda ctx, param, value: value.lower(),
        show_default = True,
        default = 'all',
        required=False
    )
@click.option(
        "-n",
        "--no_cache",
        help="Don't use 16S rRNA cache information.",
        required=False,
        is_flag=True
    )
@click.option(
        "-c",
        "--cache",
        help="FOR DEV MODE ONLY: Path to create new cache file.",
        type=str,
        required=False,
)

#########################################

############ define commands ############

@click.version_option(get_version(), "--version", "-v")
@click.help_option("--help", "-h")
@click.pass_context

def run(ctx, input, output_path, threads, force, level, dev_mode, sequence_download, no_cache, cache, mmseqs_16s_check): 

    """
    Run the ratatoskr pipeline
    """

    set_up_logger(output_path, force, debug=dev_mode)
    if sequence_download == '16s' or sequence_download == 'all':
        fetch_SILVA_seqs(threads)
    email, api_key = get_genbank_api_info(dev_mode) 
    lpsn_client, bacdive_client = initialise_clients(dev_mode)
    lpsn_types = retrieve_LPSN_type_info(input, output_path, threads, level, lpsn_client, dev_mode, no_cache, cache, api_key=api_key, email=email)
    lpsn_types = retrieve_extra_info_from_bacdive(lpsn_types, bacdive_client)
    lpsn_types = retrieve_info_from_genbank(lpsn_types, output_path, threads, dev_mode, input, sequence_download, email=email, api_key=api_key, mmseqs_16s_check=mmseqs_16s_check)
    
####################################################################################################


#########################################
main_cli.add_command(run)

def main():
    main_cli()

if __name__ == '__main__':
    main()
