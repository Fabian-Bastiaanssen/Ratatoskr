import re
import sys
from typing_extensions import final
import os
import pickle as pkl
from async_dsmz import bacdive_async
from loguru import logger
from tqdm import tqdm
import asyncio

import traceback

from ratatoskr.utils import suppress_stdout

assembly_levels = ("complete", "chromosome", "scaffold", "contig")

type_strain_attributes = ["genome_acc", "isolation_sample_type", "culture_pH", "culture_temp", 
 "culture_medium", "morphology_pigmentation", "morphology_cell", 
 "morphology_colony", "API_results", "fatty_acid_profile", 
 "oxygen_tolerance", "spore_formation", "compound_production", 
 "metabolite_production", "metabolite_utilization"]

def set_bacdive_client(email, password):
    """
    Set up the BacDive client.
    """
    logger.debug(f"Setting up BacDive client with email {email}.")
    try:
        with suppress_stdout():
            bacdive_client = bacdive_async(email, password)
        logger.success("Bacdive client set up.\n")
    except Exception as e:
        logger.error(f"Bacdive client could not be set up: {e}")
        sys.exit(1)
    return bacdive_client


def split_lpsn_types_list_to_genus_dict(lpsn_types):
    genus_dict = {}
    for type_strain in lpsn_types:
        genus = type_strain.parent_genus
        if genus == None:
            logger.warning(f"Type strain {type_strain} has no assigned genus. Skipping BacDive retrieval for this strain.")
            continue
        if genus not in genus_dict:
            genus_dict[genus] = []
        genus_dict[genus].append(type_strain)
    return genus_dict

def get_sample_isolation_info(record):
    
    isolation_sample = record.get('Isolation, sampling and environmental information', {}).get('isolation', {})

    isolation_sample_type = set()
    
    if isinstance(isolation_sample, dict):
        if 'sample type' in isolation_sample:
            isolation_sample_type.add(isolation_sample.get('sample type'))
    elif isinstance(isolation_sample, list):
        isolation_sample_type.update([x.get('sample type') for x in isolation_sample if 'sample type' in x])
    else:
        isolation_sample_type = None

    if isolation_sample_type:
        isolation_sample_type = ';'.join(isolation_sample_type)
    
    return isolation_sample_type


def get_strain_designations(record):
    
    culture_ids = [x.strip() for x in record.get("External links", {}).get("culture collection no.", "").split(',')]
    # the re below is because some strain designations have a trailing ' T' to indicate type strain for some reason
    strain_designations = [re.sub(r" T$", "", x.strip()) for x in record.get("Name and taxonomic classification", {}).get("strain designation", "").split(',') if x.strip() != '']
    final_names = set(culture_ids + strain_designations)
    
    return final_names


def get_ncbi_tax_id(record):

    taxid_info = record.get("General", {}).get("NCBI tax id")
    
    if taxid_info is not None:
        if type(taxid_info) is list:
            return [(x.get("NCBI tax id"), x.get("Matching level")) for x in taxid_info]
        else:
            return [(taxid_info.get("NCBI tax id"), taxid_info.get("Matching level"))]
    else:
        return None

def get_16S_sequence_info(record):

    rRNA_info = record.get("Sequence information", {}).get("16S sequences")
    
    if rRNA_info is not None:
        if type(rRNA_info) is list:
            rRNA_list = [{"accession": x.get("accession").split(".")[0], "sequence": x.get("sequence"), "length": x.get("length")} for x in rRNA_info if x.get("database") == "nuccore" and x.get("length") is not None and x.get("length") <= 2000 and len(x.get("accession")) < 12]
            if len(rRNA_list) > 0:
                rRNA_list = sorted(rRNA_list, key=lambda x: x["length"], reverse=True)
                rRNA_acc = rRNA_list[0]
            else:
                rRNA_acc = None
        else:
            if rRNA_info.get("database") == "nuccore" and rRNA_info.get("length") is not None and rRNA_info.get("length") <= 2000 and len(rRNA_info.get("accession")) < 12:
                rRNA_acc = {"accession": rRNA_info.get("accession").split(".")[0], "sequence": rRNA_info.get("sequence"), "length": rRNA_info.get("length")}
            else:
                rRNA_acc = None
    else:
        rRNA_acc = None

    return rRNA_acc

def get_genome_sequence_info(record):

    genome_info = record.get("Sequence information", {}).get("Genome sequences")
    
    if genome_info is not None:
        if type(genome_info) is list:
            genome_list = [{"accession": x.get("accession"), "assembly level": x.get("assembly level")} for x in genome_info if x.get("database") == "ncbi" and x.get("assembly level") in assembly_levels]
            if len(genome_list) > 0:
                genome_list = sorted(genome_list, key=lambda x: assembly_levels.index(x.get("assembly level")))
                genome_acc = genome_list[0]
            else:
                genome_acc = None
        else:
            if genome_info.get("database") == "ncbi" and genome_info.get("assembly level") in assembly_levels:
                genome_acc = {"accession": genome_info.get("accession"), "assembly level": genome_info.get("assembly level")}
            else:
                genome_acc = None
    else:
        genome_acc = None

    return genome_acc

def get_culture_pH(record):

    culture_pH = record.get("Culture and growth conditions", {}).get("culture pH")

    if isinstance(culture_pH, dict) and culture_pH.get("ability") in ["positive"]:
        try:
            culture_pH = str(culture_pH.get("pH").rstrip("."))
            if '<' in culture_pH:
                culture_pH = f"<{float(culture_pH.replace('<',''))}"
            if '>' in culture_pH:
                culture_pH = f">{float(culture_pH.replace('>',''))}"
            
            if '-' in culture_pH:
                culture_pH = re.split(r'(?<=[0-9.])-(?=-?\d)', culture_pH)
                culture_pH = sorted([float(x.rstrip(".")) for x in culture_pH])
                culture_pH = "-".join([str(x) for x in culture_pH])
            elif '<' not in culture_pH and '>' not in culture_pH:
                culture_pH = float(culture_pH)
        except Exception as e:
            logger.warning(f"Could not parse pH value {culture_pH} in record {record} as entry appears malformed. Setting pH to None.")
            logger.debug(f"Error details: {e}")
            culture_pH = None 
    elif isinstance(culture_pH, list):
        min = None
        max = None
        optimum = []
        for i in culture_pH:
            if isinstance(i, dict) and i.get("type") in ["growth", "optimum", "minimum", "maximum"] and i.get("pH") is not None and i.get("ability") == "positive":
                try:
                    pH_value = i.get("pH").rstrip(".")
                    if '<' in str(pH_value):
                        pH_value = float(str(pH_value).replace('<',''))
                    if '>' in str(pH_value):
                        pH_value = float(str(pH_value).replace('>',''))
                    try:
                        pH_value = [float(x.rstrip(".")) for x in re.split(r'(?<=[0-9.])-(?=-?\d)', str(pH_value), maxsplit=1)]
                    except:
                        logger.warning(f"Could not parse pH value '{pH_value}' in {record.get("Name and taxonomic classification", {}).get("LPSN", {}).get("species")} record as entry appears malformed. Skipping this pH value.")
                        continue
                    for pH in pH_value:
                        if min is None or pH < min:
                            min = pH
                        if max is None or pH > max:
                            max = pH
                        if i.get("type") == "optimum":
                            optimum.append(pH)
                except Exception as e:
                    logger.warning(f"Could not parse pH value '{pH_value}' in {record.get("Name and taxonomic classification", {}).get("LPSN", {}).get("species")} record as entry appears malformed. Skipping this pH value.")
                    logger.debug(f"Error details: {e}")
        optimum = list(set(optimum) - set([min, max]))
        if len(optimum) > 0:
            culture_pH = f'{min}-({';'.join([str(x) for x in sorted(optimum)])})*-{max}'
        elif min is not None and max is not None and min != max:
            culture_pH = f"{min}-{max}"
        else:
            culture_pH = min       
    else:
        culture_pH = None

    return culture_pH


def get_culture_temp(record):

    culture_temp = record.get("Culture and growth conditions", {}).get("culture temp", {})
    
    if isinstance(culture_temp, dict) and culture_temp != {} and culture_temp.get("growth") in ["positive"]:
        if '-' in str(culture_temp.get("temperature")):
            culture_temp = re.split(r'(?<=[0-9.])-(?=-?\d)', culture_temp.get("temperature"), maxsplit=1)
            culture_temp = sorted([float(x.rstrip(".")) for x in culture_temp])
            culture_temp = "-".join([str(x) for x in culture_temp])
        else:
            culture_temp = float(culture_temp.get("temperature").rstrip(".")) 
    elif isinstance(culture_temp, list):
        min = None
        max = None
        optimum = []
        for i in culture_temp:
            if isinstance(i, dict) and i.get("type") in ["growth", "optimum"] and i.get("temperature") is not None and i.get("growth") == "positive":
                temperature_value = i.get("temperature")
                temperature_value = [float(x.rstrip(".")) for x in re.split(r'(?<=[0-9.])-(?=-?\d)', str(temperature_value), maxsplit=1)]
                for temp in temperature_value:
                    if min is None or temp < min:
                        min = temp
                    if max is None or temp > max:
                        max = temp
                    if i.get("type") == "optimum":
                        optimum.append(temp)
        optimum = list(set(optimum) - set([min, max]))
        if len(optimum) > 0:
            culture_temp = f'{min}-({';'.join([str(x) for x in sorted(optimum)])})*-{max}'
        elif min is not None and max is not None and min != max:
            culture_temp = f"{min}-{max}"
        else:
            culture_temp = min       
    else:
        culture_temp = None

    return culture_temp


def get_culture_conditions(record):

    culture_medium = record.get("Culture and growth conditions", {}).get("culture medium")

    culture_medium = set()

    if isinstance(culture_medium, dict) and culture_medium.get("growth") == "yes":
        culture_medium.add(culture_medium.get("name"))
    elif isinstance(culture_medium, list):
        culture_medium.update([x.get("name") for x in culture_medium if isinstance(x, dict) and x.get("growth") == "yes" and x.get("name") is not None])
    else:
        culture_medium = None

    if culture_medium:
        culture_medium = ';'.join(culture_medium)

    return culture_medium

def get_morphology_pigmentation(record):

    pigmentation = record.get('Morphology', {}).get('pigmentation')
    pigmentation_val = []
    if pigmentation:
        if isinstance(pigmentation, dict) and pigmentation.get('production') == 'yes' and (pigmentation.get('name') is not None or pigmentation.get('colour') is not None):
            name = pigmentation.get('name', 'ND')
            colour = pigmentation.get('colour', 'ND')
            pigmentation_val.append(f"name:{name}~colour:{colour}")
        elif isinstance(pigmentation, list):
            for i in pigmentation:
                if isinstance(i, dict):
                    if i.get('production') == 'yes' and (i.get('name') is not None or i.get('colour') is not None):
                        name = i.get('name', 'ND')
                        colour = i.get('colour', 'ND')
                        pigmentation_val.append(f"name:{name}~colour:{colour}")
    
    if pigmentation_val == []:
        pigmentation_val = None
    else:
        pigmentation_val = ';'.join(pigmentation_val)

    return pigmentation_val

def get_morphology_cell(record):

    cell_morphology = record.get('Morphology', {}).get('cell morphology')
    
    cell_width = []
    cell_length = []
    gram_stain = []
    cell_shape = []
    motility = []
    flagellum_arrangement = []

    if cell_morphology:
        if isinstance(cell_morphology, dict):
            cell_width.append(cell_morphology.get('cell width'))
            cell_length.append(cell_morphology.get('cell length'))
            gram_stain.append(cell_morphology.get('gram stain'))
            cell_shape.append(cell_morphology.get('cell shape'))
            motility.append(cell_morphology.get('motility'))
            flagellum_arrangement.append(cell_morphology.get('flagellum arrangement'))
        elif isinstance(cell_morphology, list):
            for i in cell_morphology:            
                if isinstance(i, dict):
                    cell_length.append(i.get('cell length'))
                    cell_width.append(i.get('cell width'))
                    gram_stain.append(i.get('gram stain'))
                    cell_shape.append(i.get('cell shape'))
                    motility.append(i.get('motility'))
                    flagellum_arrangement.append(i.get('flagellum arrangement'))

    for var in [cell_width, cell_length, gram_stain, cell_shape, motility, flagellum_arrangement]:
        var[:] = list(set([x for x in var if x is not None]))

    cell_morphology_val = {
        "cell_width": ';'.join([str(x) for x in cell_width]) if len(cell_width) > 0 else None,
        "cell_length": ';'.join([str(x) for x in cell_length]) if len(cell_length) > 0 else None,
        "gram_stain": ';'.join([str(x) for x in gram_stain]) if len(gram_stain) > 0 else None,
        "cell_shape": ';'.join([str(x) for x in cell_shape]) if len(cell_shape) > 0 else None,
        "motility": ';'.join([str(x) for x in motility]) if len(motility) > 0 else None,
        "flagellum_arrangement": ';'.join([str(x) for x in flagellum_arrangement]) if len(flagellum_arrangement) > 0 else None,
    }    

    return cell_morphology_val

def get_morphology_colony(record):

    colony_morphology = record.get('Morphology', {}).get('colony morphology', {})

    possible_morphs = {}
    description = None

    if colony_morphology:
        if isinstance(colony_morphology, dict):
                description = {'colony_shape':[], 'type_of_hemolysis':[], 'colony_size':[], 'colony_color':[], 'incubation_period':[]}
                description['colony_shape'].append(colony_morphology.get('colony shape'))
                description['type_of_hemolysis'].append(colony_morphology.get('type of hemolysis')) 
                description['colony_size'].append(colony_morphology.get('colony size'))
                description['colony_color'].append(colony_morphology.get('colony color'))
                description['incubation_period'].append(colony_morphology.get('incubation period'))
                if colony_morphology.get('medium used') is not None:
                    possible_morphs[colony_morphology.get('medium used')] = description
                else:                   
                    possible_morphs['ND'] = description
        
        elif isinstance(colony_morphology, list):
            for i in colony_morphology:            
                if isinstance(i, dict):
                    description = {'colony_shape':[], 'type_of_hemolysis':[], 'colony_size':[], 'colony_color':[], 'incubation_period':[]}
                    description['colony_shape'].append(i.get('colony shape'))
                    description['type_of_hemolysis'].append(i.get('type of hemolysis')) 
                    description['colony_size'].append(i.get('colony size'))
                    description['colony_color'].append(i.get('colony color'))
                    description['incubation_period'].append(i.get('incubation period'))
                if description:
                    medium = i.get('medium used', 'ND')
                    if medium is not None:
                        if medium not in possible_morphs:
                            possible_morphs[medium] = description
                        elif description and medium in possible_morphs:
                            for key in description:
                                possible_morphs[medium][key].extend(description[key])
    
    if possible_morphs == {}:
        possible_morphs = None

    return possible_morphs

def get_API_results(record): 

    phenotype = record.get('Physiology and metabolism', {})
    API_keys = [k for k in record.get('Physiology and metabolism', {}) if isinstance(k, str) and k.lower().startswith('api')]
    
    API_results = {}

    for key in API_keys:
        if isinstance(phenotype.get(key), dict):
            API_results[key] = {k: [v] for k, v in phenotype.get(key, {}).items() if not str(k).startswith('@')}
        elif isinstance(phenotype.get(key), list):
            API_results[key] = {}
            for i in phenotype.get(key, []):
                if isinstance(i, dict):
                    for k, v in i.items():
                        if not str(k).startswith('@'):
                            if k not in API_results[key]:
                                API_results[key][k] = set(v)
                            else:
                                API_results[key][k].add(v)
            for k in API_results[key]:
                API_results[key][k] = list(API_results[key][k])
    
    if API_results == {}:
        return None

    return API_results

def get_fatty_acid_profile(record):
    
    profiles = record.get('Physiology and metabolism', {}).get("fatty acid profile", {})

    if profiles is None:
        return None

    fatty_acid_profile_results = {}

    if isinstance(profiles, dict) and profiles != {} and 'fatty acids' in profiles:
        if isinstance(profiles.get('fatty acids'), dict):
            name = profiles.get('fatty acids').get('fatty acid')
            percentage = profiles.get('fatty acids').get('percentage')
            fatty_acid_profile_results[name] = str(percentage)
        elif isinstance(profiles.get('fatty acids'), list):
            for profile in profiles.get('fatty acids'):
                if isinstance(profile, dict):
                    name = profile.get('fatty acid')
                    percentage = profile.get('percentage')
                    if name and percentage is not None:
                        if name not in fatty_acid_profile_results:
                            fatty_acid_profile_results[name] = str(percentage)
                        else:
                            fatty_acid_profile_results[name] += f";{str(percentage)}"
                else:
                    logger.warning(f"Unexpected fatty acid profile format: {profiles.get('fatty acids')}. Skipping this entry.")
        else:
            fatty_acid_profile_results = None

    if fatty_acid_profile_results == {}:
        fatty_acid_profile_results = None
    
    return fatty_acid_profile_results

def get_oxygen_tolerance(record):
    oxygen_tolerance = record.get('Physiology and metabolism', {}).get('oxygen tolerance')

    oxygen_tolerance_result = set()
    
    if oxygen_tolerance is None:
        return None
    
    if isinstance(oxygen_tolerance, dict) and oxygen_tolerance.get('oxygen tolerance') is not None:
        oxygen_tolerance_result.add(oxygen_tolerance.get('oxygen tolerance'))
    elif isinstance(oxygen_tolerance, list):
        for i in oxygen_tolerance:
            if isinstance(i, dict) and i.get('oxygen tolerance') is not None:
                oxygen_tolerance_result.add(i.get('oxygen tolerance'))
            else:
                logger.warning(f"Unexpected oxygen tolerance format: {i}. Skipping this entry.")
    else:
        oxygen_tolerance_result = None

    if oxygen_tolerance_result == set() or oxygen_tolerance_result is None:
        oxygen_tolerance_result = None
    else:
        oxygen_tolerance_result = ';'.join(oxygen_tolerance_result)

    return oxygen_tolerance_result

def get_spore_formation(record):

    spore_formation = record.get('Physiology and metabolism', {}).get('spore formation')

    spore_formation_result = set()

    if spore_formation is None:
        return None

    if isinstance(spore_formation, dict):
        spore_formation_result.add(spore_formation.get('spore formation'))
    elif isinstance(spore_formation, list):
        for i in spore_formation:
            if isinstance(i, dict):
                spore_formation_result.add(i.get('spore formation'))
            else:
                logger.warning(f"Unexpected spore formation format: {i}. Skipping this entry.")
    else:
        spore_formation_result = None
    
    if spore_formation_result == set():
        spore_formation_result = None
    else:
        spore_formation_result = ';'.join(spore_formation_result)
        
    return spore_formation_result

def get_compound_production(record):

    compound_production = record.get('Physiology and metabolism', {}).get('compound production', {})

    compound_production_result = set()

    if compound_production is None or compound_production == {}:
        return None

    if isinstance(compound_production, dict):
        compound_production_result.add(compound_production.get('compound'))
    elif isinstance(compound_production, list):
        for i in compound_production:
            if isinstance(i, dict):
                compound_production_result.add(i.get('compound'))
            else:
                logger.warning(f"Unexpected compound production format: {i}. Skipping this entry.")
    else:
        compound_production_result = None
    
    if compound_production_result == set() or compound_production_result is None or compound_production_result == {None}:
        compound_production_result = None
    else:
        compound_production_result = ';'.join(compound_production_result)
    return compound_production_result

def get_metabolite_production(record):

    metabolite_production = record.get('Physiology and metabolism', {}).get('metabolite production', {})

    if metabolite_production is None or metabolite_production == {}:
        return None

    metabolite_production_result = {}

    if isinstance(metabolite_production, dict):
        metabolite_production_result[metabolite_production.get('metabolite')] = metabolite_production.get('production')
    elif isinstance(metabolite_production, list):
        for i in metabolite_production:
            if isinstance(i, dict):
                if i.get('metabolite') is not None and i.get('production') is not None:
                    if i.get('metabolite') not in metabolite_production_result:
                        metabolite_production_result[i.get('metabolite')] = set([i.get('production')])
                    else:
                        metabolite_production_result[i.get('metabolite')].add(i.get('production'))
            else:
                logger.warning(f"Unexpected metabolite production format: {i}. Skipping this entry.")
    else:
        metabolite_production_result = None

    if metabolite_production_result == {} or metabolite_production_result is None or metabolite_production_result == {None}:
        return None

    for k in metabolite_production_result:
        if isinstance(metabolite_production_result[k], set):
            metabolite_production_result[k] = ';'.join(metabolite_production_result[k])
    
    return metabolite_production_result

def get_metabolite_utilization(record):

    metabolite_utilization = record.get('Physiology and metabolism', {}).get('metabolite utilization', {})

    if metabolite_utilization is None or metabolite_utilization == {}:
        return None

    metabolite_utilization_result = {}

    if isinstance(metabolite_utilization, dict):
        metabolite = metabolite_utilization.get('metabolite')
        activity = metabolite_utilization.get('utilization activity')
        utilisation_type = metabolite_utilization.get('kind of utilization tested')
        if metabolite is not None and activity is not None:
            if utilisation_type is not None:
                metabolite_utilization_result[metabolite] = {"activity": [activity], "utilisation_type": [utilisation_type]}
            else:
                metabolite_utilization_result[metabolite] = {"activity": [activity], "utilisation_type": ["ND"]}
    elif isinstance(metabolite_utilization, list):
        for i in metabolite_utilization:
            if isinstance(i, dict):
                metabolite = i.get('metabolite')
                activity = i.get('utilization activity')
                utilisation_type = i.get('kind of utilization tested')
                if metabolite is not None and activity is not None:
                    if utilisation_type is not None:
                        if metabolite not in metabolite_utilization_result:
                            metabolite_utilization_result[metabolite] = {"activity": [activity], "utilisation_type": [utilisation_type]}
                        else:
                            metabolite_utilization_result[metabolite]["activity"].append(activity)
                            metabolite_utilization_result[metabolite]["utilisation_type"].append(utilisation_type)
                    else:
                        if metabolite not in metabolite_utilization_result:
                            metabolite_utilization_result[metabolite] = {"activity": [activity], "utilisation_type": ["ND"]}
                        else:
                            metabolite_utilization_result[metabolite]["activity"].append(activity)
                            metabolite_utilization_result[metabolite]["utilisation_type"].append("ND")
            else:
                logger.warning(f"Unexpected metabolite utilization format: {i}. Skipping this entry.")
    else:
        metabolite_utilization_result = None

    if metabolite_utilization_result == {} or metabolite_utilization_result is None or metabolite_utilization_result == {None}:
        return None
    
    return metabolite_utilization_result


def parse_bacdive_records(records):

    parsed_records = []
    
    for record in records:

        parsed_record = {}
       
        parsed_record["type_names"] = get_strain_designations(record)
        parsed_record['ncbi_tax_id'] = get_ncbi_tax_id(record)
        parsed_record["16S_rRNA"] = get_16S_sequence_info(record)
        parsed_record['isolation_sample_type'] = get_sample_isolation_info(record)
        parsed_record["genome_acc"] = get_genome_sequence_info(record)
        parsed_record["culture_pH"] = get_culture_pH(record)
        parsed_record["culture_temp"] = get_culture_temp(record)
        parsed_record["culture_medium"] = get_culture_conditions(record)
        parsed_record["morphology_pigmentation"] = get_morphology_pigmentation(record)
        parsed_record["morphology_cell"] = get_morphology_cell(record)
        parsed_record["morphology_colony"] = get_morphology_colony(record)
        parsed_record["API_results"] = get_API_results(record)
        parsed_record["fatty_acid_profile"] = get_fatty_acid_profile(record)
        parsed_record["oxygen_tolerance"] = get_oxygen_tolerance(record)
        parsed_record["spore_formation"] = get_spore_formation(record)
        parsed_record["compound_production"] = get_compound_production(record)
        parsed_record["metabolite_production"] = get_metabolite_production(record)
        parsed_record["metabolite_utilization"] = get_metabolite_utilization(record)
                                                                      
        # 'Physiology and metabolism' values that could be added:
        ### ['antibiogram', 'antibiotic resistance', 'observation', 'tolerance', 'nutrition type', 'halophily', 'metabolite tests', 'enzymes', 'murein']
        
        parsed_records.append(parsed_record)
    
    return parsed_records


async def fetch_genus_data(genus, type_strains, bacdive_client, pbar):
    
    lpsn_types = []
    bacdive_client.setSearchType("exact")
    
    records = [x for x in await bacdive_client.retrieve_async(f"taxon/{genus}")
                if x.get("Name and taxonomic classification", {}).get("type strain") == "yes"]
    
    if len(records) == 0:
        logger.info(f"No BacDive type strain data found for genus {genus}. This may because it is only recently described.")
        pbar.update(1)
        return type_strains
    
    try:
        thing = set()
        parsed_records = parse_bacdive_records(records)

    except Exception as e:
        logger.error(f"Error parsing BacDive records for genus {genus}: {e}")
        pbar.update(1)
        sys.exit(1)
    
    for type_strain in type_strains:
        if type_strain.type_names is None:
            type_strain.type_names = [] # in these cases need to allow search based on binomial name
        for record in parsed_records:
            if set(type_strain.type_names) & set(record["type_names"]):
                try:
                    type_strain.type_names = list(set(type_strain.type_names) | set(record["type_names"]))
                    if type_strain.rRNA_acc is None and record.get("16S_rRNA") is not None:
                        type_strain.rRNA_acc = record.get("16S_rRNA").get("accession")
                    if record.get('ncbi_tax_id') is not None:
                        type_strain.species_ncbi_tax_id = [x[0] for x in record.get("ncbi_tax_id") if x[1].lower() == "species"]
                        type_strain.strain_ncbi_tax_id = [x[0] for x in record.get("ncbi_tax_id") if x[1].lower() == "strain"]
                    for i in type_strain_attributes:
                        if record.get(i) is not None:
                            setattr(type_strain, i, record.get(i))
                    
                except Exception as e:
                    logger.error(f"Error updating type strain {type_strain.parent_species} with BacDive data: {e}")
                    traceback.print_exc()
                break
        lpsn_types.append(type_strain)
    pbar.update(1)
    
    return lpsn_types


async def gather_all_genera(bacdive_client, type_strain_dict):
    pbar = tqdm(total=len(type_strain_dict), desc="Genus groups", unit="genera", ncols=100, colour="magenta")
    tasks = [
        fetch_genus_data(genus, type_strains, bacdive_client, pbar)
        for genus, type_strains in type_strain_dict.items()
    ]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    results = [type_strain for sublist in results for type_strain in sublist]
    await bacdive_client.close()
    pbar.close()
    return results



def get_genus_type_strains(bacdive_client, type_strain_dict):

    lpsn_types = []

    logger.info("Retrieving BacDive data for type strains in genus groups.")
    lpsn_types = asyncio.run(gather_all_genera(bacdive_client, type_strain_dict))

    return lpsn_types


def retrieve_extra_info_from_bacdive(lpsn_types, bacdive_client):

    logger.info("Step 2 of 4: Retrieving extra info from BacDive for LPSN type strains.")
    lpsn_types = get_genus_type_strains(bacdive_client, split_lpsn_types_list_to_genus_dict(lpsn_types))
    logger.success("Extra info retrieval from BacDive complete.\n")

    return lpsn_types