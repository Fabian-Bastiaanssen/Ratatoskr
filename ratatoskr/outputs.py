import polars as pl
from pathlib import Path
from ratatoskr.utils import make_dir

phenotypic_attrs = ["morphology_cell", "spore_formation", "oxygen_tolerance", "isolation_sample_type", 
                    "culture_medium", "culture_temp", "culture_pH", "morphology_colony", "morphology_pigmentation",  
                    "compound_production", "metabolite_production"]

genomic_attrs = ["rRNA_acc", "genome_acc"]


taxonomic_attrs = [
        "type_names", "authority", "parent_domain", "parent_kingdom", "parent_phylum", "parent_class", 
        "parent_order", "parent_family", "parent_genus", "parent_species", "parent_subspecies", 
        "binomial_synonyms", "parent_domain_id", "parent_kingdom_id", "parent_phylum_id", "parent_class_id",
        "parent_order_id", "parent_family_id", "parent_genus_id", "parent_species_id", "parent_subspecies_id",
        "species_ncbi_tax_id", "strain_ncbi_tax_id"
    ]


def get_genomic_data(type_strain, type_name):
    values = [type_name]
    for a in genomic_attrs:
        if a == "rRNA_acc":
            rRNA_info = getattr(type_strain, "rRNA_info") or {}
            values.append(
                    f"{type_strain.rRNA_acc if type_strain.rRNA_acc is not None else ''}\t"
                    f"{rRNA_info.get('source','') if rRNA_info.get('source') is not None else ''}\t"
                    f"{rRNA_info.get('97.2%_match','') if rRNA_info.get('97.2%_match') is not None else ''}\t"
                    f"{rRNA_info.get('90.1%_match','') if rRNA_info.get('90.1%_match') is not None else ''}\t"
                    f"{rRNA_info.get('80.1%_match','') if rRNA_info.get('80.1%_match') is not None else ''}\t"
                    f"{rRNA_info.get('72.9%_match','') if rRNA_info.get('72.9%_match') is not None else ''}\t"
                    f"{rRNA_info.get('note','') if rRNA_info.get('note') is not None else ''}"
                )
        elif a == "genome_acc":
            bg = getattr(type_strain, a) or {}
            values.append(
                    f"{bg.get('accession','') if bg.get('accession') is not None else ''}\t"
                    f"{bg.get('assembly level','') if bg.get('assembly level') is not None else ''}\t"
                    f"{bg.get('checkm_completeness','') if bg.get('checkm_completeness') is not None else ''}\t"
                    f"{bg.get('checkm_contamination','') if bg.get('checkm_contamination') is not None else (0 if bg.get('checkm_completeness') is not None and (bg.get('checkm_contamination') is None) else '')}\t"
                    f"{bg.get('ncbis_species_match_assembly','') if bg.get('ncbis_species_match_assembly') is not None else ''}\t"
                    f"{bg.get('ncbis_species_match_name','') if bg.get('ncbis_species_match_name') is not None else ''}\t"
                    f"{bg.get('ncbis_species_match_ani','') if bg.get('ncbis_species_match_ani') is not None else ''}\t"
                    f"{bg.get('ncbis_species_match_coverage','') if bg.get('ncbis_species_match_coverage') is not None else ''}\t"
                    f"{bg.get('ncbis_best_match_assembly','') if bg.get('ncbis_best_match_assembly') is not None else ''}\t"
                    f"{bg.get('ncbis_best_match_name','') if bg.get('ncbis_best_match_name') is not None else ''}\t"
                    f"{bg.get('ncbis_best_match_ani','') if bg.get('ncbis_best_match_ani') is not None else ''}\t"
                    f"{bg.get('ncbis_best_match_coverage','') if bg.get('ncbis_best_match_coverage') is not None else ''}"
                )

    return "\t".join(values)


def get_taxonomic_data(type_strain, type_name):
    values = [type_name]
    for a in taxonomic_attrs:
        if a == "binomial_synonyms":
            correct = type_strain.parent_subspecies if type_strain.parent_subspecies is not None else type_strain.parent_species
            synonyms = [x for x in type_strain.binomial_synonyms if correct not in x] if type_strain.binomial_synonyms is not None else []
            values.append(';'.join([str(x) for x in synonyms] or ""))
        else:
            if type(getattr(type_strain, a)) is list:
                values.append(';'.join([str(x) for x in getattr(type_strain, a)] or ""))
            else:
                values.append(str(getattr(type_strain, a) or ""))
    return "\t".join(values)


def write_tsv(data, output_file):
    with open(output_file, "w") as f:
        for row in data:
            f.write(row + "\n")


def format_colony_morph_dict(d):
    sections = []
    for medium, attrs in d.items():
        parts = []
        for key, values in attrs.items():
            cleaned = [v for v in values if v is not None]
            if cleaned:
                parts.append("/".join(cleaned))
            else:
                parts.append("")
        section_str = f"{medium}=({ '~'.join(parts) })"
        sections.append(section_str)
    return ";".join(sections)


def get_general_phenotypic_data(type_strain, type_name):
    values = [type_name]
    for a in phenotypic_attrs:
        attr_val = getattr(type_strain, a)

        if attr_val is None or attr_val == set():
            attr_val = None

        if a == "morphology_cell":
            if attr_val is None:
                cell_width = cell_length = gram_stain = cell_shape = motility = flagellum_arrangement = ""
            else:
                cell_data = attr_val
                cell_width = cell_data.get('cell_width') or ""
                cell_length = cell_data.get('cell_length') or ""
                gram_stain = cell_data.get('gram_stain') or ""
                cell_shape = cell_data.get('cell_shape') or ""
                motility = cell_data.get('motility') or ""
                flagellum_arrangement = cell_data.get('flagellum_arrangement') or ""
            values.append("\t".join([
                gram_stain,
                motility,
                cell_shape,
                cell_width,
                cell_length,
                flagellum_arrangement,
            ]))

        elif a == "morphology_colony":
            if attr_val is None:
                morphology_colony = ""
            else:
                morphology_colony = format_colony_morph_dict(attr_val)
            values.append(morphology_colony)

        elif a == "metabolite_production":
            if attr_val is None:
                values.append("")
            else:
                metabolites = [f"{k}={v}" for k, v in attr_val.items()]
                values.append(";".join(metabolites))

        else:
            values.append("" if attr_val is None else str(attr_val))

    return "\t".join(values)

def output_metabolite_utilization_data(lpsn_types, output_path):

    all_metabolites = set()

    for ts in lpsn_types:
        if ts.metabolite_utilization is not None:
            for metabolite in ts.metabolite_utilization.keys():
                all_metabolites.add(metabolite)

    all_dict = {x: [] for x in ['Name'] + sorted((list(all_metabolites)))}

    for ts in lpsn_types:
        type_name = " ".join([ts.parent_subspecies, ts.type_names[0]] if ts.parent_subspecies is not None else [ts.parent_species, ts.type_names[0]])
        all_dict['Name'].append(type_name)
        for metabolite in all_metabolites:
            if ts.metabolite_utilization is not None and metabolite in ts.metabolite_utilization:
                activity = '/'.join(set(ts.metabolite_utilization.get(metabolite).get('activity')))
                utilisation = '/'.join(set(ts.metabolite_utilization.get(metabolite).get('utilisation_type', 'ND')))
                all_dict[metabolite].append(f"{activity} ({utilisation})")
            else:
                all_dict[metabolite].append("")

    metabolite_utilization_tsv = pl.DataFrame(all_dict)
    metabolite_utilization_tsv.write_csv(output_path / "characteristics" / "metabolite_utilisation_metadata.tsv", separator="\t")
    

def output_fatty_acid_profile_data(lpsn_types, output_path):
    
    all_fatty_acids = set()

    for ts in lpsn_types:
        if ts.fatty_acid_profile is not None:
            for fatty_acid in ts.fatty_acid_profile.keys():
                all_fatty_acids.add(fatty_acid)

    all_dict = {x: [] for x in ['Name'] + sorted((list(all_fatty_acids)))}

    for ts in lpsn_types:
        if ts.type_names is None:
            ts.type_names = [""]
        if len(ts.type_names) == 0:
            ts.type_names = [""]
        type_name = " ".join([ts.parent_subspecies, ts.type_names[0]] if ts.parent_subspecies is not None else [ts.parent_species, ts.type_names[0]])
        all_dict['Name'].append(type_name)
        for fatty_acid in all_fatty_acids:
            if ts.fatty_acid_profile is not None and fatty_acid in ts.fatty_acid_profile:
                percentage = '/'.join(set(ts.fatty_acid_profile.get(fatty_acid)))
            else:
                percentage = ""
            all_dict[fatty_acid].append(percentage)

    fatty_acid_profile_tsv = pl.DataFrame(all_dict)
    fatty_acid_profile_tsv.write_csv(output_path / "characteristics" / "fatty_acid_profile_metadata.tsv", separator="\t")


def output_API_results_data(lpsn_types, output_path):

    make_dir( output_path / "characteristics" / "API_results" )

    apis = set()
    tests = {}

    for i in lpsn_types:
        if i.API_results is not None:
            for api, reactions in i.API_results.items():
                apis.add(api)
                if api not in tests:
                    tests[api] = set(i.API_results.get(api).keys())
                else:
                    tests[api].update(set(i.API_results.get(api).keys()))

    
    all_dfs = []

    for api in apis:
        all_dict = {api: {x: [] for x in ['Name'] + sorted(tests[api])}}
        for ts in lpsn_types:
            if ts.type_names is None:
                ts.type_names = [""]
            if len(ts.type_names) == 0:
                ts.type_names = [""]
            type_name = " ".join([ts.parent_subspecies, ts.type_names[0]] if ts.parent_subspecies is not None else [ts.parent_species, ts.type_names[0]])
            all_dict[api]['Name'].append(type_name)
            for test in tests[api]:
                if ts.API_results is not None and ts.API_results.get(api, {}).get(test) is not None:
                    result = '|'.join(set(ts.API_results.get(api).get(test)))
                    all_dict[api][test].append(result)
                else:
                    all_dict[api][test].append("")
        all_dfs.append(all_dict)

    for x in all_dfs:
        api = list(x.keys())[0]

        df = pl.DataFrame(x[api])
        df.write_csv(output_path / "characteristics" / "API_results" / f"{api}_results.tsv", separator = "\t")


def output_taxonomy(lpsn_types, output_path):

    taxonomic_tsv = ["\t".join(["Name", "Type_names", "Authority", "Domain", "Kingdom", "Phylum", "Class", 
                     "Order", "Family", "Genus", "Species", "Subspecies", "Synonyms", "Domain_LPSN_ID", 
                     "Kingdom_LPSN_ID", "Phylum_LPSN_ID", "Class_LPSN_ID", "Order_LPSN_ID", "Family_LPSN_ID", 
                     "Genus_LPSN_ID", "Species_LPSN_ID", "Subspecies_LPSN_ID", "species_ncbi_tax_id", 
                     "strain_ncbi_tax_id"])]
    
    for ts in lpsn_types:
        if ts.type_names is None:
            ts.type_names = [""]
        if len(ts.type_names) == 0:
            ts.type_names = [""]
        type_name = " ".join([ts.parent_subspecies, ts.type_names[0]] if ts.parent_subspecies is not None else [ts.parent_species, ts.type_names[0]])
        taxonomic_tsv.append(get_taxonomic_data(ts, type_name))
    
    write_tsv(taxonomic_tsv, output_path / "taxonomy" / "taxonomic_metadata.tsv")
    
    
def output_general_characteristics(lpsn_types, output_path):

    general_phenotypic_tsv = ["\t".join(["Name", "Gram_stain", "Motility", "Cell_shape", "Cell_width", "Cell_length", 
                                         "Flagellum_arrangement", "Spore_formation", "Oxygen_tolerance", "Isolation_sample", 
                                         "Culture_medium", "Culture_temp", "Culture_pH", "Morphology_colony", "Morphology_pigmentation",
                                         "Compound_production", "Metabolite_production", "Metabolite_utilization"])]
    
    for ts in lpsn_types:
        if ts.type_names is None:
            ts.type_names = [""]
        if len(ts.type_names) == 0:
            ts.type_names = [""]
        type_name = " ".join([ts.parent_subspecies, ts.type_names[0]] if ts.parent_subspecies is not None else [ts.parent_species, ts.type_names[0]])
        general_phenotypic_tsv.append(get_general_phenotypic_data(ts, type_name))
    
    write_tsv(general_phenotypic_tsv, output_path / "characteristics" / "general_characteristics.tsv")


def output_sequence_metadata(lpsn_types, output_path):
    make_dir( output_path / "sequences" )
    genomic_tsv = ["\t".join(["Name", 
                              "rRNA_accession", 
                              "rRNA_source",
                              "rRNA_97.2%_match_status",
                              "rRNA_90.1%_match_status",
                              "rRNA_80.1%_match_status",
                              "rRNA_72.9%_match_status",
                              "rRNA_note",
                              "Genome_accession",
                              "Assembly_level", 
                              "CheckM_completeness", 
                              "CheckM_contamination",
                              "NCBIs_species_match_accession",
                              "NCBIs_species_match_name",
                              "NCBIs_species_match_ANI",
                              "NCBIs_species_match_assembly_coverage",
                              "NCBIs_best_match_accession",
                              "NCBIs_best_match_name",
                              "NCBIs_best_match_ANI",
                              "NCBIs_best_match_assembly_coverage",])]
    
    
    for ts in lpsn_types:
        type_name = " ".join([ts.parent_subspecies, ts.type_names[0]] if ts.parent_subspecies is not None else [ts.parent_species, ts.type_names[0]])
        genomic_tsv.append(get_genomic_data(ts, type_name))

    write_tsv(genomic_tsv, output_path / "sequences" / "sequence_metadata.tsv")


def output_metadata(lpsn_types, output_path):

    lpsn_types = sorted(lpsn_types, key=lambda x: (x.parent_species or "", x.parent_subspecies or ""))

    make_dir( output_path / "characteristics" )
    make_dir( output_path / "taxonomy" )

    output_taxonomy(lpsn_types, output_path)
    output_general_characteristics(lpsn_types, output_path)
    output_sequence_metadata(lpsn_types, output_path)
    output_metabolite_utilization_data(lpsn_types, output_path)
    output_fatty_acid_profile_data(lpsn_types, output_path)
    output_API_results_data(lpsn_types, output_path)
    
