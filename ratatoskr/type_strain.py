from dataclasses import dataclass, replace
from typing import Optional



@dataclass
class TypeStrain:
    parent_domain:    Optional[str] = None
    parent_domain_id:    Optional[int] = None
    parent_kingdom:   Optional[str] = None
    parent_kingdom_id:   Optional[int] = None
    parent_phylum:    Optional[str] = None
    parent_phylum_id:    Optional[int] = None
    parent_class:     Optional[str] = None
    parent_class_id:     Optional[int] = None
    parent_subclass: Optional[str] = None
    parent_subclass_id: Optional[int] = None
    parent_order:     Optional[str] = None
    parent_order_id:     Optional[int] = None
    parent_suborder: Optional[str] = None
    parent_suborder_id: Optional[int] = None
    parent_family:    Optional[str] = None
    parent_family_id:    Optional[int] = None
    parent_genus:     Optional[str] = None
    parent_genus_id:     Optional[int] = None
    parent_species:   Optional[str] = None
    parent_species_id:   Optional[int] = None
    parent_subspecies: Optional[str] = None
    parent_subspecies_id: Optional[int] = None

    last_id:          Optional[int] = None 
    last_level:       Optional[str] = None

    authority:        Optional[str] = None

    binomial_synonyms: list = None
    lpsn_id:           Optional[str] = None
    rRNA_acc:           Optional[str] = None
    rRNA_info:          Optional[str] = None
    list_ref:          Optional[str] = None
    pub_ref:           Optional[int] = None
    type_names:        list = None
    authority: Optional[str] = None,
    species_ncbi_tax_id:      Optional[int] = None
    strain_ncbi_tax_id:       Optional[int] = None
    genome_acc:        Optional[str] = None

    isolation_sample_type: list = None
    culture_pH:        Optional[str] = None
    culture_temp:      Optional[str] = None
    culture_medium:    Optional[str] = None
    morphology_pigmentation: Optional[str] = None
    morphology_cell: Optional[str] = None
    morphology_colony: Optional[str] = None
    API_results: Optional[str] = None
    fatty_acid_profile: Optional[str] = None
    oxygen_tolerance: Optional[str] = None
    spore_formation: Optional[str] = None
    compound_production: Optional[str] = None
    metabolite_production: Optional[str] = None
    metabolite_utilization: Optional[str] = None

    def copy(self) -> "TypeStrain":
        """Return a shallow copy of this instance."""
        return replace(self)

    def set_parent(self, rank: str, name: str, parent_id: Optional[int] = None) -> "TypeStrain":
        attr = f"parent_{rank}"
        attr_id = f"{attr}_id"
        if not hasattr(self, attr) or not hasattr(self, attr_id):
            raise ValueError(f"Unsupported rank: {rank!r}")
        setattr(self, attr, name)
        if parent_id is not None:
            setattr(self, attr_id, parent_id)
        return self 


    def set_last(self, level: str, last_id: Optional[int] = None) -> "TypeStrain":
        self.last_level = level
        if last_id is not None:
            self.last_id = last_id
        return self 

    def set_metadata(self, 
                     rRNA_acc: Optional[str] = None, 
                     list_ref: Optional[str] = None, 
                     pub_ref: Optional[int] = None, 
                     type_names: list = None,
                     authority: Optional[str] = None,
                     binomial_synonyms: list = None) -> "TypeStrain":
        
        if rRNA_acc is not None:
            self.rRNA_acc = rRNA_acc
        if list_ref is not None:
            self.list_ref = list_ref
        if pub_ref is not None:
            self.pub_ref = pub_ref
        if type_names is not None:
            self.type_names = type_names
        if authority is not None:
            self.authority = authority
        if binomial_synonyms is not None:
            self.binomial_synonyms = binomial_synonyms
        
        return self
        # Could add more metadata fields here as needed

    def __repr__(self) -> str:
        return (
            "ReferenceStrain("
            f"domain={self.parent_domain!r}, "
            f"domain_id={self.parent_domain_id!r}, "
            f"kingdom={self.parent_kingdom!r}, "
            f"kingdom_id={self.parent_kingdom_id!r}, "
            f"phylum={self.parent_phylum!r}, "
            f"phylum_id={self.parent_phylum_id!r}, "
            f"class={self.parent_class!r}, "
            f"class_id={self.parent_class_id!r}, "
            f"subclass={self.parent_subclass!r}, "
            f"subclass_id={self.parent_subclass_id!r}, "
            f"order={self.parent_order!r}, "
            f"order_id={self.parent_order_id!r}, "
            f"suborder={self.parent_suborder!r}, "
            f"suborder_id={self.parent_suborder_id!r}, "
            f"family={self.parent_family!r}, "
            f"family_id={self.parent_family_id!r}, "
            f"genus={self.parent_genus!r}, "
            f"genus_id={self.parent_genus_id!r}, "
            f"species={self.parent_species!r}, "
            f"species_id={self.parent_species_id!r}, "
            f"subspecies={self.parent_subspecies!r}, "
            f"subspecies_id={self.parent_subspecies_id!r}, "
            f"rRNA_acc={self.rRNA_acc!r}, "
            f"rRNA_info={self.rRNA_info!r}, "
            f"list_ref={self.list_ref!r}, "
            f"pub_ref={self.pub_ref!r}, "
            f"type_names={self.type_names!r}, "
            f"authority={self.authority!r}, "
            f"binomial_synonyms={self.binomial_synonyms!r},"
            f"species_ncbi_tax_id={self.species_ncbi_tax_id!r}, "
            f"strain_ncbi_tax_id={self.strain_ncbi_tax_id!r}, "
            f"genome_acc={self.genome_acc!r}, "
            f"culture_pH={self.culture_pH!r}, "
            f"culture_temp={self.culture_temp!r}, "
            f"culture_medium={self.culture_medium!r}, "
            f"morphology_pigmentation={self.morphology_pigmentation!r}, "
            f"morphology_cell={self.morphology_cell!r}, "
            f"morphology_colony={self.morphology_colony!r}, "
            f"API_results={self.API_results!r}, "
            f"fatty_acid_profile={self.fatty_acid_profile!r}, "
            f"oxygen_tolerance={self.oxygen_tolerance!r}, "
            f"spore_formation={self.spore_formation!r}, "
            f"compound_production={self.compound_production!r}, "
            f"metabolite_production={self.metabolite_production!r}, "
            f"metabolite_utilization={self.metabolite_utilization!r}, "
            ")"
        )