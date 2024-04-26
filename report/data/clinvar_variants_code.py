import time
import pandas as pd
from evepidr.variants.clinvar_variants import *

uniprot_ids = ["O00571", "Q06787", "P17600", "P35222", "Q00839", "Q96PV0", "P15884", "O60741", "O60741", "Q02548", "Q9H165", "P15056", "Q9Y6K1", "Q9H334"]

genes_to_sequence = get_canonical_sequence_from_uniprot(uniprot_ids)

clinvar_ids = []
for gene in genes_to_sequence.keys():
    clinvar_ids += clinvar_snp_missense_variants_id_list(gene)
    time.sleep(0.333333334)

clinvar_xml = clinvar_variant_info(clinvar_ids)
clinvar_xml.write("ClinVar_data_all.xml")

clinvar_df = clean_clinvar_xml_variants(genes_to_sequence, clinvar_xml)
clinvar_df.to_csv("ClinVar_data_patho_benign.csv", index=False)
