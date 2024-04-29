import time
import pandas as pd

from evepidr.variants.uniprot_proteins import *
from evepidr.variants.clinvar_variants import *
from evepidr.variants.add_labels_to_variants import label_if_substitution_in_idr


uniprot_ids = ["O00571", "Q06787", "P17600", "P35222", "Q00839", "P15884", "O60741", "O60741", "Q02548", "Q9H165", "P15056", "Q9Y6K1", "Q9H334", "P49711", "P60484", "Q06124"]

gene_to_sequence, gene_to_uniprot_id = get_canonical_sequence_from_uniprot(uniprot_ids)

save_sequences_as_fasta(gene_to_sequence, "report/data/asd_linked_idps.fasta")

clinvar_ids = []
for gene in gene_to_uniprot_id.keys():
    clinvar_ids += clinvar_snp_missense_variants_id_list(gene)
    time.sleep(0.333333334)

clinvar_xml = clinvar_variant_info(clinvar_ids)
clinvar_xml.write("report/data/clincar_data.xml")

clinvar_df = clean_clinvar_xml_variants(gene_to_uniprot_id, clinvar_xml)

# List of IDRs sourced from MobiDB
idrs = {
    'BCL11A': [(17,49), (68, 173), (192, 382), (421, 741), (821, 835)],
    'BRAF': [(1, 43), (102, 157), (274, 454), (721, 766)],
    'CTNNB1': [(1, 222), (549, 562), (662, 781)],
    'DDX3X': [(1, 167), (576, 662)],
    'DNMT3A': [(1, 39), (51, 183), (198, 282), (385, 397), (431, 473), (610, 620), (832, 847)],
    'FMR1': [(200, 213), (281, 429), (437, 632)],
    'FOXP1': [(1, 130), (154, 179), (198, 300), (357, 447), (549, 677)],
    'HCN1': [(1, 93), (635, 890)],
    'HNRNPU': [(41, 281), (671, 825)],
    'PAX5': [(1, 83), (143, 391)],
    'SYN1': [(1, 111), (418, 705)],
    'TCF4': [(1, 10), (28, 569), (626, 667)],
    'PTEN': [(286, 309), (352, 403)],
    'PTPN11': [(535, 593)],
    'CTCF': [(1, 269), (572, 727)]
}

clinvar_df = label_if_substitution_in_idr(clinvar_df, idrs)
clinvar_df.to_csv("report/data/clinvar_data_patho_benign.csv", index=False)
