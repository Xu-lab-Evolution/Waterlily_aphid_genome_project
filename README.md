# Waterlily_aphid_genome_project
This repository archives codes that were used in the waterlily aphid genome assembly and gene annotation project.

Codes were used in the MOGON cluster (Mainz University). 

Explanation of Major Steps (in brief):

1. BRAKER1/2: Predict genes using RNA-seq reads (BRAKER1) and protein evidence (BRAKER2).
2. Cupcake Collapse: Merge redundant isoforms from Iso-Seq data.
3. GeneMarkS-T: Identify protein-coding regions in transcripts.
4. TSEBRA: Merge and filter BRAKER results for a refined gene set.
5. BUSCO: Assess completeness of gene models.
6. Miniprot Alignment: Align external protein references to the genome for validation.
7. PfamScan & DOGMA: Characterize protein domains and architectures.
8. Filtering: Remove incomplete, intronless, and TE-contaminated gene models.
9. Final GFF3: Format, rename, and extract final protein sequences for downstream functional annotation.