######################################################################################
######################################################################################
######################################################################################
######################################################################################
# Annotation Pipeline

######################################################## Running BRAKER1
#!/bin/bash
# Script for BRAKER1 gene prediction using short RNA-seq reads

#SBATCH --export=ALL               # Start with a clean environment
#SBATCH -A m2_jgu-EvolTroph        # Specify project/account
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --partition parallel       # Partition to submit to
#SBATCH --job-name=braker1         # Job name
#SBATCH --output=/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/0.braker1/braker1.out
#SBATCH --error=/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/0.braker1/braker1.err
#SBATCH --time=24:00:00            # Walltime
#SBATCH --cpus-per-task=24         # Number of CPU cores
#SBATCH --mem=16G                  # Memory per node

ml bio/BRAKER/2.1.5-foss-2020b-Python-3.8.6  # Load BRAKER module

# Define paths
braker_dir="/home/ywang/software/long_reads_anno/BRAKER_long_reads_branch/BRAKER/scripts/"
shortRNA_bam="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/RNA-seq/clean_Reads/Rn_RNA-seq.sort.bam"
iso_seq_bam="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/iso/Rn_0228.iso_seq.sort.bam"
uni_prot_fa="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/uniprot_protein_fa/uniprot-compressed_true_download_true_format_fasta_query__28Aphidida-2023.05.11-14.59.28.21.fasta.name_fixed.fa"
fa_soft_mask="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/repeat_anno.per_chr/per_chr_fa.dir/merge_masked_res/soft.masked.fa_dir/Rn_01.asm.bp.hap1.p_ctg.soft_masked.fasta"
work_dir="/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/"

# Uncomment this section if you'd like to run the first attempt. 
# $braker_dir/braker.pl \
#   --genome=$fa_soft_mask \
#   --useexisting \
#   --species=Rn_01 \
#   --softmasking --cores=24 \
#   --bam=$shortRNA_bam  \
#   --workingdir=$work_dir/0.braker1/ \
#   --AUGUSTUS_CONFIG_PATH=/home/ywang/software/long_reads_anno/Augustus/config/ \
#   --CDBTOOLS_PATH=/gpfs/fs1/home/ywang/software/long_reads_anno/cdbfasta/ \
#   --GENEMARK_PATH=/home/ywang/software/genemark-ex/gmes_linux_64/ \
#   --DIAMOND_PATH=/home/ywang/software/long_reads_anno/diamond/ \
#   --PROTHINT_PATH=/home/ywang/software/long_reads_anno/ProtHint/bin/ \
#   2> $work_dir/0.braker1/braker1.log

# Actual BRAKER1 command
$braker_dir/braker.pl \
  --genome=$fa_soft_mask \
  --useexisting \
  --species=Rn_01 \
  --softmasking --cores=24 \
  --bam=$shortRNA_bam  \
  --workingdir=$work_dir/0.braker1/ \
  --AUGUSTUS_CONFIG_PATH=/home/ywang/software/long_reads_anno/Augustus/config/ \
  --AUGUSTUS_SCRIPTS_PATH=/home/ywang/software/long_reads_anno/Augustus/scripts/ \
  --AUGUSTUS_BIN_PATH=/home/ywang/software/long_reads_anno/Augustus/bin/ \
  --CDBTOOLS_PATH=/gpfs/fs1/home/ywang/software/long_reads_anno/cdbfasta/ \
  --DIAMOND_PATH=/home/ywang/software/long_reads_anno/diamond/ \
  --PROTHINT_PATH=/home/ywang/software/long_reads_anno/ProtHint/bin/ \
  --GENEMARK_PATH=/home/ywang/software/long_reads_anno/GeneMark-ETP/bin/gmes/ \
  2> $work_dir/0.braker1/braker1.log


######################################################## Running BRAKER2
#!/bin/bash
# Script for BRAKER2 gene prediction using external protein evidence

#SBATCH --export=ALL               # Start with a clean environment
#SBATCH -A m2_jgu-EvolTroph        # Specify project/account
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --partition parallel       # Partition to submit to
#SBATCH --job-name=braker2         # Job name
#SBATCH --output=/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/1.braker2/braker2.out
#SBATCH --error=/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/1.braker2/braker2.err
#SBATCH --time=24:00:00            # Walltime
#SBATCH --cpus-per-task=24         # Number of CPU cores
#SBATCH --mem=24G                  # Memory per node

ml bio/BRAKER/2.1.5-foss-2020b-Python-3.8.6  # Load BRAKER module

# Define paths
braker_dir="/home/ywang/software/long_reads_anno/BRAKER_long_reads_branch/BRAKER/scripts/"
shortRNA_bam="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/RNA-seq/clean_Reads/Rn_RNA-seq.sort.bam"
iso_seq_bam="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/iso/Rn_0228.iso_seq.sort.bam"
uni_prot_fa="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/uniprot_protein_fa/uniprot-compressed_true_download_true_format_fasta_query__28Aphidida-2023.05.11-14.59.28.21.fasta.name_fixed.fa"
fa_soft_mask="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/repeat_anno.per_chr/per_chr_fa.dir/merge_masked_res/soft.masked.fa_dir/Rn_01.asm.bp.hap1.p_ctg.soft_masked.fasta"
work_dir="/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/"

# BRAKER2 command with protein sequences
$braker_dir/braker.pl \
  --genome=$fa_soft_mask \
  --useexisting \
  --species=Rn_01 \
  --softmasking --cores=24 \
  --epmode --prot_seq=$uni_prot_fa  \
  --workingdir=$work_dir/1.braker2/ \
  --AUGUSTUS_CONFIG_PATH=/home/ywang/software/long_reads_anno/Augustus/config/ \
  --AUGUSTUS_SCRIPTS_PATH=/home/ywang/software/long_reads_anno/Augustus/scripts/ \
  --AUGUSTUS_BIN_PATH=/home/ywang/software/long_reads_anno/Augustus/bin/ \
  --CDBTOOLS_PATH=/gpfs/fs1/home/ywang/software/long_reads_anno/cdbfasta/ \
  --DIAMOND_PATH=/home/ywang/software/long_reads_anno/diamond/ \
  --PROTHINT_PATH=/home/ywang/software/long_reads_anno/ProtHint/bin/ \
  --GENEMARK_PATH=/home/ywang/software/long_reads_anno/GeneMark-ETP/bin/gmes/ \
  2> $work_dir/1.braker2/braker2.log


######################### Cupcake collapse
#!/bin/bash
# Script to collapse isoforms using Cupcake

#SBATCH --export=ALL               # Start with a clean environment
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --partition parallel       # Partition to submit to
#SBATCH --job-name=Cupcake_collapse
#SBATCH --output=/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/2.GeneMarkS-T/cupcake_collapse/cupcake_collapse.out
#SBATCH --error=/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/2.GeneMarkS-T/cupcake_collapse/cupcake_collapse.err
#SBATCH --time=4:00:00             # Walltime
#SBATCH --cpus-per-task=8          # Number of CPU cores
#SBATCH --mem=12G                  # Memory per node

ml bio/BRAKER/2.1.5-foss-2020b-Python-3.8.6  # Load necessary modules

# Define Iso-Seq directory
isoseq_dir="/lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/iso"

# Convert BAM to FASTQ if needed (commented out)
# samtools fastq -@ 4 Rn_0228.iso_seq.sort.bam > Rn_0228.iso_seq.sort.fastq

# Collapse isoforms
collapse_isoforms_by_sam.py --input $isoseq_dir/Rn_0228.iso_seq.sort.fastq --fq \
 -b $isoseq_dir/Rn_0228.iso_seq.sort.bam --dun-merge-5-shorter \
 -o Rn_0228.iso_seq.Cupcake.collapse --cpus 8

################## Run GeneMarkS-T to predict protein-coding regions in transcripts
# Extract fasta sequences from the reference using transcripts GFF
stringtie2fa.py -g /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/Rn_01.asm.bp.hap1.p_ctg.fasta \
 -f Rn_0228.iso_seq.Cupcake.collapse.collapsed.gff \
 -o Rn_0228.iso_seq.Cupcake.collapse.collapsed.fa

# Run GeneMarkS-T on the transcript fasta
gmst.pl --strand direct Rn_0228.iso_seq.Cupcake.collapse.collapsed.fa.mrna \
 --output gmst.out --format GFF

################## Convert GeneMarkS-T coordinates to genomic global coordinates
/home/ywang/software/long_reads_anno/BRAKER_long_reads_branch/BRAKER/scripts/gmst2globalCoords.py \
 -t Rn_0228.iso_seq.Cupcake.collapse.collapsed.gff \
 -p gmst.out \
 -o gmst.global.gtf \
 -g /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/Rn_01.asm.bp.hap1.p_ctg.fasta

################# Run TSEBRA to select best gene models from BRAKER outputs
/home/ywang/software/long_reads_anno/TSEBRA_long_reads_branch/TSEBRA/bin/tsebra.py \
 -g braker1.augustus.hints.gtf,braker2.augustus.hints.gtf,braker3.augustus.hints.gtf \
 -e braker1.hintsfile.gff,braker2.hintsfile.gff,braker3.hintsfile.gff \
 -l gmst.global.gtf \
 -c /gpfs/fs1/home/ywang/software/long_reads_anno/TSEBRA_long_reads_branch/TSEBRA/config/long_reads.cfg \
 -o tsebra.gtf --score_tab tsebra.gtf.score

# Example TSEBRA usage with an additional GeneMark file (commented out)
# /home/ywang/software/long_reads_anno/TSEBRA/bin/tsebra.py --gtf \
#    /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/3.braker3/augustus.hints.gtf,\
#    /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/3.braker3/GeneMark-ETP/genemark.gtf \
#    --keep_gtf /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/3.braker3/GeneMark-ETP/training.gtf \
#    --hintfiles /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/3.braker3/hintsfile.gff \
#    --filter_single_exon_genes \
#    --cfg /home/ywang/software/long_reads_anno/TSEBRA/bin/../config/braker3.cfg \
#    --out /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/3.braker3/braker.gtf -q 2>...

# TSEBRA filter single-exon genes (commented out)
# ~/software/long_reads_anno/TSEBRA/bin/tsebra.py -g braker1.augustus.hints.gtf,braker2.augustus.hints.gtf,braker3.augustus.hints.gtf \
#  -e braker1.hintsfile.gff,braker2.hintsfile.gff,braker3.hintsfile.gff \
#  -l gmst.global.gtf \
#  --filter_single_exon_genes \
#  -c /gpfs/fs1/home/ywang/software/long_reads_anno/TSEBRA_long_reads_branch/TSEBRA/config/long_reads.cfg \
#  -o tsebra.filter_single_exon.gtf


########################################################### BUSCO evaluation
#!/bin/bash
# BUSCO assessment of TSEBRA-predicted proteins

#SBATCH --export=ALL               # Start with a clean environment
#SBATCH -A m2_jgu-EvolTroph        # Specify project/account
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --cpus-per-task=24         # Number of CPU cores
#SBATCH --mem=8G                   # Memory per node
#SBATCH --partition parallel       # Partition
#SBATCH --time=02:00:00            # Walltime
#SBATCH --job-name=busco.tsebra.aa
#SBATCH -o /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/busco.tsebra.aa.out
#SBATCH -e /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/busco.tsebra.aa.err

ml bio/BUSCO/5.4.3-foss-2021b  # Load BUSCO module

# Run BUSCO on predicted protein set (tsebra.aa.fa)
busco -f -m proteins \
 -i /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/tsebra.aa.fa \
 -o Rn_01.tsebra.aa \
 -l hemiptera_odb10 \
 --out_path /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval \
 --offline \
 --download_path /lustre/project/m2_jgu-evoltroph/ywang/busco_download/v5/data \
 -c 24

#!/bin/bash
# BUSCO for TSEBRA proteins after TE filtering

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=8G
#SBATCH --partition parallel
#SBATCH --time=02:00:00
#SBATCH --job-name=busco.tsebra.filter_TEs.aa
#SBATCH -o /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/busco.tsebra.filter_TEs.aa.out
#SBATCH -e /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/busco.tsebra.filter_TEs.aa.err

ml bio/BUSCO/5.4.3-foss-2021b  # Load BUSCO module

# Run BUSCO on TSEBRA proteins (TE-filtered)
busco -f -m proteins \
 -i /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/tsebra.filter_TEs.aa.fa \
 -o RRn_01.tsebra.filter_TEs.aa \
 -l hemiptera_odb10 \
 --out_path /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval \
 --offline \
 --download_path /lustre/project/m2_jgu-evoltroph/ywang/busco_download/v5/data \
 -c 24

#!/bin/bash
# BUSCO run on genome assembly for comparison

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=8G
#SBATCH --partition parallel
#SBATCH --time=02:00:00
#SBATCH --job-name=busco
#SBATCH -o /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/busco.out
#SBATCH -e /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval/busco.err

ml bio/BUSCO/5.4.3-foss-2021b  # Load BUSCO module

# Run BUSCO on the raw genome assembly
busco -f -m proteins \
 -i /lustre/project/m2_jgu-evolpest/xus/Rhopalosiphum_nymphaeae/02.annotation/Rn_01.asm.bp.hap1.p_ctg.fasta \
 -o Rn_01.busco \
 -l hemiptera_odb10 \
 --out_path /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/busco_eval \
 --offline \
 --download_path /lustre/project/m2_jgu-evoltroph/ywang/busco_download/v5/data \
 -c 24


##################################################### Additional evaluation steps
cd /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation

# Extract protein sequences from GFF annotation
agat_sp_extract_sequences.pl -p \
 -g tsebra.longest_isoform.gff \
 -f /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/repeat_anno.per_chr/per_chr_fa.dir/merge_masked_res/soft.masked.fa_dir/Rn_01.asm.bp.hap1.p_ctg.soft_masked.fasta \
 -o tsebra.longest_isoform.aa.fa &

agat_sp_extract_sequences.pl -p \
 -g tsebra.filter.gff \
 -f /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/repeat_anno.per_chr/per_chr_fa.dir/merge_masked_res/soft.masked.fa_dir/Rn_01.asm.bp.hap1.p_ctg.soft_masked.fasta \
 -o tsebra.filter.aa.fa &

# Evaluate BUSCO on tsebra.longest_isoform.aa.fa
#! /bin/bash
# busco.tsebra.aa

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=8G
#SBATCH --partition parallel
#SBATCH --time=02:00:00
#SBATCH --job-name=busco.tsebra.aa
#SBATCH -o /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.longest_isoform.aa_busco/busco.tsebra.aa.out
#SBATCH -e /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.longest_isoform.aa_busco/busco.tsebra.aa.err

ml bio/BUSCO/5.4.3-foss-2021b

busco -f -m proteins \
 -i /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.longest_isoform.aa.fa \
 -o Rn_01.tsebra.longest_isoform.aa \
 -l hemiptera_odb10 \
 --out_path /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.longest_isoform.aa_busco \
 --offline \
 --download_path /lustre/project/m2_jgu-evoltroph/ywang/busco_download/v5/data \
 -c 24

# Evaluate BUSCO on tsebra.filter.aa.fa
#! /bin/bash
# busco.tsebra.filter.aa

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=8G
#SBATCH --partition parallel
#SBATCH --time=02:00:00
#SBATCH --job-name=busco.tsebra.filter_TEs.aa
#SBATCH -o /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.filter.aa_busco/busco.tsebra.filter.aa.out
#SBATCH -e /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.filter.aa_busco/busco.tsebra.filter.aa.err

ml bio/BUSCO/5.4.3-foss-2021b

busco -f -m proteins \
 -i /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.filter.aa.fa \
 -o Rn_01.tsebra.filter.aa \
 -l hemiptera_odb10 \
 --out_path /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/tsebra.filter.aa_busco \
 --offline \
 --download_path /lustre/project/m2_jgu-evoltroph/ywang/busco_download/v5/data \
 -c 24


############################################################ Evaluation using genome browser
# Align protein sequences to the reference genome using miniprot

#! /bin/bash
# miniprot_index: build miniprot index

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --partition parallel
#SBATCH --time=02:00:00
#SBATCH --job-name=miniprot_index
#SBATCH -o /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/miniprot_index.out
#SBATCH -e /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/miniprot_index.err

# Build index for the reference genome
miniprot -t24 -d /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/Rn_01.asm.bp.hap1.p_ctg.fasta.mpi \
 /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/Rn_01.asm.bp.hap1.p_ctg.fasta

#! /bin/bash
# miniprot_align: align external protein sequences with the reference

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --partition parallel
#SBATCH --time=02:00:00
#SBATCH --job-name=miniprot_align
#SBATCH -o /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/miniprot_align.out
#SBATCH -e /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/5.evaluation/miniprot_align.err

# Run protein alignment with miniprot
miniprot -t24 /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/Rn_01.asm.bp.hap1.p_ctg.fasta.mpi \
 /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/uniprot_protein_fa/uniprot-compressed_true_download_true_format_fasta_query__28Aphidida-2023.05.11-14.59.28.21.fasta.name_fixed.fa \
 --gff > protein_align.gff


###########################################################################################
################################## DOGMA evaluation #######################################
# Use PfamScan to identify protein domains and run DOGMA to classify domain architectures

#! /bin/bash
# PfamScan

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --partition parallel
#SBATCH --time=3:00:00
#SBATCH --job-name=PfamScan
#SBATCH -o ./PfamScan.out
#SBATCH -e ./PfamScan.err

ml lang/Perl/5.26.1-foss-2017a  # Load Perl environment

# Run PfamScan on protein fasta
/gpfs/fs1/home/ywang/software/PfamScan/pfam_scan.pl \
 -cpu 32 \
 -fasta Rn01.aa.fa \
 -dir /gpfs/fs1/home/ywang/software/PfamScan/Pfam_data \
 -out Rn01.pfam_scan

# Run DOGMA on PfamScan results
/home/ywang/software/bin/dogma.py proteome \
 -a Rn01.pfam_scan \
 -r insects \
 -o Rn01.DOGMA


########################################################
########################################################
########################################################
######################################################## Running TSEBRA on BRAKER1 & BRAKER2 only
# TSEBRA first run (using only braker1, braker2 results)
# cd /lustre/project/m2_jgu-evoltroph/ywang/03.annotation_insect/braker_run/4.tsebra/bk12
ml bio/BRAKER/2.1.5-foss-2020b-Python-3.8.6

/home/ywang/software/long_reads_anno/TSEBRA_long_reads_branch/TSEBRA/bin/tsebra.py \
 -g ../braker1.augustus.hints.gtf,../braker2.augustus.hints.gtf \
 -e ../braker1.hintsfile.gff,../braker2.hintsfile.gff \
 -l ../gmst.global.gtf \
 -c /gpfs/fs1/home/ywang/software/long_reads_anno/TSEBRA_long_reads_branch/TSEBRA/config/long_reads.cfg \
 -o tsebra.gtf


###################################################### Filtration Steps
############################## Reformat GTF to GFF
ml bio/AGAT/1.1.0-GCC-11.2.0

# Convert GTF to GFF
agat_convert_sp_gxf2gxf.pl -g tsebra.gtf -o tsebra.gff 2> cvt_gtf2gff.log

# Rename UTR features to the standard format
sed -i "s/5'-UTR/five_prime_UTR/g"  tsebra.gff
sed -i "s/3'-UTR/three_prime_UTR/g"  tsebra.gff

# Keep only the longest isoform per gene
agat_sp_keep_longest_isoform.pl -gff tsebra.gff -o tsebra.longest_isoform.gff

############################### Remove incomplete single-exon gene models
fa_soft_mask="/lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/repeat_anno.per_chr/per_chr_fa.dir/merge_masked_res/soft.masked.fa_dir/Rn_01.asm.bp.hap1.p_ctg.soft_masked.fasta"

# Tag incomplete coding models
agat_sp_filter_incomplete_gene_coding_models.pl \
 --gff tsebra.longest_isoform.gff \
 --fasta $fa_soft_mask --ad \
 -o tsebra.longest_isoform.start_stop_codon_tag.gff

# Filter out intronless genes
agat_sp_filter_gene_by_intron_numbers.pl \
 --gff tsebra.longest_isoform.start_stop_codon_tag.gff \
 -o tsebra.longest_isoform.start_stop_codon_tag.intronless.gff

# List incomplete and intronless genes
grep "incomplete" tsebra.longest_isoform.start_stop_codon_tag.intronless.gff | grep -oP "Parent=[^;]+" | sed 's/Parent=//' > intronless_and_incomplete_gene.list

# Remove those genes from the annotation
agat_sp_filter_feature_from_kill_list.pl \
 --gff tsebra.longest_isoform.gff \
 --kill_list intronless_and_incomplete_gene.list \
 -o tsebra.longest_isoform.rm_intronless_incomp.gff

############################ Extract merged CDS sequences for TE analysis
agat_sp_extract_sequences.pl -t cds --merge \
 -g tsebra.longest_isoform.rm_intronless_incomp.gff \
 -f $fa_soft_mask \
 -o tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa

# Run RepeatMasker on extracted CDS
#! /bin/bash
# RepeatMasker on transcripts

#SBATCH --export=ALL
#SBATCH -A m2_jgu-EvolTroph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --partition parallel
#SBATCH --time=1:00:00
#SBATCH --job-name=rm_on_transcripts.merge_cds
#SBATCH -o ./rm_on_transcripts.merge_cds.out
#SBATCH -e ./rm_on_transcripts.merge_cds.err

ml bio/RepeatMasker/4.1.2-p1-foss-2020b

/gpfs/fs1/home/ywang/software/repeatMasker/RepeatMasker/RepeatMasker \
 -e rmblast -xsmall -pa 32 -gff \
 -lib /lustre/miifs01/project/m2_jgu-evoltroph/ywang/03.annotation_insect/ref/repeat_anno.per_chr/merged.consensi.fa.classified \
 ./tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa

# Check how many sequences are masked
grep -v "#" tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa.out.gff | awk '{print $1}' | sort | uniq | wc -l

# Calculate TE proportion for each gene
calc_TE_proportion.pl tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa.out.gff \
 tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa \
 > tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa.out.gff.TE_prop.out

# Count how many genes exceed various TE proportion thresholds
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
  echo -n "TH $i: "
  awk -v awk_var="$i" '{if ($5> awk_var){print}}' \
  tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa.out.gff.TE_prop.out | wc -l
done > TH.1-.9_of_re_on_genes.count.list

##################################################### Remove gene models with TE overlap > 0.2
awk '{if ($5 > 0.2)print $2}' tsebra.longest_isoform.rm_intronless_incomp.merge_cds.fa.out.gff.TE_prop.out \
 > gene_list.TE_gt_20perc.list

agat_sp_filter_feature_from_kill_list.pl \
 --gff tsebra.longest_isoform.rm_intronless_incomp.gff \
 --kill_list gene_list.TE_gt_20perc.list \
 -o tsebra.longest_isoform.rm_intronless_incomp.filter_20perc_TE.gff

##################################################### Remove genes with length < 100 bp
agat_sp_filter_gene_by_length.pl \
 --gff tsebra.longest_isoform.rm_intronless_incomp.filter_20perc_TE.gff \
 --size 100 --test ">" \
 -o tsebra.longest_isoform.rm_intronless_incomp.filter_20perc_TE.rm_lt_100.gff

##################################################### Final formatting for functional annotation
agat_sp_manage_IDs.pl \
 --gff tsebra.longest_isoform.rm_intronless_incomp.filter_20perc_TE.rm_lt_100.gff \
 --prefix Rn01_ --tair \
 -o Rn01.gff3

agat_sp_manage_attributes.pl \
 -gff Rn01.gff3 \
 -att gene_id,transcript_id \
 -p level2,level3 \
 -o Rn01.f.gff3

# Rename final GFF3
mv Rn01.f.gff3 Rn01.gff3

# Extract final protein sequences
agat_sp_extract_sequences.pl \
 -p \
 -g Rn01.gff3 \
 -f $fa_soft_mask \
 -o Rn01.aa.fa
