#!/bin/bash
#SBATCH --job-name=split
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --error=slurm-%J.err
#SBATCH --mem 128G
#SBATCH --array=0-7
#SBATCH --time=72:00:00

# activate conda env:
source /opt/apps/anaconda/2020.07/etc/profile.d/conda.sh
conda activate splitseq

# load samtools and star on HPC
module load samtools
module load star

# move to the project directory
cd /dfs3b/swaruplab/smorabit/data/FIRE_mouse_2021/

# set data and reference directories
fastqs="./expdata_combined/"
ref_dir="./genomes/mm10_mRNA/"

let index="$SLURM_ARRAY_TASK_ID"

# find current sublibrary
sublibraries=($(ls $fastqs | cut -d '-' -f 1 | uniq))
sample=${sublibraries[$index]}

# make output dir for this sample
mkdir ./analysis/$sample

splitpipe/split-pipe \
  --mode all \
  --nthreads 32 \
  --genome_dir $ref_dir \
  --fq1 $fastqs$sample-READ1.fastq.gz \
  --fq2 $fastqs$sample-READ2.fastq.gz \
  --output_dir analysis/$sample \
  --sample D-1 A1-A1 \
  --sample D-7 A2-A2 \
  --sample A-1 A3-A3 \
  --sample C-2 A4-A4 \
  --sample E-2 A5-A5 \
  --sample B-2 A6-A6 \
  --sample F-8 A7-A7 \
  --sample E-7 A8-A8 \
  --sample F-7 A9-A9 \
  --sample B-6 A10-A10 \
  --sample D-5 A11-A11 \
  --sample A-2 A12-A12 \
  --sample D-4 B1-B1 \
  --sample F-4 B2-B2 \
  --sample C-5 B3-B3 \
  --sample A-8 B4-B4 \
  --sample B-8 B5-B5 \
  --sample E-6 B6-B6 \
  --sample C-1 B7-B7 \
  --sample C-7 B8-B8 \
  --sample A-9 B9-B9 \
  --sample E-1 B10-B10 \
  --sample D-3 B11-B11 \
  --sample F-5 B12-B12 \
  --sample F-2 C1-C1 \
  --sample E-4 C2-C2 \
  --sample B-1 C3-C3 \
  --sample B-4 C4-C4 \
  --sample A-6 C5-C5 \
  --sample F-6 C6-C6 \
  --sample C-6 C7-C7 \
  --sample B-7 C8-C8 \
  --sample B-3 C9-C9 \
  --sample C-3 C10-C10 \
  --sample F-3 C11-C11 \
  --sample A-7 C12-C12 \
  --sample D-8 D1-D1 \
  --sample B-5 D2-D2 \
  --sample E-8 D3-D3 \
  --sample A-4 D4-D4 \
  --sample D-6 D5-D5 \
  --sample C-8 D6-D6 \
  --sample F-1 D7-D7 \
  --sample C-4 D8-D8 \
  --sample E-3 D9-D9 \
  --sample E-5 D10-D10 \
  --sample A-5 D11-D11 \
  --sample D-2 D12-D12


################################################################################
# format reference transcriptome
################################################################################

# make mRNA reference using split-pipe mkref
# splitpipe/split-pipe \
# --mode mkref \
# --genome mm10 \
# --fasta /data/homezvol1/smorabit/swaruplab/smorabit/data/FIRE_mouse_2021/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
# --genes /data/homezvol1/smorabit/swaruplab/smorabit/data/FIRE_mouse_2021/genomes/Mus_musculus.GRCm38.93.gtf.gz \
# --output_dir ./genomes/mm10_mRNA
