#!/usr/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=SV.Demultiplexing_pt1_2_.script
#SBATCH --output=SV.Demultiplexing_pt1_2.output
#SBATCH --error=SV.Demultiplexing_pt1_2.error
#SBATCH --time=0-012:00:00
#SBATCH --nodes=1

conda deactivate
conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3

cd /projects/bgmp/smv/Bi621
/usr/bin/time -v python Demultiplexing_pt1.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r 8 -H '1294_S1_L008_R2_001' > output.1294_S1_L008_R2_001.tsv
