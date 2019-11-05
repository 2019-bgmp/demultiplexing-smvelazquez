#!/usr/bin/bash

#ssh hpc
#cd /projects/bgmp/shared/Bi623/assign2

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=Demultiplexing_pt2
#SBATCH --output=Demultiplexing_pt2.output
#SBATCH --error=Demultiplexing_pt2.error
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1

conda deactivate
conda activate bgmp_py3

/usr/bin/time -v python Demultiplexing_pt2.py --FILEPATH_R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz --FILEPATH_R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz --FILEPATH_I1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz  --FILEPATH_I2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz --FILEPATH_INDEXES /projects/bgmp/shared/2017_sequencing/indexes.txt
