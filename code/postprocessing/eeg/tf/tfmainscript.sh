#!/bin/bash
#SBATCH --job-name=TFjob        	# job name
#SBATCH --nodes=1                	# node count
#SBATCH --ntasks=1               	# total number of tasks across all nodes
#SBATCH --mail-type=end          	# send email when job ends
#SBATCH --mail-user=cknowlto@fiu.edu  # email address
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --account=iacc_gbuzzell	# SLURM account name (delete these 3 lines if not running a highmem job)
#SBATCH --partition=highmem1       # partition name (use high memory nodes)
#SBATCH --qos=highmem1             # QOS
#SBATCH --output=%x-%j.out

module load matlab-2023a;
pwd; hostname; date

cd /home/data/NDClab/analyses/read-study2-alpha/code/postprocessing/eeg/tf

module load matlab-2023a;

matlab -nodisplay -nosplash -r "run('MainScript_Calculate_TF_ITPS_ICPS.m'); exit;"
echo "Job completed at: $(date)"