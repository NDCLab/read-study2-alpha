#!/bin/bash
#SBATCH --job-name=TFjob
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end          	# send email when job ends
#SBATCH --mail-user=cknowlto@fiu.edu  # email address
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --account=iacc_gbuzzell	# SLURM account name (delete these 3 lines if not running a highmem job)
#SBATCH --partition=highmem1       # partition name (use high memory nodes)
#SBATCH --qos=highmem1             # QOS
#SBATCH --output=%x-%j.out

module load matlab-2023a;

mkdir -p /home/data/NDClab/analyses/read-study2-alpha/derivatives/tf/logs
 
echo "Starting parallel TF analysis with $SLURM_CPUS_PER_TASK workers"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

matlab -nodisplay -nosplash -r "cd('/home/data/NDClab/analyses/read-study2-alpha/code/postprocessing/eeg/tf'); MainScript_Calculate_TF_ITPS_ICPS_para; exit;"
 
echo "End time: $(date)"
echo "Time-frequency analysis completed"