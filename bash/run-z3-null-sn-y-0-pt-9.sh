#!/bin/sh

#!/bin/bash
#SBATCH --job-name 103
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000m
#SBATCH --time=12:00:00
#SBATCH --account=jzelner1
#SBATCH --partition=standard
#SBATCH --mail-type=NONE
#SBATCH --export=ALL
#SBATCH --output=%x-%j.log
#SBATCH --array=1-100

#i=$(($SLURM_ARRAY_TASK_ID + 50))
i=$SLURM_ARRAY_TASK_ID
cd $SLURM_SUBMIT_DIR
Rscript --verbose run-server-design.R "5k_design_Z_3_sparse_null_sp_Y_0_pt_9.RDS" $i "covariate-meas-error-a-looser-priors-sn-sp-theta.RDS" "5k_Z_3_cov_a_sparse_looser_priors_theta_sp_Y_0_pt_9_null" 
Rscript --verbose run-server-design.R "10k_design_Z_3_sparse_null_sp_Y_0_pt_9.RDS" $i "covariate-meas-error-a-looser-priors-sn-sp-theta.RDS" "10k_Z_3_cov_a_sparset_looser_sn_sp_theta_sp_Y_0_pt_9_null" 
Rscript --verbose run-server-design.R "5k_design_Z_3_sparse_noisy_null_sp_Y_0_pt_9.RDS" $i "covariate-meas-error-a-looser-priors-sn-sp-theta.RDS" "5k_Z_3_cov_a_sparse_noisy_looser_priors_theta_sp_Y_0_pt_9_null" 
Rscript --verbose run-server-design.R "10k_design_Z_3_sparse_noisy_null_sp_Y_0_pt_9.RDS" $i "covariate-meas-error-a-looser-priors-sn-sp-theta.RDS" "10k_Z_3_cov_a_sparse_noisy_looser_priors_theta_sp_Y_0_pt_9_null" 
Rscript --verbose run-server-design.R "15k_design_Z_3_sparse_null_sp_Y_0_pt_9.RDS" $i "covariate-meas-error-a-looser-priors-sn-sp-theta.RDS" "15k_Z_3_cov_a_sparse_looser_priors_theta_sp_Y_0_pt_9_null" 
Rscript --verbose run-server-design.R "15k_design_Z_3_sparse_noisy_null_sp_Y_0_pt_9.RDS" $i "covariate-meas-error-a-looser-priors-sn-sp-theta.RDS" "15k_Z_3_cov_a_sparse_noisy_looser_priors_theta_sp_Y_0_pt_9_null" 
# Run rscript on each node
