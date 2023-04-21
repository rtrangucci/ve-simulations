#!/bin/sh

#!/bin/bash
#SBATCH --job-name 103
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=6:00:00
#SBATCH --account=jzelner1
#SBATCH --partition=standard
#SBATCH --mail-type=NONE
#SBATCH --export=ALL
#SBATCH --output=%x-%j.log
#SBATCH --array=1-100

#i=$(($SLURM_ARRAY_TASK_ID + 50))
i=$SLURM_ARRAY_TASK_ID
cd $SLURM_SUBMIT_DIR
Rscript --verbose run-server-design.R "1k_design_Z_2_transmission_sparse.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "1k_Z_2_transmission_cov_a_sparse_looser_theta" 
Rscript --verbose run-server-design.R "1k_design_Z_2_transmission_sparse_noisy.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "1k_Z_2_transmission_cov_a_sparse_noisy_looser_theta" 
Rscript --verbose run-server-design.R "5k_design_Z_2_transmission_sparse.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "5k_Z_2_transmission_cov_a_sparse_looser_theta" 
Rscript --verbose run-server-design.R "5k_design_Z_2_transmission_sparse_noisy.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "5k_Z_2_transmission_cov_a_sparse_noisy_looser_theta" 
Rscript --verbose run-server-design.R "10k_design_Z_2_transmission_sparse.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "10k_Z_2_transmission_cov_a_sparse_looser_theta" 
Rscript --verbose run-server-design.R "10k_design_Z_2_transmission_sparse_noisy.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "10k_Z_2_transmission_cov_a_sparse_noisy_looser_theta" 
Rscript --verbose run-server-design.R "20k_design_Z_2_transmission_sparse.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "20k_Z_2_transmission_cov_a_sparse_looser_theta" 
Rscript --verbose run-server-design.R "20k_design_Z_2_transmission_sparse_noisy.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "20k_Z_2_transmission_cov_a_sparse_noisy_looser_theta" 

Rscript --verbose run-server-design.R "1k_design_Z_2_transmission_sparse_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "1k_Z_2_transmission_cov_a_sparse_looser_theta_null" 
Rscript --verbose run-server-design.R "1k_design_Z_2_transmission_sparse_noisy_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "1k_Z_2_transmission_cov_a_sparse_noisy_looser_theta_null" 
Rscript --verbose run-server-design.R "5k_design_Z_2_transmission_sparse_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "5k_Z_2_transmission_cov_a_sparse_looser_theta_null" 
Rscript --verbose run-server-design.R "5k_design_Z_2_transmission_sparse_noisy_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "5k_Z_2_transmission_cov_a_sparse_noisy_looser_theta_null" 
Rscript --verbose run-server-design.R "10k_design_Z_2_transmission_sparse_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "10k_Z_2_transmission_cov_a_sparse_looser_theta_null" 
Rscript --verbose run-server-design.R "10k_design_Z_2_transmission_sparse_noisy_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "10k_Z_2_transmission_cov_a_sparse_noisy_looser_theta_null" 
Rscript --verbose run-server-design.R "20k_design_Z_2_transmission_sparse_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "20k_Z_2_transmission_cov_a_sparse_looser_theta_null" 
Rscript --verbose run-server-design.R "20k_design_Z_2_transmission_sparse_noisy_null.RDS" $i "covariate-meas-error-a-looser-priors-equal-sn-sp-theta.RDS" "20k_Z_2_transmission_cov_a_sparse_noisy_looser_theta_null" 
# Run rscript on each node
