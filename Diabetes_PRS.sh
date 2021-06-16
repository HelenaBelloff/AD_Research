#!/usr/bin/env sh
# Polygenic risk score calculation for Diabetes Genetic-Wide Association Study (GWAS)
module load prsice

Rscript /hpc/packages/minerva-centos7/prsice/2.2.6/src/PRSice/PRSice.R \
    --prsice /hpc/packages/minerva-centos7/prsice/2.2.6/bin/PRSice \
    --A1 RISK_ALLELE \
    --A2 OTHER_ALLELE \
    --bar-levels 0.001000,0.050000,0.100000,0.200000,0.300000,0.400000,0.500000,1.000000 \
    --base /hpc/users/belloh02/GWAS_Data/Diabetes/DIAGRAMv3.2012DEC172.txt \
    --binary-target T \
    --bp POSITION \
    --chr CHROMOSOME \
    --interval 0.000050 \
    --lower 0.000100 \
    --model add \
    --clump-kb 1M \
    --missing MEAN_IMPUTE \
    --out /hpc/users/belloh02/PRS_OUTPUTS/Diabetes_PRS_Output/Diabetes_PRS \
    --perm 50 \
    --pvalue P_VALUE \
    --seed 1234567890 \
    --snp SNP \
    --stat OR \
    --target /sc/hydra/projects/zhangb03a/shared/msbb-wgs/PLINK/AMPADWGS_WGS_all_rsIDs_dbSNP153_020520 \
    --thread 1 \
    --upper 0.500000 \
    --score std