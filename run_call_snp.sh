#!/bin/bash
#SBATCH --job-name=wsm_callsnp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56  # 使用 56 线程
#SBATCH --nodes=1
#SBATCH -p Cnode
#SBATCH --time=48:00:00  # 根据实际需要调整
#SBATCH --mem=0  # 让 SLURM 自动分配所有内存
#SBATCH --output=slurm-%j.out

# 激活 conda
source ~/software/miniforge/etc/profile.d/conda.sh
conda activate snakemake
# 运行 snakemake
snakemake --snakefile call_snp.smk --cores 56 --use-conda
