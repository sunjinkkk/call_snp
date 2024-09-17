# Define configuration
# 从外部文件读取样本列表
with open('samples.txt', 'r') as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

REF_GENOME = "reference/sorghum_reference.fa"
THREADS = 20  # 设置线程数，适配节点

rule all:
    input:
        "reference/sorghum_reference.fa.pac",
        "reference/sorghum_reference.dict",
        "reference/sorghum_reference.fa.fai",
        expand("03.vcf/{sample}.g.vcf", sample=SAMPLES),
        "03.vcf/final_output.all.vcf",
        "03.vcf/final_output.snp.vcf",
        "03.vcf/final_output.indel.vcf"

# 参考基因组索引 (使用 bwa 和 samtools)
rule index_reference:
    input:
        REF_GENOME
    output:
        "reference/sorghum_reference.fa.pac",
        "reference/sorghum_reference.fa.fai"
    conda:
        "mapping"
    shell:
        """
        bwa-mem2.avx512bw index {input}
        samtools index {input}
        """

rule create_ref_dict:
    input:
        REF_GENOME
    output:
        "reference/sorghum_reference.dict",
    conda:
        "gatk4"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output[0]}"

rule fastp:
    input:
        r1="00.rawdata/{sample}_1.fq.gz",
        r2="00.rawdata/{sample}_2.fq.gz"
    output:
        filtered_r1="01.cleandata/{sample}_1_clean.fq.gz",
        filtered_r2="01.cleandata/{sample}_2_clean.fq.gz",
        report_html="01.cleandata/{sample}_fastp_report.html",
        report_json="01.cleandata/{sample}_fastp_report.json"
    conda:
        "fastp"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.filtered_r1} -O {output.filtered_r2} \
              --html {output.report_html} --json {output.report_json} -q 20 -u 30
        """

# Mapping reads 到参考基因组 (使用 bwa 和 samtools)
rule map_reads:
    input:
        ref=REF_GENOME,
        r1="01.cleandata/{sample}_1_clean.fq.gz",
        r2="01.cleandata/{sample}_2_clean.fq.gz"
    output:
        "02.bam/{sample}_sorted.bam"
    threads: THREADS
    conda:
        "mapping"
    shell:
        """
        bwa-mem2.avx512bw mem -t {threads} -R "@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA" {input.ref} {input.r1} {input.r2} | samtools view -Sb | samtools sort -o {output} --output-fmt BAM
        """

rule mark_duplicates:
    input:
        "02.bam/{sample}_sorted.bam"
    output:
        bam="02.bam/{sample}_dedup.bam",
        metrics="02.bam/{sample}_metrics.txt"
    conda:
        "gatk4"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}"

rule index_bam:
    input:
        dedup_bam="02.bam/{sample}_dedup.bam"
    output:
        bai="02.bam/{sample}_dedup.bam.bai"
    conda:
        "mapping"
    shell:
        "samtools index {input.dedup_bam}"

rule haplotype_caller:
    input:
        bam="02.bam/{sample}_dedup.bam",
        bai="02.bam/{sample}_dedup.bam.bai",
        ref=REF_GENOME
    output:
        "03.vcf/{sample}.g.vcf"
    conda:
        "gatk4"
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -ERC GVCF"

rule combine_gvcfs:
    input:
        gvcfs=expand("03.vcf/{sample}.g.vcf", sample=SAMPLES),
        ref=REF_GENOME
    output:
        "03.vcf/combined.g.vcf"
    params:
        gvcf_args=lambda wildcards, input: " ".join(f"-V {gvcf}" for gvcf in input.gvcfs)
    conda:
        "gatk4"
    shell:
        """
        gatk CombineGVCFs -R {input.ref} {params.gvcf_args} -O {output}
        """

rule genotype_gvcfs:
    input:
        gvcf="03.vcf/combined.g.vcf",
        ref=REF_GENOME
    output:
        "03.vcf/final_output.all.vcf"
    conda:
        "gatk4"
    shell:
       "gatk GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"

rule split_vcf:
    input:
        vcf="03.vcf/final_output.all.vcf"
    output:
        snp_vcf="03.vcf/final_output.snp.vcf",
        indel_vcf="03.vcf/final_output.indel.vcf"
    conda:
        "qtlseq"
    shell:
        """
        # 提取 SNP
        bcftools view -v snps {input.vcf} -o {output.snp_vcf} -O v
        
        # 提取 InDel
        bcftools view -v indels {input.vcf} -o {output.indel_vcf} -O v
        """

