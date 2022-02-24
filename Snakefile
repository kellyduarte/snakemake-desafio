samples = ["A", "B", "C"]

rule all:
    input:
        "results/calls/all.vcf",
        "results/plots/quals.svg"

"""        
rule all:
    input:
        expand("results/{sample}.bam", sample = samples)
"""

rule map_reads:
    input:
       refgene = "data/genome.fa",
       amostra = "data/samples/{sample}.fastq"
    output:
        "results/mapped/{sample}.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa mem {input.refgene} {input.amostra} | samtools view -Sb - > {output}"

rule sort_alignments:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"

rule call:
    input:
        fa="data/genome.fa",
        bam=expand("results/mapped/{sample}.sorted.bam", sample=samples)
    output:
        "results/calls/all.vcf"
    conda:
        "envs/calling.yaml"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output}"
    
rule plot_quals:
    input:
        "results/calls/all.vcf"
    output:
        "results/plots/quals.svg"
    conda:
        "envs/stats.yaml"
    notebook:
        "notebooks/plot-quals.py.ipynb"

    #não está geraando o output "results/plots/quals.svg"