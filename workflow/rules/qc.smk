# -----------------------------------------------------------------------------
# Read QC (FastQC before/after trimming) and trimming (fastp)
# Guide sections 4 and 5
# -----------------------------------------------------------------------------


rule fastqc_raw:
    input:
        "data/raw/{sample}_{read}.fastq.gz",
    output:
        html="results/qc/fastqc_raw/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc_raw/{sample}_{read}_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_{read}.log",
    container:
        "docker://saharbagheri/bacterial-qc:1.0"
    shell:
        "fastqc --quiet --outdir results/qc/fastqc_raw {input} > {log} 2>&1"


rule fastp:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    output:
        r1="data/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="data/trimmed/{sample}_R2.trimmed.fastq.gz",
        html="results/qc/fastp/{sample}.fastp.html",
        json="results/qc/fastp/{sample}.fastp.json",
    params:
        extra=config["fastp"].get("extra", ""),
    threads: THREADS
    log:
        "logs/fastp/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-qc:1.0"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.html} -j {output.json} \
              -w {threads} {params.extra} > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        "data/trimmed/{sample}_{read}.trimmed.fastq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.html",
        zip="results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.zip",
    log:
        "logs/fastqc_trimmed/{sample}_{read}.log",
    container:
        "docker://saharbagheri/bacterial-qc:1.0"
    shell:
        "fastqc --quiet --outdir results/qc/fastqc_trimmed {input} > {log} 2>&1"


rule multiqc_raw:
    """MultiQC of the raw-read FastQC reports only."""
    input:
        expand("results/qc/fastqc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
    output:
        "results/qc/multiqc_raw/multiqc_report.html",
    log:
        "logs/multiqc/raw.log",
    container:
        "docker://saharbagheri/bacterial-qc:1.0"
    shell:
        "multiqc --force -o results/qc/multiqc_raw {input} > {log} 2>&1"


rule multiqc_final:
    """One MultiQC report with raw + trimmed FastQC, fastp, SAMtools and BCFtools stats."""
    input:
        expand("results/qc/fastqc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        expand("results/alignment/{sample}.{stat}.txt", sample=SAMPLES, stat=["flagstat", "idxstats", "stats"]),
        expand("results/variants/{sample}.{kind}_stats.txt", sample=SAMPLES, kind=["raw", "filtered"]),
    output:
        "results/qc/multiqc_final/multiqc_report.html",
    log:
        "logs/multiqc/final.log",
    container:
        "docker://saharbagheri/bacterial-qc:1.0"
    shell:
        "multiqc --force -o results/qc/multiqc_final {input} > {log} 2>&1"
