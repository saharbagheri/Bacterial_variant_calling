# -----------------------------------------------------------------------------
# Get the data: reads (SRA or local) and the reference genome + annotation
# -----------------------------------------------------------------------------


rule get_reads:
    """prefetch + fasterq-dump --split-files, then gzip (guide section 3)."""
    output:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    params:
        acc=lambda wc: SAMPLE_INFO[wc.sample]["sra"],
        tmp=lambda wc: f"data/raw/tmp_{wc.sample}",
    wildcard_constraints:
        sample=regex_for(SRA_SAMPLES),
    threads: THREADS
    log:
        "logs/get_reads/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-sra:1.0"
    retries: 2
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp}
        prefetch {params.acc} --output-directory {params.tmp} > {log} 2>&1
        fasterq-dump {params.tmp}/{params.acc} --split-files --threads {threads} \
            --outdir {params.tmp} --temp {params.tmp} >> {log} 2>&1
        pigz -p {threads} -c {params.tmp}/{params.acc}_1.fastq > {output.r1}
        pigz -p {threads} -c {params.tmp}/{params.acc}_2.fastq > {output.r2}
        rm -rf {params.tmp}
        """


rule link_local_fastq:
    """Use FASTQ files already on disk (fq1/fq2 columns of the sample sheet)."""
    input:
        r1=lambda wc: workflow.source_path(SAMPLE_INFO[wc.sample]["fq1"]),
        r2=lambda wc: workflow.source_path(SAMPLE_INFO[wc.sample]["fq2"]),
    output:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    wildcard_constraints:
        sample=regex_for(LOCAL_SAMPLES),
    localrule: True
    run:
        import os

        for src, dst in [(input.r1, output.r1), (input.r2, output.r2)]:
            if not src.endswith(".gz"):
                raise ValueError(f"{src}: local FASTQ files must be gzipped (.fastq.gz)")
            os.symlink(os.path.abspath(src), dst)


rule get_reference:
    """Download (or copy) the reference FASTA and GFF3 (guide sections 3 and 12)."""
    input:
        fasta=lambda wc: reference_input("fasta"),
        gff=lambda wc: reference_input("gff"),
    output:
        fasta=REF,
        gff=GFF,
    params:
        fasta_cmd=lambda wc, input: fetch_cmd(
            input.fasta or config["reference"]["fasta"]
        ),
        gff_cmd=lambda wc, input: fetch_cmd(
            input.gff or config["reference"]["gff"]
        ),
    log:
        "logs/reference/get_reference.log",
    retries: 2
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        """
        ({params.fasta_cmd} > {output.fasta}) 2> {log}
        ({params.gff_cmd} > {output.gff}) 2>> {log}
        echo "FASTA sequence IDs:" >> {log}; grep '^>' {output.fasta} >> {log}
        """


rule bwa_index:
    input:
        REF,
    output:
        BWA_INDEX,
    log:
        "logs/reference/bwa_index.log",
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        "bwa index {input} > {log} 2>&1"


rule faidx:
    input:
        REF,
    output:
        REF_FAI,
    log:
        "logs/reference/faidx.log",
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        "samtools faidx {input} > {log} 2>&1"
