# -----------------------------------------------------------------------------
# Alignment (BWA-MEM) -> sorted, indexed BAM -> alignment QC
# Guide sections 6, 7 and 8
#
# Instead of writing SAM -> BAM -> sorted BAM to disk (as in the guide),
# BWA output is streamed straight into `samtools sort`. The result is the same
# sorted BAM, without the large intermediate SAM/BAM files.
# A read group (@RG) is added so the sample name ends up in the VCF.
# -----------------------------------------------------------------------------


rule bwa_mem:
    input:
        r1="data/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="data/trimmed/{sample}_R2.trimmed.fastq.gz",
        ref=REF,
        idx=BWA_INDEX,
    output:
        bam="results/alignment/{sample}.sorted.bam",
    params:
        rg=lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ILLUMINA",
        tmp=lambda wc: f"results/alignment/{wc.sample}.sorttmp",
    threads: THREADS
    log:
        "logs/bwa_mem/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        """
        (bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -T {params.tmp} -o {output.bam} -) 2> {log}
        """


rule samtools_index:
    input:
        "results/alignment/{sample}.sorted.bam",
    output:
        "results/alignment/{sample}.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        "samtools index {input} > {log} 2>&1"


rule alignment_stats:
    """flagstat, idxstats and stats (guide section 8)."""
    input:
        bam="results/alignment/{sample}.sorted.bam",
        bai="results/alignment/{sample}.sorted.bam.bai",
        ref=REF,
    output:
        flagstat="results/alignment/{sample}.flagstat.txt",
        idxstats="results/alignment/{sample}.idxstats.txt",
        stats="results/alignment/{sample}.stats.txt",
    log:
        "logs/alignment_stats/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        samtools stats -r {input.ref} {input.bam} > {output.stats} 2>> {log}
        """


rule coverage:
    """
    Per-position depth (samtools depth -a, i.e. including zero-coverage
    positions, so the mean is over the whole genome), a coverage summary,
    and a BED of low-coverage regions.
    """
    input:
        bam="results/alignment/{sample}.sorted.bam",
        bai="results/alignment/{sample}.sorted.bam.bai",
    output:
        depth="results/alignment/{sample}.depth.txt.gz",
        summary="results/alignment/{sample}.coverage.tsv",
        lowcov="results/alignment/{sample}.low_coverage.bed",
    params:
        awk=workflow.source_path("../scripts/coverage_summary.awk"),
        low=config["coverage"]["low_depth"],
    log:
        "logs/coverage/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-alignment:1.0"
    shell:
        """
        samtools depth -a {input.bam} 2> {log} | gzip -c > {output.depth}
        gzip -cd {output.depth} \
            | awk -v sample={wildcards.sample} -v low={params.low} -v bed={output.lowcov} \
                  -f {params.awk} > {output.summary} 2>> {log}
        """
