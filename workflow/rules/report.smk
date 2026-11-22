# -----------------------------------------------------------------------------
# One summary table across samples: the numbers the guide asks you to record
# (reads before/after trimming, mapping rate, coverage, raw/filtered variants)
# -----------------------------------------------------------------------------


rule summary:
    input:
        fastp=expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        flagstat=expand("results/alignment/{sample}.flagstat.txt", sample=SAMPLES),
        coverage=expand("results/alignment/{sample}.coverage.tsv", sample=SAMPLES),
        raw=expand("results/variants/{sample}.raw.vcf.gz", sample=SAMPLES),
        filtered=expand("results/variants/{sample}.filtered.norm.vcf.gz", sample=SAMPLES),
        table=expand("results/annotation/{sample}.variant_table.tsv", sample=SAMPLES),
    output:
        "results/summary/summary.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/summary.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    script:
        "../scripts/summary.py"
