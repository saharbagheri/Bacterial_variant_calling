# =============================================================================
# Bacterial variant calling: FASTQ -> QC -> trimming -> BWA-MEM -> BAM QC
#                            -> BCFtools calling -> filtering -> annotation
#
# Run from the project directory:
#   snakemake -np            # dry run: show what will be done
#   snakemake                # run (uses workflow/profiles/default)
# =============================================================================

from snakemake.utils import min_version

min_version("8.0")


configfile: "config/config.yaml"


include: "workflow/rules/common.smk"
include: "workflow/rules/data.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/variants.smk"
include: "workflow/rules/annotate.smk"
include: "workflow/rules/report.smk"


rule all:
    input:
        # QC
        "results/qc/multiqc_raw/multiqc_report.html",
        "results/qc/multiqc_final/multiqc_report.html",

        # Alignment + alignment QC
        expand("results/alignment/{sample}.sorted.bam.bai", sample=SAMPLES),
        expand("results/alignment/{sample}.coverage.tsv", sample=SAMPLES),

        # Variants
        expand("results/variants/{sample}.{kind}.vcf.gz", sample=SAMPLES, kind=["snps", "indels"]),
        expand("results/variants/{sample}.filtered.norm.vcf.gz.csi", sample=SAMPLES),

        # Annotation + interpretation table
        expand("results/annotation/{sample}.annotated.vcf.gz.csi", sample=SAMPLES),
        expand("results/annotation/{sample}.variant_table.tsv", sample=SAMPLES),

        # One-page summary across samples
        "results/summary/summary.tsv",
