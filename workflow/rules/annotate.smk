# -----------------------------------------------------------------------------
# Consequence annotation with `bcftools csq` + an interpretation table
# Guide sections 12 and 13
#
# bcftools csq expects an Ensembl-style GFF3 (gene -> transcript -> exon/CDS,
# with ID=gene:..., ID=transcript:..., biotype=...). NCBI's RefSeq GFF3 for
# bacteria uses a different layout (gene -> CDS, gene_biotype=...), which is
# the "format or feature error" the guide warns about. The rule below converts
# it, and first checks that the FASTA and GFF sequence IDs agree.
# -----------------------------------------------------------------------------


rule prepare_csq_gff:
    input:
        gff=GFF,
        fai=REF_FAI,
    output:
        gff=CSQ_GFF,
        genes=f"reference/{REF_NAME}.genes.tsv",
    log:
        "logs/reference/prepare_csq_gff.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    script:
        "../scripts/ncbi_gff_to_csq.py"


rule annotate_variants:
    input:
        vcf="results/variants/{sample}.filtered.norm.vcf.gz",
        ref=REF,
        fai=REF_FAI,
        gff=CSQ_GFF,
    output:
        vcf="results/annotation/{sample}.annotated.vcf.gz",
        csi="results/annotation/{sample}.annotated.vcf.gz.csi",
    log:
        "logs/annotate_variants/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        """
        bcftools csq -f {input.ref} -g {input.gff} --phase a \
            {input.vcf} -Oz -o {output.vcf} 2> {log}
        bcftools index -f {output.vcf} 2>> {log}
        """


rule variant_table:
    """
    Readable table for interpretation: position, REF/ALT, type, QUAL, depth,
    allele frequency, gene / product, csq consequence and amino-acid change,
    or the flanking genes for intergenic variants.
    """
    input:
        vcf="results/annotation/{sample}.annotated.vcf.gz",
        genes=f"reference/{REF_NAME}.genes.tsv",
    output:
        query=temp("results/annotation/{sample}.query.tsv"),
        table="results/annotation/{sample}.variant_table.tsv",
    params:
        script=workflow.source_path("../scripts/variant_table.py"),
    log:
        "logs/variant_table/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        """
        bcftools query \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%TYPE\\t%QUAL\\t%INFO/DP\\t%INFO/MQ\\t[%AD]\\t%INFO/BCSQ\\n' \
            {input.vcf} > {output.query} 2> {log}
        python {params.script} \
            {output.query} {input.genes} > {output.table} 2>> {log}
        """

