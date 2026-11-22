# -----------------------------------------------------------------------------
# Variant calling, summarising, filtering and normalisation (BCFtools)
# Guide sections 10 and 11
# -----------------------------------------------------------------------------


rule call_variants:
    """bcftools mpileup | bcftools call -mv (haploid by default, see config)."""
    input:
        bam="results/alignment/{sample}.sorted.bam",
        bai="results/alignment/{sample}.sorted.bam.bai",
        ref=REF,
        fai=REF_FAI,
    output:
        vcf="results/variants/{sample}.raw.vcf.gz",
        csi="results/variants/{sample}.raw.vcf.gz.csi",
    params:
        mq=config["mpileup"]["min_mapq"],
        bq=config["mpileup"]["min_baseq"],
        maxdp=config["mpileup"]["max_depth"],
        mpileup_extra=config["mpileup"].get("extra", ""),
        ploidy=config["call"]["ploidy"],
        call_extra=config["call"].get("extra", ""),
    threads: THREADS
    log:
        "logs/call_variants/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        """
        (bcftools mpileup -Ou -f {input.ref} \
            -q {params.mq} -Q {params.bq} -d {params.maxdp} \
            -a FORMAT/AD,FORMAT/DP {params.mpileup_extra} {input.bam} \
         | bcftools call -mv --ploidy {params.ploidy} {params.call_extra} \
            --threads {threads} -Oz -o {output.vcf}) 2> {log}
        bcftools index -f {output.vcf} 2>> {log}
        """


rule split_variant_types:
    """Separate SNP and indel VCFs from the raw calls (guide section 10)."""
    input:
        "results/variants/{sample}.raw.vcf.gz",
    output:
        "results/variants/{sample}.{kind}.vcf.gz",
    wildcard_constraints:
        kind="snps|indels",
    log:
        "logs/split_variant_types/{sample}.{kind}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        "bcftools view -v {wildcards.kind} {input} -Oz -o {output} 2> {log}"


rule filter_variants:
    """Project-specific demonstration filter (default QUAL>=30 && INFO/DP>=10)."""
    input:
        "results/variants/{sample}.raw.vcf.gz",
    output:
        vcf="results/variants/{sample}.filtered.vcf.gz",
        csi="results/variants/{sample}.filtered.vcf.gz.csi",
    params:
        expr=config["filter"]["include"],
    log:
        "logs/filter_variants/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        """
        bcftools filter -i {params.expr:q} {input} -Oz -o {output.vcf} 2> {log}
        bcftools index -f {output.vcf} 2>> {log}
        """


rule normalize_variants:
    """Left-align and normalise indels against the reference (guide section 11)."""
    input:
        vcf="results/variants/{sample}.filtered.vcf.gz",
        ref=REF,
        fai=REF_FAI,
    output:
        vcf="results/variants/{sample}.filtered.norm.vcf.gz",
        csi="results/variants/{sample}.filtered.norm.vcf.gz.csi",
    log:
        "logs/normalize_variants/{sample}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        """
        bcftools norm -f {input.ref} {input.vcf} -Oz -o {output.vcf} 2> {log}
        bcftools index -f {output.vcf} 2>> {log}
        """


rule variant_stats:
    input:
        "results/variants/{sample}.{kind}.vcf.gz",
    output:
        "results/variants/{sample}.{kind}_stats.txt",
    wildcard_constraints:
        kind="raw|filtered",
    log:
        "logs/variant_stats/{sample}.{kind}.log",
    container:
        "docker://saharbagheri/bacterial-variants:1.0"
    shell:
        "bcftools stats {input} > {output} 2> {log}"
