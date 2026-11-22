"""
Build a readable interpretation table from `bcftools query` output.

Input 1 (query TSV, no header):
  CHROM POS REF ALT TYPE QUAL INFO/DP INFO/MQ FORMAT/AD INFO/BCSQ
Input 2: gene table written by ncbi_gff_to_csq.py

For each variant this reports read support (depth, ALT reads, ALT fraction),
the gene it falls in (with product description), the bcftools csq
consequence + amino-acid change, or - for intergenic variants - the two
flanking genes and their distances (guide section 13).

Usage: python variant_table.py query.tsv genes.tsv > variant_table.tsv
"""

import bisect
import sys
from collections import defaultdict


def load_genes(path):
    genes = defaultdict(list)
    with open(path) as fh:
        next(fh)
        for line in fh:
            chrom, start, end, strand, tag, name, biotype, product = line.rstrip("\n").split("\t")
            genes[chrom].append((int(start), int(end), strand, tag, name, biotype, product))
    for chrom in genes:
        genes[chrom].sort()
    starts = {c: [g[0] for g in gl] for c, gl in genes.items()}
    return genes, starts


def locate(genes, starts, chrom, pos, ref_len):
    """Return (region, locus_tags, names, biotypes, products, flank description)."""
    gl = genes.get(chrom, [])
    end = pos + ref_len - 1
    hits = [g for g in gl if g[0] <= end and g[1] >= pos]
    if hits:
        return ("genic", ",".join(h[3] for h in hits), ",".join(h[4] for h in hits),
                ",".join(h[5] for h in hits), " | ".join(h[6] for h in hits), "")
    i = bisect.bisect_right(starts.get(chrom, []), pos)
    left = max((g for g in gl[:i] if g[1] < pos), key=lambda g: g[1], default=None)
    right = gl[i] if i < len(gl) else None
    parts = []
    if left:
        parts.append(f"left: {left[4]} ({left[2]} strand, {pos - left[1]} bp away)")
    if right:
        parts.append(f"right: {right[4]} ({right[2]} strand, {right[0] - pos} bp away)")
    return ("intergenic", "", "", "", "", " / ".join(parts))


def parse_bcsq(bcsq):
    """First csq record -> (consequence, gene, aa_change, dna_change)."""
    if bcsq in ("", "."):
        return "", "", "", ""
    first = bcsq.split(",")[0]
    if first.startswith("@"):  # consequence reported at another position (compound)
        return f"see POS {first[1:]}", "", "", ""
    f = first.split("|")
    f += [""] * (7 - len(f))
    return f[0], f[1], f[5], f[6]


def main(query, genes_path):
    genes, starts = load_genes(genes_path)
    cols = ["chrom", "pos", "ref", "alt", "type", "qual", "depth", "ref_reads", "alt_reads",
            "alt_fraction", "mapq", "region", "locus_tag", "gene", "biotype", "product",
            "consequence", "protein_change", "dna_change", "nearest_genes", "bcsq"]
    print("\t".join(cols))
    with open(query) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt, vtype, qual, dp, mq, ad, bcsq = f[:10]
            ads = [int(x) for x in ad.split(",") if x.isdigit()]
            ref_reads = ads[0] if ads else ""
            alt_reads = sum(ads[1:]) if len(ads) > 1 else ""
            frac = f"{alt_reads / sum(ads):.2f}" if ads and sum(ads) else ""
            region, tags, names, biotypes, products, flank = locate(genes, starts, chrom, int(pos), len(ref))
            csq, csq_gene, aa, dna = parse_bcsq(bcsq)
            if len(aa) > 40:  # long frameshift strings; the full text stays in the bcsq column
                aa = aa[:40] + "..."
            if region == "intergenic" and not csq:
                csq = "intergenic"
            try:
                qual = f"{float(qual):.1f}"
            except ValueError:
                pass
            print("\t".join(map(str, [
                chrom, pos, ref, alt, vtype, qual, dp, ref_reads, alt_reads, frac, mq,
                region, tags, names or csq_gene, biotypes, products, csq, aa, dna, flank, bcsq,
            ])))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
