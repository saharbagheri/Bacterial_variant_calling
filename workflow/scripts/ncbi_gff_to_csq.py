"""
Convert an NCBI RefSeq (prokaryotic) GFF3 into the Ensembl-style GFF3 layout
that `bcftools csq` understands, and write a simple gene table.

NCBI bacterial layout             ->  csq layout
  gene  ID=gene-X;gene_biotype=..     gene  ID=gene:X;biotype=...;Name=...
  CDS   Parent=gene-X                 mRNA  ID=transcript:X;Parent=gene:X;biotype=...
                                      exon  Parent=transcript:X
                                      CDS   Parent=transcript:X
  tRNA/rRNA Parent=gene-X             ncRNA ID=transcript:X;Parent=gene:X;biotype=...
  exon  Parent=rna-X                  exon  Parent=transcript:X

Pseudogenes and CDS flagged pseudo=true are written as non-coding, so csq
does not predict protein changes for them.

Before converting, the FASTA (.fai) and GFF sequence IDs are compared; the
script stops with a clear message if they do not match (guide section 12).

Usage: run by Snakemake (rule prepare_csq_gff), or standalone:
  python ncbi_gff_to_csq.py in.gff ref.fasta.fai out.csq.gff3 out.genes.tsv
"""

import sys
from collections import defaultdict
from urllib.parse import unquote

# NCBI gene_biotype -> biotype names bcftools csq recognises
# (csq silently ignores transcripts with biotypes it does not know, such as
#  "tRNA" or "ncRNA", so those are mapped to misc_RNA; the real NCBI biotype
#  is kept in the gene table)
BIOTYPE_MAP = {
    "protein_coding": "protein_coding",
    "tRNA": "misc_RNA",
    "rRNA": "rRNA",
    "ncRNA": "misc_RNA",
    "tmRNA": "misc_RNA",
    "SRP_RNA": "misc_RNA",
    "RNase_P_RNA": "misc_RNA",
    "antisense_RNA": "misc_RNA",
    "misc_RNA": "misc_RNA",
    "pseudogene": "pseudogene",
    "transcribed_pseudogene": "pseudogene",
}
GENE_TYPES = {"gene", "pseudogene"}


def parse_attrs(s):
    out = {}
    for part in s.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def clean(v):
    """Make a value safe to put back into a GFF attribute."""
    v = unquote(v or "")
    for ch in ";=,\t\n":
        v = v.replace(ch, "_")
    return v


def clean_text(v):
    """Make free text (e.g. a product description) safe for a TSV column."""
    return " ".join(unquote(v or "").split())


def read_gff(path):
    feats = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            feats.append(
                dict(
                    seqid=f[0], source=f[1], type=f[2], start=int(f[3]), end=int(f[4]),
                    strand=f[6], phase=f[7], attrs=parse_attrs(f[8]),
                )
            )
    return feats


def check_ids(feats, fai, log):
    fasta_ids = []
    with open(fai) as fh:
        for line in fh:
            fasta_ids.append(line.split("\t")[0])
    gff_ids = sorted({f["seqid"] for f in feats})
    shared = set(fasta_ids) & set(gff_ids)
    print(f"FASTA sequence IDs: {', '.join(fasta_ids)}", file=log)
    print(f"GFF sequence IDs:   {', '.join(gff_ids)}", file=log)
    if not shared:
        sys.exit(
            "ERROR: no sequence IDs shared between the reference FASTA and the GFF.\n"
            f"  FASTA: {fasta_ids[:5]}\n  GFF:   {gff_ids[:5]}\n"
            "Use a FASTA and GFF from the same assembly (e.g. both GCF_000017985.1)."
        )
    for missing in sorted(set(fasta_ids) - shared):
        print(f"WARNING: FASTA sequence {missing} has no annotation in the GFF", file=log)
    for missing in sorted(set(gff_ids) - shared):
        print(f"WARNING: GFF sequence {missing} is not in the FASTA; its features are skipped", file=log)
    return shared


def convert(gff_in, fai, gff_out, genes_out, log):
    feats = read_gff(gff_in)
    shared = check_ids(feats, fai, log)
    feats = [f for f in feats if f["seqid"] in shared]

    by_id = {}
    for f in feats:
        fid = f["attrs"].get("ID")
        if fid and fid not in by_id:
            by_id[fid] = f

    def gene_of(f):
        """Walk Parent links up to the gene-like feature."""
        seen = set()
        cur = f
        while True:
            parent = cur["attrs"].get("Parent", "").split(",")[0]
            if not parent or parent in seen or parent not in by_id:
                return None
            seen.add(parent)
            cur = by_id[parent]
            if cur["type"] in GENE_TYPES:
                return parent

    cds = defaultdict(list)       # gene ID -> CDS segments
    rnas = defaultdict(list)      # gene ID -> RNA features (tRNA, rRNA, ...)
    exons = defaultdict(list)     # RNA ID  -> exon segments
    for f in feats:
        if f["type"] in GENE_TYPES:
            continue
        gid = gene_of(f)
        if gid is None:
            continue
        if f["type"] == "CDS":
            cds[gid].append(f)
        elif f["type"] == "exon":
            exons[f["attrs"].get("Parent", "").split(",")[0]].append(f)
        elif f["type"] not in {"mRNA"}:
            rnas[gid].append(f)

    out = ["##gff-version 3"]
    gene_rows = []
    n = defaultdict(int)
    genes = sorted((f for f in feats if f["type"] in GENE_TYPES), key=lambda g: (g["seqid"], g["start"]))
    for g in genes:
        gid = g["attrs"]["ID"]
        a = g["attrs"]
        tag = clean(a.get("locus_tag") or gid.replace("gene-", ""))
        name = clean(a.get("gene") or a.get("Name") or tag)
        ncbi_biotype = a.get("gene_biotype", "pseudogene" if g["type"] == "pseudogene" else "")
        is_pseudo = g["type"] == "pseudogene" or a.get("pseudo") == "true" or "pseudogene" in ncbi_biotype

        segs = [c for c in cds.get(gid, []) if c["attrs"].get("pseudo") != "true"]
        if segs and not is_pseudo:
            # keep the segments of a single CDS ID (handles ribosomal-slippage genes)
            first_id = segs[0]["attrs"].get("ID")
            segs = sorted((c for c in segs if c["attrs"].get("ID") == first_id), key=lambda c: c["start"])
            biotype = "protein_coding"
        else:
            segs = []
            biotype = BIOTYPE_MAP.get(ncbi_biotype, "pseudogene" if is_pseudo else "misc_RNA")

        products = [r["attrs"].get("product") for r in cds.get(gid, []) + rnas.get(gid, [])]
        product = clean_text(next((p for p in products if p), "pseudogene" if is_pseudo else ""))
        gene_rows.append([g["seqid"], g["start"], g["end"], g["strand"], tag, name,
                          ncbi_biotype or biotype, product])

        col = lambda t, s, e, ph, attrs: "\t".join(
            [g["seqid"], "ncbi2csq", t, str(s), str(e), ".", g["strand"], ph, attrs]
        )
        gene_attr = f"ID=gene:{tag};Name={name};biotype={biotype};locus_tag={tag}"
        tx_attr = f"ID=transcript:{tag};Parent=gene:{tag};Name={name};biotype={biotype}"

        if segs:
            s, e = min(c["start"] for c in segs), max(c["end"] for c in segs)
            out.append(col("gene", s, e, ".", gene_attr))
            out.append(col("mRNA", s, e, ".", tx_attr))
            for c in segs:
                out.append(col("exon", c["start"], c["end"], ".", f"Parent=transcript:{tag}"))
            for c in segs:
                out.append(col("CDS", c["start"], c["end"], c["phase"], f"Parent=transcript:{tag}"))
            n["protein_coding"] += 1
        else:
            out.append(col("gene", g["start"], g["end"], ".", gene_attr))
            out.append(col("ncRNA" if not is_pseudo else "pseudogenic_transcript",
                           g["start"], g["end"], ".", tx_attr))
            rna_exons = [e for r in rnas.get(gid, []) for e in exons.get(r["attrs"].get("ID"), [])]
            for ex in sorted(rna_exons, key=lambda x: x["start"]) or [g]:
                out.append(col("exon", ex["start"], ex["end"], ".", f"Parent=transcript:{tag}"))
            n[biotype] += 1

    with open(gff_out, "w") as fh:
        fh.write("\n".join(out) + "\n")
    with open(genes_out, "w") as fh:
        fh.write("chrom\tstart\tend\tstrand\tlocus_tag\tgene\tbiotype\tproduct\n")
        for r in gene_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    print("Genes written by biotype: " + ", ".join(f"{k}={v}" for k, v in sorted(n.items())), file=log)
    if n["protein_coding"] == 0:
        sys.exit("ERROR: no protein-coding genes found in the GFF - is this an NCBI GFF3?")


if __name__ == "__main__":
    if "snakemake" in globals():
        with open(snakemake.log[0], "w") as log:  # noqa: F821
            convert(snakemake.input.gff, snakemake.input.fai,  # noqa: F821
                    snakemake.output.gff, snakemake.output.genes, log)  # noqa: F821
    else:
        convert(*sys.argv[1:5], sys.stderr)
