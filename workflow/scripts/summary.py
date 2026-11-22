"""
Collect the numbers the guide asks you to record into one table
(results/summary/summary.tsv), one row per sample:
reads before/after trimming, % Q30, mapping rate, coverage,
raw / filtered variant counts (SNPs, indels) and csq consequence counts.
"""

import gzip
import json
import re
from collections import Counter


def fastp_stats(path):
    with open(path) as fh:
        s = json.load(fh)["summary"]
    b, a = s["before_filtering"], s["after_filtering"]
    return {
        "reads_raw": b["total_reads"],
        "reads_trimmed": a["total_reads"],
        "pct_reads_kept": round(100 * a["total_reads"] / b["total_reads"], 2) if b["total_reads"] else 0,
        "pct_q30_raw": round(100 * b["q30_rate"], 2),
        "pct_q30_trimmed": round(100 * a["q30_rate"], 2),
    }


def flagstat_stats(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"(\d+) \+ \d+ (.*)", line)
            if not m:
                continue
            n, what = int(m.group(1)), m.group(2)
            pct = re.search(r"\(([\d.]+)%", what)
            if what.startswith("in total"):
                out["alignments_total"] = n
            elif what.startswith("primary mapped"):
                out["primary_mapped"] = n
                out["pct_primary_mapped"] = float(pct.group(1)) if pct else ""
            elif what.startswith("mapped") and "pct_mapped" not in out:
                out["pct_mapped"] = float(pct.group(1)) if pct else ""
            elif what.startswith("properly paired"):
                out["pct_properly_paired"] = float(pct.group(1)) if pct else ""
    return out


def coverage_stats(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        values = fh.readline().rstrip("\n").split("\t")
    d = dict(zip(header, values))
    d.pop("sample", None)
    d.pop("genome_positions", None)
    return d


def vcf_counts(path, prefix):
    c = Counter()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            ref, alts = f[3], f[4].split(",")
            c["total"] += 1
            if all(len(a) == len(ref) == 1 for a in alts):
                c["snps"] += 1
            elif any(len(a) != len(ref) for a in alts):
                c["indels"] += 1
            else:
                c["other"] += 1
    return {f"{prefix}_{k}": c[k] for k in ("total", "snps", "indels", "other")}


def consequence_counts(path):
    c = Counter()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index("consequence")
        for line in fh:
            csq = line.rstrip("\n").split("\t")[idx] or "unannotated"
            c[csq] += 1
    return {"consequences": "; ".join(f"{k}={v}" for k, v in c.most_common())}


rows = []
inp = snakemake.input  # noqa: F821
for i, sample in enumerate(snakemake.params.samples):  # noqa: F821
    row = {"sample": sample}
    row.update(fastp_stats(inp.fastp[i]))
    row.update(flagstat_stats(inp.flagstat[i]))
    row.update(coverage_stats(inp.coverage[i]))
    row.update(vcf_counts(inp.raw[i], "raw"))
    row.update(vcf_counts(inp.filtered[i], "filtered"))
    row.update(consequence_counts(inp.table[i]))
    rows.append(row)

cols = []
for r in rows:
    cols += [k for k in r if k not in cols]
with open(snakemake.output[0], "w") as fh:  # noqa: F821
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(k, "")) for k in cols) + "\n")
