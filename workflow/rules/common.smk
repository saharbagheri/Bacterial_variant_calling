import csv
import re
import shlex


def load_samples(path):
    """Read the tab-separated sample sheet.

    Columns: sample, sra, fq1, fq2
      - give `sra` to download reads from NCBI SRA, or
      - give `fq1` + `fq2` (gzipped FASTQ) to use reads already on disk.
    """
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    samples = {}
    for row in csv.DictReader(lines, delimiter="\t"):
        row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
        name = row.get("sample", "")
        if not name:
            continue
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            raise ValueError(f"Sample name '{name}' may only contain letters, digits, '_', '.', '-'")
        has_local = bool(row.get("fq1") and row.get("fq2"))
        if not has_local and not row.get("sra"):
            raise ValueError(f"Sample '{name}': give either an SRA accession or both fq1 and fq2")
        row["source"] = "local" if has_local else "sra"
        samples[name] = row
    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


SAMPLE_INFO = load_samples(config["samples"])
SAMPLES = list(SAMPLE_INFO)
SRA_SAMPLES = [s for s, r in SAMPLE_INFO.items() if r["source"] == "sra"]
LOCAL_SAMPLES = [s for s, r in SAMPLE_INFO.items() if r["source"] == "local"]

REF_NAME = config["reference"]["name"]
REF = f"reference/{REF_NAME}.fasta"
GFF = f"reference/{REF_NAME}.gff"
CSQ_GFF = f"reference/{REF_NAME}.csq.gff3"

REF_FAI = f"{REF}.fai"
BWA_INDEX = [f"{REF}.{ext}" for ext in ("amb", "ann", "bwt", "pac", "sa")]

THREADS = int(config.get("threads", 4))


def fetch_cmd(source):
    """Shell snippet that writes a plain or gzipped source to stdout."""
    if re.match(r"^(https?|ftp)://", source):
        download = f"curl -fsSL --retry 3 '{source}'"
        return f"{download} | gzip -cdf" if source.split("?", 1)[0].lower().endswith(".gz") else download
    if source.endswith(".gz"):
        return f"gzip -cdf {shlex.quote(str(source))}"
    return f"cat {shlex.quote(str(source))}"


def reference_input(key):
    """Expose local reference files to Apptainer; URLs are fetched by the rule."""
    source = config["reference"][key]
    if re.match(r"^(https?|ftp)://", source):
        return []
    return workflow.source_path(source)


def regex_for(names):
    """Regex matching exactly these sample names (or nothing, if empty)."""
    return "|".join(re.escape(n) for n in names) if names else "(?!x)x"


wildcard_constraints:
    sample=regex_for(SAMPLES),
    read="R1|R2",
