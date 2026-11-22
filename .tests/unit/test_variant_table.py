"""
Rule test for `variant_table` (written by hand; Snakemake's generator skipped it
because the rule's temporary output, *.query.tsv, no longer existed).

It reruns only the variant_table rule on a copy of your real inputs:
  .tests/unit/variant_table/data/      annotated VCF + gene table
and compares the result with the table from your real run:
  .tests/unit/variant_table/expected/  SRR030257.variant_table.tsv

On top of the "same as last time" check, it checks a few things that must
always be true of the table, whatever the data.
"""

import csv
import os
import shutil
import sys
import tempfile
from pathlib import Path
from subprocess import check_output

sys.path.insert(0, os.path.dirname(__file__))
from tsv_checker import TsvOutputChecker  # noqa: E402

RULE = "variant_table"
TARGET = "results/annotation/SRR030257.variant_table.tsv"

EXPECTED_COLUMNS = [
    "chrom", "pos", "ref", "alt", "type", "qual", "depth", "ref_reads", "alt_reads",
    "alt_fraction", "mapq", "region", "locus_tag", "gene", "biotype", "product",
    "consequence", "protein_change", "dna_change", "nearest_genes", "bcsq",
]


def run_rule(workdir, conda_prefix):
    """Run just this one rule in `workdir`, like the generated tests do."""
    cmd = [
        "python", "-m", "snakemake", TARGET,
        "--snakefile", "Snakefile",
        "-f", "--notemp", "--show-failed-logs", "-j1",
        "--target-files-omit-workdir-adjustment",
        "--allowed-rules", RULE,
        "--configfile", "config/config.yaml",
        "--directory", workdir,
    ]
    if shutil.which("conda"):  # use the rule's conda env (ecoli_variant) when conda is available
        cmd += ["--software-deployment-method", "conda"]
    check_output(cmd + conda_prefix)


def test_variant_table(conda_prefix):
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(f".tests/unit/{RULE}/config")
        data_path = Path(f".tests/unit/{RULE}/data")
        expected_path = Path(f".tests/unit/{RULE}/expected")

        shutil.copytree(config_path, workdir)
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        run_rule(workdir, conda_prefix)

        # 1. Same output as the saved run (row-by-row, with a readable diff)
        TsvOutputChecker(data_path, expected_path, workdir).check()

        # 2. Things that must always hold
        with open(workdir / TARGET) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        with open(workdir / TARGET) as fh:
            header = fh.readline().rstrip("\n").split("\t")

        assert header == EXPECTED_COLUMNS, f"Columns changed: {header}"
        assert rows, "The variant table is empty"

        # one row per variant in the annotated VCF
        n_vcf = sum(
            1 for line in check_output(["gzip", "-cd", workdir / "results/annotation/SRR030257.annotated.vcf.gz"],
                                       text=True).splitlines()
            if line and not line.startswith("#")
        )
        assert len(rows) == n_vcf, f"{len(rows)} table rows but {n_vcf} VCF records"

        for r in rows:
            where = f"{r['chrom']}:{r['pos']}"
            assert r["region"] in {"genic", "intergenic"}, where
            if r["region"] == "intergenic":
                assert r["consequence"] == "intergenic", where
                assert r["nearest_genes"], f"{where}: intergenic variant without flanking genes"
            else:
                assert r["gene"], f"{where}: genic variant without a gene name"
            if r["alt_fraction"]:
                assert 0.0 <= float(r["alt_fraction"]) <= 1.0, where
            assert len(r["protein_change"]) <= 43, f"{where}: protein_change not shortened"
