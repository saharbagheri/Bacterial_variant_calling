"""
Rule test for `summary` (written by hand; Snakemake's generator skipped it).

It reruns only the summary rule on a copy of your real inputs:
  .tests/unit/summary/data/      fastp JSON, flagstat, coverage, raw + filtered VCFs, variant table
and compares the result with the summary from your real run:
  .tests/unit/summary/expected/  summary.tsv

It also checks that the numbers in the summary agree with each other.
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

RULE = "summary"
TARGET = "results/summary/summary.tsv"


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


def test_summary(conda_prefix):
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

        # 2. The numbers must agree with each other
        with open(workdir / TARGET) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [r["sample"] for r in rows] == ["SRR030257"]

        for r in rows:
            num = {k: float(v) for k, v in r.items() if k not in ("sample", "consequences") and v != ""}
            # trimming can only remove reads
            assert num["reads_trimmed"] <= num["reads_raw"]
            # percentages are percentages
            for k, v in num.items():
                if k.startswith("pct_"):
                    assert 0 <= v <= 100, f"{k} = {v}"
            # coverage breadth can only shrink as the threshold rises
            assert num["pct_genome_ge1x"] >= num["pct_genome_ge10x"] >= num["pct_genome_ge30x"]
            # filtering can only remove variants, and the type counts add up
            assert num["filtered_total"] <= num["raw_total"]
            for stage in ("raw", "filtered"):
                assert (num[f"{stage}_snps"] + num[f"{stage}_indels"] + num[f"{stage}_other"]
                        == num[f"{stage}_total"]), f"{stage} counts do not add up"
            # every filtered variant has exactly one consequence entry
            n_csq = sum(int(part.split("=")[1]) for part in r["consequences"].split("; ") if part)
            assert n_csq == num["filtered_total"], (
                f"consequences cover {n_csq} variants, filtered_total is {num['filtered_total']:.0f}")
