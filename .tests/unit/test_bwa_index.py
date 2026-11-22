"""
Rule test code for unit testing of rules generated with Snakemake 9.27.0.
"""

import os
import sys
import shutil
import tempfile
from pathlib import Path
from subprocess import check_output

sys.path.insert(0, os.path.dirname(__file__))


def test_bwa_index(conda_prefix):

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/bwa_index/config")
        data_path = Path(".tests/unit/bwa_index/data")
        expected_path = Path(".tests/unit/bwa_index/expected")

        # Copy config to the temporary workdir.
        shutil.copytree(config_path, workdir)

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # Run the test job.
        check_output(
            [
                "python",
                "-m",
                "snakemake",
                "reference/ecoli_rel606.fasta.amb",
                "reference/ecoli_rel606.fasta.ann",
                "reference/ecoli_rel606.fasta.bwt",
                "reference/ecoli_rel606.fasta.pac",
                "reference/ecoli_rel606.fasta.sa",
                "--snakefile",
                "Snakefile",
                "-f",
                "--notemp",
                "--show-failed-logs",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--allowed-rules",
                "bwa_index",
                "--configfile",
                "config/config.yaml",
                "--directory",
                workdir,
            ]
            + conda_prefix
        )

        # Check the output byte by byte using cmp/zmp/bzcmp/xzcmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        import common
        common.OutputChecker(data_path, expected_path, workdir).check()
