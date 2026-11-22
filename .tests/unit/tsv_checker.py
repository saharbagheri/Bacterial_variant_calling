"""
Output checker for the hand-written rule tests.

Same idea as common.OutputChecker (expected files must exist and match), but
.tsv files are compared line by line, so a failure shows exactly which rows
changed instead of only "files differ". Other files fall back to the
generator's byte-by-byte comparison.
"""

import difflib

import pytest

from common import OutputChecker, cmp_cmds


class TsvOutputChecker(OutputChecker):
    def compare_files(self, expected_file, generated_file, cmp_cmds=cmp_cmds):
        if expected_file.suffix != ".tsv":
            return super().compare_files(expected_file, generated_file, cmp_cmds)
        expected = expected_file.read_text().splitlines()
        generated = generated_file.read_text().splitlines()
        if expected != generated:
            diff = "\n".join(difflib.unified_diff(
                expected, generated, "expected", "generated", lineterm="", n=0))
            pytest.fail(f"{expected_file.name} differs from the saved run:\n{diff}", pytrace=False)
