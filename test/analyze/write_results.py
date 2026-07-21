"""Test methods for analyze tool's write_results()."""

import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from reditools.compiled_position import CompiledPosition, RTResult
from reditools.logger import Logger
from reditools.tools.analyze.parse_args.parse_args import parse_args
from reditools.tools.analyze.rtchecks import RTChecks
from reditools.tools.analyze.write_results import write_results


class TestWriteResults(unittest.TestCase):
    """Test methods for analyze tool's write_results()."""

    def setUp(self) -> None:
        """Pre-flight setup."""
        self.rtresults = []

        cp = CompiledPosition(ref="A", contig="chr1", position=100)
        cp.add_base(40, "+", "A")
        cp.add_base(30, "+", "G")
        cp.add_base(30, "+", "T")
        self.rtresults.append(RTResult(cp, "+"))
        self.rtresults.append(RTResult(cp, "-"))
        self.rtresults.append(RTResult(cp, "*"))

        with NamedTemporaryFile(mode="w", dir=".", delete=False) as temp_file:
            self.fname = temp_file.name
        self.rtchecks = RTChecks(parse_args(["x.bam"]))
        self.logger = Logger(Logger.silent_level)

    def tearDown(self) -> None:
        """Post-checks cleanup."""
        Path(self.fname).unlink()

    def test_write_results(self) -> None:
        """Check basic functionality."""
        write_results(
            self.rtresults,
            self.fname,
            self.rtchecks,
            self.logger,
        )
        with Path(self.fname).open("r") as stream:
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "+", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "-", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "*", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )

    def test_write_results_numbers(self) -> None:
        """Check strand-number option."""
        write_results(
            self.rtresults,
            self.fname,
            self.rtchecks,
            self.logger,
            True,  # noqa: FBT003
        )
        with Path(self.fname).open("r") as stream:
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "1", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "0", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
            self.assertEqual(
                next(stream).strip().split("\t") ,
                [
                    "chr1", "101", "A", "2", "3", "33.33", "[1, 0, 1, 1]",
                    "AG AT", "0.50", "-", "-", "-", "-", "-",
                ],
            )
