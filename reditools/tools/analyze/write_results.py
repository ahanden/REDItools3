"""Write analysis results."""

import csv
from pathlib import Path
from typing import Callable, Iterable

from reditools.compiled_position import RTResult
from reditools.constants import (
    forward_strand_number,
    forward_strand_symbol,
    reverse_strand_number,
    reverse_strand_symbol,
    undetermined_strand_number,
    undetermined_strand_symbol,
)
from reditools.logger import Logger
from reditools.tools.analyze.rtchecks import RTChecks

_empty = "-"

class ResultWriter:
    """Helper class for writing REDItools analyze output."""

    def __init__(
        self,
        filters: RTChecks,
        logger: Callable,
        strand_numbers: bool = False,
        phred_breakdown: bool = False,
    ) -> None:
        """Construct new ResultWriter.

        Parameters
        ----------
        filters : RTChecks
            The quality control checks to apply.
        logger : Callable
            The logger function for debug messages.
        strand_numbers : bool
            If False, uses -, +, and * for the Strand column.
            If True, uses 0, 1, and 2 for the Strand column.
        phred_breakdown : bool
            If False, reports the mean PHRED score across all reads.
            If True, reports the mean PHRED score for reads that report
            each base.
        """
        self.filters = filters
        self.logger = logger

        if strand_numbers:
            self.strand_lookup = {
                forward_strand_symbol: forward_strand_number,
                reverse_strand_symbol: reverse_strand_number,
                undetermined_strand_symbol: undetermined_strand_number,
            }
        else:
            self.strand_lookup = {
                forward_strand_symbol: forward_strand_symbol,
                reverse_strand_symbol: reverse_strand_symbol,
                undetermined_strand_symbol: undetermined_strand_symbol,
            }

        if phred_breakdown:
            self._quality_fn = self._list_quality
        else:
            self._quality_fn = self._mean_quality

    def write_results(
        self,
        rtresults: Iterable[RTResult],
        filename: str,
    ) -> None:
        """Write analysis results to a file.

        Parameters
        ----------
        rtresults : Iterable[RTResult]
            The analysis results for each position.
        filename : str
            Where to save the results.
        """
        with Path(filename).open("w") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            for rt_result in rtresults:
                msg = self.filters.check(rt_result)
                if msg:
                    self.logger(Logger.debug_level, *msg)
                    continue
                if rt_result.variants:
                    variants = " ".join(rt_result.variants)
                else:
                    variants = _empty
                writer.writerow([
                    rt_result.contig,
                    rt_result.position + 1,
                    rt_result.reference,
                    self.strand_lookup[rt_result.strand],
                    len(rt_result),
                    self._quality_fn(rt_result),
                    list(rt_result),
                    variants,
                    f"{rt_result.edit_ratio:.2f}",
                    _empty, _empty, _empty, _empty, _empty,
                ])

    def _mean_quality(self, rt_result: RTResult) -> str:
        return f"{rt_result.mean_quality:.2f}"

    def _list_quality(self, rt_result: RTResult) -> str:
        insides = ", ".join([f"{_:.2f}" for _ in rt_result.quality_by_base])
        return f"[{insides}]"
