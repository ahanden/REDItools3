from reditools import file_utils, reditools
from reditools.logger import Logger


def setup_rtools_dna(options):
    """
    Create a REDItools object for DNA processing.

    Parameters:
        options (namespace): Commandline arguments from argparse

    Returns:
        A configured REDItools object
    """
    rtools = reditools.REDItools()

    if options.debug:
        rtools.log_level = Logger.debug_level
    elif options.verbose:
        rtools.log_level = Logger.info_level

    if options.bed_file:
        for fname in options.bed_file:
            regions = file_utils.read_bed_file(fname)
            rtools.add_target_regions(regions)
    if options.exclude_regions:
        for fname in options.exclude_regions:
            regions = file_utils.read_bed_file(fname)
            rtools.add_exclude_regions(regions)
    if options.reference:
        rtools.add_reference(options.reference)

    if options.splicing_file:
        rtools.splice_positions = file_utils.load_splicing_file(
            options.splicing_file,
            options.splicing_span,
        )
        rtools.add_exclude_regions(regions)

    rtools.min_base_position = options.dna_min_base_position
    rtools.max_base_position = options.dna_max_base_position
    rtools.min_base_quality = options.dna_min_base_quality

    rtools.min_column_length = options.dna_min_read_depth
    rtools.min_edits = options.dna_min_edits
    rtools.min_edits_per_nucleotide = options.dna_min_edits_per_nucleotide
    rtools.max_alts = options.dna_max_editing_nucleotides

    rtools.strand = 0

    return rtools
