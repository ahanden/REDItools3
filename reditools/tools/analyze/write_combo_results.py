import csv
from tempfile import NamedTemporaryFile


def write_combo_results(rna_rtools, dna_rtools, rna_sam_manager,
                        dna_sam_manager, file_name, region,
                        output_format, temp_dir):
    """
    Write the results from a REDItools analysis to a temporary file.

    Parameters:
        rna_rtools (REDItools): REDItools instance
        dna_rtools (REDItools): REDItools instance
        rna_sam_manager (AlignmentManager): Source of reads
        dna_sam_manager (AlignmentManager): Source of reads
        file_name (string): Input file name for analysis
        region: Region to analyze
        output_format (dict): keyword arguments for csv.writer constructor.

    Returns:
        string: Name of the temporary file.
    """
    with NamedTemporaryFile(mode='w', delete=False, dir=temp_dir) as stream:
        writer = csv.writer(stream, **output_format)
        dna_results_iter = dna_rtools.analyze(dna_sam_manager, region)
        dna_result = next(dna_results_iter, None)
        for rna_result in rna_rtools.analyze(rna_sam_manager, region):
            variants = rna_result.variants
            row = [
                rna_result.contig,
                rna_result.position,
                rna_result.reference,
                rna_result.strand,
                rna_result.depth,
                f'{rna_result.mean_quality:.2f}',
                rna_result.per_base_depth,
                ' '.join(sorted(variants)) if variants else '-',
                f'{rna_result.edit_ratio:.2f}',
            ]
            while dna_result is not None and dna_result.position < rna_result.position:
                dna_result = next(dna_results_iter, None)
            if dna_result is not None and dna_result.position == rna_result.position:
                dna_variants = dna_result.variants
                row.extend([
                    dna_result.depth,
                    f'{dna_result.mean_quality:.2f}',
                    dna_result.per_base_depth,
                    ' '.join(sorted(dna_variants)) if dna_variants else '-',
                    f'{dna_result.edit_ratio:.2f}',
                ])
            else:
                row.extend(['-', '-', '-', '-', '-'])
            writer.writerow(row)
        return stream.name
