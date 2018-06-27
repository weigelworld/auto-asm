"""
utils module

common miscellaneous and utility functions
"""

from io import StringIO
from statistics import mean
from typing import List, Union


GENOME_SIZE_SUFFIX_FACTORS = {
    'k': 1000,
    'm': 1000000,
    'g': 1000000000
}


def NXX(lengths: List[int], mass_percentage: Union[int, float]) -> int:
    """

    :param lengths: a list of lengths of contigs
    :param mass_percentage: the percentage cutoff (e.g. N50 = 50)
    :return:
    """
    if 0 > mass_percentage > 100:
        raise ValueError("NXX cutoff is outside valid range (0 <= cutoff <= 100)")

    total = sum(lengths)
    length_cutoff = (mass_percentage / 100) * total

    working_sum = 0
    for length in sorted(lengths, reverse=True):
        working_sum += length
        if working_sum >= length_cutoff:
            return length


def N50(lengths: List[int]) -> int:
    return NXX(lengths, 50)


def N90(lengths: List[int]) -> int:
    return NXX(lengths, 90)


def expand_genome_size(genome_size: str) -> int:
    suffix = genome_size[-1].lower()
    try:
        prefix = float(genome_size[:-1])
    except ValueError:
        raise ValueError("Unrecognized genome size: {}".format(genome_size))

    try:
        expanded_prefix = round(prefix * GENOME_SIZE_SUFFIX_FACTORS[suffix])
    except KeyError:
        raise ValueError("Unrecognized genome size: {}".format(genome_size))

    return expanded_prefix


def readset_summary_table(read_lengths: List[int], assembly_id: str, readset_id: str, genome_size: str) -> StringIO:
    output_buffer = StringIO()
    output_buffer.write("assembly.id\treadset.id\ttotal.reads\ttotal.length\taverage.length\tN50\tN90\tcoverage\n")

    total_length = sum(read_lengths)
    expanded_genome_size = expand_genome_size(genome_size)
    output_buffer.write(
        "{assembly_id}\t{readset_id}\t{total_reads}\t{total_length}\t{average_length}\t{N50}\t{N90}\t{coverage}\n".format(
            assembly_id=assembly_id,
            readset_id=readset_id,
            total_reads=len(read_lengths),
            total_length=total_length,
            average_length=round(mean(read_lengths)),
            N50=round(N50(read_lengths)),
            N90=round(N90(read_lengths)),
            coverage=round((total_length / expanded_genome_size))
        )
    )

    # reset the buffer position for reading
    output_buffer.seek(0)

    return output_buffer
