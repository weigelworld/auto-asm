from Bio import SeqIO

import auto_asm.plots
import auto_asm.utils


def raw_long_reads_qc(read_file_path, read_distribution_plot_path, read_statistics_path):
    reads = SeqIO.parse(str(read_file_path), 'fasta')
    read_lengths = [len(read) for read in reads]

    # plot read length distribution
    read_dist_figure = auto_asm.plots.read_length_histogram(read_lengths, assembly_id=snakemake.wildcards['asm'],
                                                            readset_id=snakemake.wildcards['long_readset_id'])
    read_dist_figure.savefig(str(read_distribution_plot_path))

    # calculate summary statistics
    with open(str(read_statistics_path), 'w') as read_statistics_file:
        read_statistics_buffer = auto_asm.utils.readset_summary_table(
            read_lengths,
            snakemake.wildcards['asm'],
            snakemake.wildcards['long_readset_id'],
            snakemake.config['assemblies'][snakemake.wildcards['asm']]['genome_size']
        )
        read_statistics_file.write(read_statistics_buffer.read())


if __name__ == '__main__':
    raw_long_reads_qc(
        snakemake.input['reads'],
        snakemake.output['read_length_dist_plot'],
        snakemake.output['summary_table']
    )
