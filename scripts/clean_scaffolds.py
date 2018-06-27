from itertools import chain

from Bio import SeqIO


def clean_scaffolds(scaffolds_path, unplaced_contigs_path, output_path):
    scaffold_records = SeqIO.parse(scaffolds_path, 'fasta')
    unplaced_contigs = SeqIO.parse(unplaced_contigs_path, 'fasta')
    cleaned_records = []
    for record in chain(scaffold_records, unplaced_contigs):
        record.id = '_'.join(record.id.split('_')[1:])
        record.description = ''
        cleaned_records.append(record)

    SeqIO.write(cleaned_records, output_path, 'fasta')


if __name__ == '__main__':
    clean_scaffolds(
        snakemake.input['scaffolds'],
        snakemake.input['unplaced_contigs'],
        snakemake.output['final_assembly']
    )
