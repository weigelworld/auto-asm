from Bio import SeqIO


def rename_contigs(contigs_path, map_output_path, renamed_output_path):
    contigs = SeqIO.parse(contigs_path, 'fasta')
    renamed_records = []
    with open(map_output_path, 'w') as map_file:
        map_file.write("old.id\tnew.id\n")
        for i, record in enumerate(contigs):
            new_id = "CTG_{}".format(i)
            map_file.write("{old}\t{new}\n".format(old=record.id, new=new_id))
            record.id = new_id
            record.description = ''
            renamed_records.append(record)

    SeqIO.write(renamed_records, renamed_output_path, 'fasta')


if __name__ == '__main__':
    rename_contigs(
        snakemake.input['contigs'],
        snakemake.output['map'],
        snakemake.output['contigs']
    )
