from pathlib import Path

from snakemake.shell import shell


def calculate_fixed_gap_size(wildcards, config):
    """ calculate a fixed gap size for scaffolding """
    readset_N50s = []
    for i in range(len(config['assemblies'][wildcards.asm]['long_read_paths'])):
        with open(str("{asm}/reports/raw_data/long_reads/{asm}.long_readset{long_readset_id}.summary.tsv").format(asm=wildcards.asm, long_readset_id=i)) as long_read_summary_file:
            columns = long_read_summary_file.readline().split("\t")
            data = long_read_summary_file.readline().split("\t")
            long_read_summary = dict(zip(columns, data))
            readset_N50s.append(int(long_read_summary['N50']))
    highest_readset_N50 = max(readset_N50s)
    fixed_gap_size = int(0.75 * (highest_readset_N50 * 2))

    return fixed_gap_size


def scaffolding(output_scaffolds):
    output_scaffolds_path = Path(output_scaffolds)
    out_basename_path = output_scaffolds_path.parent / output_scaffolds_path.stem

    fixed_gap_size = calculate_fixed_gap_size(snakemake.wildcards, snakemake.config)
    print("Using fixed gap size of {} bp for scaffolding".format(fixed_gap_size))

    shell((
        "reveal finish --64 "
        "-m {params.kmer_size} "
        "--nproc {threads} "
        "--order {params.order} "
        "--fixedgapsize "
        "--gapsize {gapsize} "
        "-o {out_basename} "
        "{input.reference} {input.assembly} "
        "2> {log.stderr} > {log.stdout}").format(
            params=snakemake.params,
            threads=snakemake.threads,
            gapsize=fixed_gap_size,
            out_basename=str(out_basename_path),
            input=snakemake.input,
            log=snakemake.log
        ),
    )


if __name__ == '__main__':
    scaffolding(
        snakemake.output['scaffolds']
    )
