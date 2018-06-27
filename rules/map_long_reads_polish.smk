rule map_long_reads_polish:
  """ map raw long reads onto the initial assembly """
  input:
    reads = rules.link_long_reads.output,
    reference = rules.initial_assembly.output.contigs_fasta
  output:
    protected("{asm}/alignments/long_reads/readset{long_readset_id}.bam")
  benchmark:
    "{asm}/benchmarks/map_long_reads_polish.readset{long_readset_id}.txt"
  log:
    stdout = "{asm}/logs/map_long_reads_polish.readset{long_readset_id}.out",
    stderr = "{asm}/logs/map_long_reads_polish.readset{long_readset_id}.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'pbalign') +
    " {input.reads} {input.reference} {output} " +
    "--nproc {threads} " +
    "2> {log.stderr} > {log.stdout}"
