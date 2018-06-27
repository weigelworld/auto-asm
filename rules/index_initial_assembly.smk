rule index_initial_assembly:
  """ index the initial assembly contigs """
  input:
    rules.initial_assembly.output.contigs_fasta
  output:
    index = rules.initial_assembly.output.contigs_fasta + '.fai'
  benchmark:
    "{asm}/benchmarks/index_initial_assembly.txt"
  log:
    stdout = "{asm}/logs/index_initial_assembly.out",
    stderr = "{asm}/logs/index_initial_assembly.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'samtools') +
    " faidx {input} " +
    "2> {log.stderr} > {log.stdout}"
