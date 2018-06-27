rule index_long_read_alignments_polish:
  """ index the merged long read mappings """
  input:
    rules.sort_long_read_alignments_polish.output.sorted_alignments
  output:
    index = rules.sort_long_read_alignments_polish.output.sorted_alignments + ".pbi"
  benchmark:
    "{asm}/benchmarks/index_long_read_alignments_polish.txt"
  log:
    stdout = "{asm}/logs/index_long_read_alignments_polish.out",
    stderr = "{asm}/logs/index_long_read_alignments_polish.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'pbindex') +
    " {input} " +
    "2> {log.stderr} > {log.stdout}"
