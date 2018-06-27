rule index_short_read_alignments_polish:
  """ index the sorted short read alignments """
  input:
    rules.sort_short_read_alignments_polish.output.sorted_alignments
  output:
    index = rules.sort_short_read_alignments_polish.output.sorted_alignments + ".bai"
  benchmark:
    "{asm}/benchmarks/index_short_read_alignments_polish.txt"
  log:
    stdout = "{asm}/logs/index_short_read_alignments_polish.out",
    stderr = "{asm}/logs/index_short_read_alignments_polish.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'samtools') +
    " index {input} " +
    "2> {log.stderr} > {log.stdout}"
