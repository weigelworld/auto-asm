rule sort_long_read_alignments_polish:
  """ sort the merged alignments of long reads to the initial assembly """
  input:
    rules.merge_long_read_alignments_polish.output.merged_alignments
  output:
    sorted_alignments = protected("{asm}/alignments/long_reads/{asm}.sorted.bam")
  benchmark:
    "{asm}/benchmarks/sort_long_read_alignments_polish.txt"
  log:
    stderr = "{asm}/logs/sort_long_read_alignments_polish.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'samtools') +
    ' sort ' +
    "--threads {threads} " +
    "{input} " +
    "2> {log.stderr} > {output.sorted_alignments}"
