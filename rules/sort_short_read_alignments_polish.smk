rule sort_short_read_alignments_polish:
  """ sort the merged alignments of short reads to the polished assembly """
  input:
    rules.merge_short_read_alignments_polish.output.merged_alignments
  conda:
    '../envs/sort_short_read_alignments_polish.yaml'
  output:
    sorted_alignments = protected("{asm}/alignments/short_reads/{asm}.sorted.bam")
  benchmark:
    "{asm}/benchmarks/sort_short_read_alignments_polish.txt"
  log:
    stderr = "{asm}/logs/sort_short_read_alignments_polish.err"
  shell:
    'samtools sort ' +
    "--threads {threads} " +
    "{input} " +
    "2> {log.stderr} > {output.sorted_alignments}"
