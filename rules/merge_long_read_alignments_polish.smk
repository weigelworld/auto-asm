rule merge_long_read_alignments_polish:
  """ combine mappings to the initial assembly for all long readsets """
  input:
    lambda wildcards: expand(rules.map_long_reads_polish.output, asm=wildcards.asm, long_readset_id=range(0, len(config['assemblies'][wildcards.asm]['long_read_paths'])))
  output:
    merged_alignments = temp("{asm}/alignments/long_reads/{asm}.bam")
  benchmark:
    "{asm}/benchmarks/merge_long_read_alignments_polish.txt"
  log:
    stdout = "{asm}/logs/merge_long_read_alignments_polish.out",
    stderr = "{asm}/logs/merge_long_read_alignments_polish.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'samtools') +
    ' merge ' +
    "{output.merged_alignments} " +
    "{input} " +
    "2> {log.stderr} > {log.stdout}"
