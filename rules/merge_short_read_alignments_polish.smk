rule merge_short_read_alignments_polish:
  """ combine mappings to the polished assembly for all paired short readsets """
  input:
    alignments = lambda wildcards: expand(rules.map_short_reads_polish.output, asm=wildcards.asm, paired_short_readset_id=range(0, len(config['assemblies'][wildcards.asm]['paired_short_read_paths'])))
  conda:
    '../envs/merge_short_read_alignments_polish.yaml'
  output:
    merged_alignments = temp("{asm}/alignments/short_reads/{asm}.bam")
  benchmark:
    "{asm}/benchmarks/merge_short_read_alignments_polish.txt"
  log:
    stdout = "{asm}/logs/merge_short_read_alignments_polish.out",
    stderr = "{asm}/logs/merge_short_read_alignments_polish.err"
  shell:
    'samtools merge ' +
    "{output.merged_alignments} {input.alignments} "
    "2> {log.stderr} > {log.stdout}"
