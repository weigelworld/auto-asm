rule map_short_reads_polish:
  """ map short reads to polished assembly """
  input:
    short_reads = rules.trim_short_reads.output.trimmed_short_reads,
    reference = rules.first_polish.output.consensus_fasta,
    reference_index = rules.index_first_polish.output
  conda:
    "../envs/map_short_reads_polish.yaml"
  output:
    alignments = temp("{asm}/alignments/short_reads/readset{paired_short_readset_id}.bam")
  benchmark:
    "{asm}/benchmarks/map_short_reads_polish.readset{paired_short_readset_id}.txt"
  log:
    stderr = "{asm}/logs/map_short_reads_polish.readset{paired_short_readset_id}.err"
  shell:
    "bwa mem " +
    "-t {threads} " +
    "{input.reference} {input.short_reads} " +
    "2> {log.stderr} | " +
    "samtools view -u " +
    "2> {log.stderr} > {output.alignments}"
