rule trim_short_reads:
  """ trim and filter short reads to desired qualities """
  input:
    rules.link_short_reads.output
  params:
    prefix = lambda wildcards, input: str(Path(input[0]).parent / Path(input[0]).stem.split('-')[0]),
    mode='pe'
  conda:
    '../envs/trim_short_reads.yaml'
  log:
    default = "{asm}/raw_data/short_reads/readset{paired_short_readset_id}-trimmed.log",
    stdout = "{asm}/logs/trim_short_reads.readset{paired_short_readset_id}.out",
    stderr = "{asm}/logs/trim_short_reads.readset{paired_short_readset_id}.err"
  output:
    trimmed_short_reads = (
      protected('{asm}/raw_data/short_reads/readset{paired_short_readset_id}-trimmed-pair1.fastq.gz'),
      protected('{asm}/raw_data/short_reads/readset{paired_short_readset_id}-trimmed-pair2.fastq.gz'))
  priority: 50
  benchmark:
    "{asm}/benchmarks/trim_short_reads.readset{paired_short_readset_id}.txt"
  shell:
    "skewer -z " +
    "-f {params.FASTQ_qual_format} " +
    "-t {threads} " +
    "-q {params.min_3prime_end_qual} " +
    "-Q {params.min_mean_pretrim_read_qual} " +
    "-l {params.min_trimmed_len} " +
    "-o {params.prefix} " +
    "-m {params.mode} " +
    "{input} " +
    "2> {log.stderr} > {log.stdout}"
