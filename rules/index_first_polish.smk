rule index_first_polish:
  """ index the polished assembly for short read alignment """
  input:
    reference = rules.first_polish.output.consensus_fasta
  conda:
    '../envs/index_first_polish.yaml'
  output:
    [rules.first_polish.output.consensus_fasta + suffix for suffix in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
  benchmark:
    "{asm}/benchmarks/index_first_polish.txt"
  log:
    stdout = "{asm}/logs/index_first_polish.out",
    stderr = "{asm}/logs/index_first_polish.err"
  shell:
    "bwa index " +
    "{input.reference} " +
    "2> {log.stderr} > {log.stdout}"
