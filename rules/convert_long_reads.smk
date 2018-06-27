rule convert_long_reads:
  """ generate a FASTA file from a raw PacBio BAM file """
  input:
    rules.link_long_reads.output
  params:
    prefix = lambda wildcards, output: str(Path(output['read_fasta']).parent / Path(output['read_fasta']).stem)
  output:
    read_fasta = protected("{asm}/raw_data/long_reads/readset{long_readset_id}.fasta")
  priority: 50
  benchmark:
    "{asm}/benchmarks/convert_long_reads.readset{long_readset_id}.txt"
  log:
    stdout = "{asm}/logs/convert_long_reads.readset{long_readset_id}.out",
    stderr = "{asm}/logs/convert_long_reads.readset{long_readset_id}.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'bam2fasta') +
    " --output {params.prefix} " +
    "-u {input} " +
    "2> {log.stderr} > {log.stdout}"
