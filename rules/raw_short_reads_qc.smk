rule raw_short_reads_qc:
  """ report statistics and generate figures for the raw short readsets """
  input:
    rules.link_short_reads.output
  conda:
    '../envs/raw_short_reads_qc.yaml'
  params:
    out_dir = lambda wildcards, output: str(Path(output['reports'][0]).parent)
  log:
    stdout = "{asm}/logs/raw_short_reads_qc.readset{paired_short_readset_id}.out",
    stderr = "{asm}/logs/raw_short_reads_qc.readset{paired_short_readset_id}.err"
  output:
    reports = (
      protected("{asm}/reports/raw_data/short_reads/readset{paired_short_readset_id}-raw-pair1_fastqc.html"),
      protected("{asm}/reports/raw_data/short_reads/readset{paired_short_readset_id}-raw-pair2_fastqc.html")
    )
  priority: 50
  benchmark:
    "{asm}/benchmarks/raw_short_reads_qc.readset{paired_short_readset_id}.txt"
  shell:
    'fastqc ' +
    '--extract ' +
    "--threads {threads} " +
    "--outdir {params.out_dir} " +
    "{input}  " +
    "2> {log.stderr} > {log.stdout}"
