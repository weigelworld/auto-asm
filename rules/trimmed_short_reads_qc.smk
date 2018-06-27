rule trimmed_short_reads_qc:
  """ report statistics and generate figures for the trimmed short readsets """
  input:
    rules.trim_short_reads.output.trimmed_short_reads
  conda:
    '../envs/trimmed_short_reads_qc.yaml'
  params:
    out_dir = lambda wildcards, output: str(Path(output['reports'][0]).parent)
  log:
    stdout = "{asm}/logs/trimmed_short_reads_qc.readset{paired_short_readset_id}.out",
    stderr = "{asm}/logs/trimmed_short_reads_qc.readset{paired_short_readset_id}.err"
  output:
    reports = (
      protected("{asm}/reports/raw_data/short_reads/readset{paired_short_readset_id}-trimmed-pair1_fastqc.html"),
      protected("{asm}/reports/raw_data/short_reads/readset{paired_short_readset_id}-trimmed-pair2_fastqc.html")
    )
  priority: 50
  benchmark:
    "{asm}/benchmarks/trimmed_short_reads_qc.readset{paired_short_readset_id}.txt"
  shell:
    'fastqc ' +
    '--extract ' +
    "--threads {threads} " +
    "--outdir {params.out_dir} " +
    "{input} " +
    "2> {log.stderr} > {log.stdout}"
