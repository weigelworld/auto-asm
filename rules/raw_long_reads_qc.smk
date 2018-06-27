rule raw_long_reads_qc:
  """ report statistics and generate figures for the long readsets """
  input:
    reads = rules.convert_long_reads.output
  conda:
    '../envs/raw_long_reads_qc.yaml'
  output:
    read_length_dist_plot = report("{asm}/reports/raw_data/long_reads/{asm}.long_readset{long_readset_id}.length_dist.svg", caption="../report/read_length_dist_plot.rst", category="Raw Data"),
    summary_table = "{asm}/reports/raw_data/long_reads/{asm}.long_readset{long_readset_id}.summary.tsv"
  priority: 50
  benchmark:
    "{asm}/benchmarks/raw_long_reads_qc.readset{long_readset_id}.txt"
  log:
    stdout = "{asm}/logs/raw_long_reads_qc.readset{long_readset_id}.out",
    stderr = "{asm}/logs/raw_long_reads_qc.readset{long_readset_id}.err"
  script:
    '../scripts/raw_long_reads_qc.py'
