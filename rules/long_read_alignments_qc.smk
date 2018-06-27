rule long_read_alignments_qc:
  """ report statistics and generate figures for a long readset alignment """
  input:
    rules.map_long_reads_polish.output
  output:
    touch("{asm}/reports/alignments/long_reads/readset{long_readset_id}.qc")
  priority: 50
