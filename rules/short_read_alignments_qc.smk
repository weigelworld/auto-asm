rule short_read_alignments_qc:
  """ report statistics and generate figures for a short readset alignment """
  input:
    rules.map_short_reads_polish.output.alignments
  output:
    touch("{asm}/reports/alignments/short_reads/readset{paired_short_readset_id}.qc")
  priority: 50
