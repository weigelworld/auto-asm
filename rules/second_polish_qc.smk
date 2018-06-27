rule second_polish_qc:
  """ extract information from the second polish log files """
  input:
    rules.second_polish.output.consensus
  output:
    touch("{asm}/reports/assembly/second_polish.qc")
  priority: 50
