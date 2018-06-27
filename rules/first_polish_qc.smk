rule first_polish_qc:
  """ extract information from the first polish log files """
  input:
    rules.first_polish.output.consensus_fasta
  output:
    touch("{asm}/reports/assembly/first_polish.qc")
  priority: 50
