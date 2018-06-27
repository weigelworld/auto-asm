rule scaffolding_qc:
  """ Extract information from the second polish log files. """
  input:
    rules.clean_scaffolds.output.final_assembly
  output:
    touch("{asm}/reports/assembly/scaffolding.qc")
  priority: 50
