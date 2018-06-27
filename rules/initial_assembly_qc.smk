rule initial_assembly_qc:
  """ extract information from the de novo assembler log files """
  input:
    rules.initial_assembly.output.contigs_fasta
  output:
    touch("{asm}/reports/assembly/initial.qc")
  priority: 50
