rule clean_scaffolds:
  """ rename and format scaffolds and merge with unplaced contigs """
  input:
    scaffolds = rules.scaffolding.output.scaffolds,
    unplaced_contigs = rules.scaffolding.output.unplaced_contigs
  conda:
    '../envs/clean_scaffolds.yaml'
  output:
    final_assembly = protected("{asm}/assembly/scaffolded/{asm}.final.fasta")
  script:
    '../scripts/clean_scaffolds.py'
