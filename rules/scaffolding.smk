rule scaffolding:
  """ Scaffold the contigs based on a reference genome """
  input:
    assembly = rules.rename_contigs.output.contigs,
    reference = lambda wildcards: config['assemblies'][wildcards.asm]['reference']
  conda:
    '../envs/scaffolding.yaml'
  output:
    scaffolds = protected("{asm}/assembly/scaffolded/{asm}.fasta"),
    unplaced_contigs = protected("{asm}/assembly/scaffolded/{asm}.unplaced.fasta")
  log:
    stdout = "{asm}/logs/scaffolding.out",
    stderr = "{asm}/logs/scaffolding.err"
  script:
    '../scripts/scaffolding.py'
