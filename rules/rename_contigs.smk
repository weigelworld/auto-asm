def get_rename_contigs_input(wildcards):
  if 'paired_short_read_paths' in config['assemblies'][wildcards.asm]:
    return {
      'contigs': rules.second_polish.output.consensus
    }
  else:
    return {
      'contigs': rules.first_polish.output.consensus_fasta
    }

rule rename_contigs:
  """ standarize contig names and generate a mapping between new and old names """
  input:
    unpack(get_rename_contigs_input)
  conda:
    '../envs/rename_contigs.yaml'
  output:
    contigs = "{asm}/assembly/{asm}.contigs.fasta",
    map = "{asm}/assembly/{asm}.contigs.map"
  script:
    '../scripts/rename_contigs.py'
