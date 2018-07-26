def get_assembly_comparison_input(wildcards):
  """ determine which assemblies to perform based on provided data """
  requested_input = {
    'initial': rules.initial_assembly.output.contigs_fasta,
    'first_polish': rules.first_polish.output.consensus_fasta,
  }

  if 'paired_short_read_paths' in config['assemblies'][wildcards.asm]:
    # paired-end short read data provided, request second polish
    requested_input['second_polish'] = rules.second_polish.output.consensus
  if 'reference' in config['assemblies'][wildcards.asm]:
    # reference genome provided, request scaffolding
    requested_input['scaffolded'] = rules.clean_scaffolds.output.final_assembly

  return requested_input

def get_assembly_labels(wildcards, input):
  """ generate a label list depending on which assemblies are being generated """
  labels = ['Initial', 'First Polish']

  if 'second_polish' in input.keys():
    labels.append('Second Polish')
  if 'scaffolded' in input.keys():
    labels.append('Scaffolded')

  return ', '.join(labels)

rule assembly_comparison:
  """ compare statistics across assembly steps and to a reference """
  input:
    unpack(get_assembly_comparison_input)
  conda:
    '../envs/assembly_comparison.yaml'
  params:
    out_dir = lambda wildcards, output: str(Path(output['report_tsv']).parent),
    labels = get_assembly_labels,
    eukaryote = lambda wildcards: '--eukaryote' if config['assemblies'][wildcards.asm]['eukaryotic'] else ''
  output:
    report_tsv = protected("{asm}/reports/assembly/comparison/transposed_report.tsv")
  benchmark:
    "{asm}/benchmarks/assembly_comparison.txt"
  log:
    stdout = "{asm}/logs/assembly_comparison.out",
    stderr = "{asm}/logs/assembly_comparison.err"
  shell:
    'quast ' +
    "-o {params.out_dir} " +
    ("-R {} ".format(config['assemblies'][wildcards.asm]['reference']) if hasattr(input, 'scaffolded') else '') +
    "--labels {params.labels:q} " +
    "--threads {threads} " +
    "{params.eukaryote} " +
    "{input} " +
    "2> {log.stderr} > {log.stdout}"
