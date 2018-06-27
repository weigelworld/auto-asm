rule initial_assembly:
  """ assemble reads into contigs """
  input:
    lambda wildcards: ["{asm}/raw_data/long_reads/readset{long_readset_id}.fasta".format(long_readset_id=i, asm=wildcards.asm) for i in range(0, len(config['assemblies'][wildcards.asm]['long_read_paths']))]
  conda:
    "../envs/initial_assembly.yaml"
  params:
    genomeSize = lambda wildcards: config['assemblies'][wildcards.asm]['genome_size'],
    outDir = lambda wildcards, output: str(Path(output['contigs_fasta']).parent),
    useGrid='false',
    ovsMethod='sequential'
  output:
    contigs_fasta = protected("{asm}/assembly/initial/{asm}.contigs.fasta")
  benchmark:
    "{asm}/benchmarks/initial_assembly.txt"
  log:
    stdout = "{asm}/logs/initial_assembly.out",
    stderr = "{asm}/logs/initial_assembly.err"
  shell:
    'canu ' +
    "-pacbio-raw {input} " +
    "genomeSize={params.genomeSize} " +
    "-p {wildcards.asm} " +
    "-d {params.outDir} " +
    "useGrid={params.useGrid} " +
    "ovsMethod={params.ovsMethod} " +
    ("maxThreads={} ".format(config['rules']['initial_assembly']['resources']['cores']) if 'initial_assembly' in config['rules'] and 'resources' in config['rules']['initial_assembly'] and 'cores' in config['rules']['initial_assembly']['resources'] else '') +
    ("maxMemory={}m ".format(config['rules']['initial_assembly']['resources']['mem_mb']) if 'initial_assembly' in config['rules'] and 'resources' in config['rules']['initial_assembly'] and 'mem_mb' in config['rules']['initial_assembly']['resources'] else '') +
    "2> {log.stderr} > {log.stdout}"
