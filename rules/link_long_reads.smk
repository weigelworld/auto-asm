def get_long_read_file_path(wildcards):
  """ fetch the paths to the long read files from the configuration """
  return config['assemblies'][wildcards.asm]['long_read_paths'][int(wildcards.long_readset_id)]

rule link_long_reads:
  """ link long read datasets to working directory """
  input:
    get_long_read_file_path
  output:
    "{asm}/raw_data/long_reads/readset{long_readset_id}.bam"
  priority: 50
  shell:
    "ln -s {input} {output} && " +
    "ln -s {input}.pbi {output}.pbi"
