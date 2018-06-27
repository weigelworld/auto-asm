def get_short_read_file_paths(wildcards):
  """ fetch the paths to the paired-end short read files from the configuration """
  return config['assemblies'][wildcards.asm]['paired_short_read_paths'][int(wildcards.paired_short_readset_id)]

rule link_short_reads:
  """ link short read datasets to working directory """
  input:
    get_short_read_file_paths
  output:
    pair1 = "{asm}/raw_data/short_reads/readset{paired_short_readset_id}-raw-pair1.fastq.gz",
    pair2 = "{asm}/raw_data/short_reads/readset{paired_short_readset_id}-raw-pair2.fastq.gz"
  priority: 50
  shell:
    "ln -s {input[0]} {output.pair1} && " +
    "ln -s {input[1]} {output.pair2}"
