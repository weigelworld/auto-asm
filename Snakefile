from inspect import cleandoc
from pathlib import Path

import snakemake.utils


##############
# PARAMETERS #
##############

snakemake.utils.min_version('5.1')

snakemake.utils.validate(config, 'schemas/config.schema.yaml')

report: "report/workflow.rst"

localrules:
  all,
  link_long_reads,
  raw_long_reads_qc,
  initial_assembly_qc,
  first_polish_qc,
  link_short_reads,
  second_polish_qc,
  rename_contigs,
  clean_scaffolds,
  scaffolding_qc,
  final_assembly_report

onstart:
  # create cluster log directories (SGE doesn't create them automatically)
  for asm in config['assemblies']:
    snakemake.utils.makedirs("{asm}/logs/cluster".format(asm=asm))

onsuccess:
  onsuccess_message = cleandoc("""
  auto-asm has finished successfully.
  Final contigs are in 'assembly/' under each assembly's subfolder in the output directory.
  Final scaffolds, if a reference was provided, are in 'assembly/scaffolded/'.
  All reports are in 'reports/' under each assembly subfolder.
  """)
  print(onsuccess_message)

onerror:
  onerror_message = cleandoc("""
  auto-asm has failed.
  Check the following log files for more information:
    - snakemake log at {log}
    - job stdout/stderr in 'logs/' under each assembly subfolder in the output directory
    - cluster stdout/stderr in 'logs/cluster' under each assembly subfolder
  """.format(log=log))
  print(onerror_message)


###########
# TARGETS #
###########

rule all:
  input:
    ["{asm}/reports/final.touch".format(asm=asm) for asm in config['assemblies']]


#########
# RULES #
#########

include: 'rules/link_long_reads.smk'
include: 'rules/convert_long_reads.smk'
include: 'rules/raw_long_reads_qc.smk'
include: 'rules/initial_assembly.smk'
include: 'rules/initial_assembly_qc.smk'
include: 'rules/map_long_reads_polish.smk'
include: 'rules/long_read_alignments_qc.smk'
include: 'rules/merge_long_read_alignments_polish.smk'
include: 'rules/sort_long_read_alignments_polish.smk'
include: 'rules/index_long_read_alignments_polish.smk'
include: 'rules/index_initial_assembly.smk'
include: 'rules/first_polish.smk'
include: 'rules/first_polish_qc.smk'
include: 'rules/index_first_polish.smk'
include: 'rules/link_short_reads.smk'
include: 'rules/raw_short_reads_qc.smk'
include: 'rules/trim_short_reads.smk'
include: 'rules/trimmed_short_reads_qc.smk'
include: 'rules/map_short_reads_polish.smk'
include: 'rules/short_read_alignments_qc.smk'
include: 'rules/merge_short_read_alignments_polish.smk'
include: 'rules/sort_short_read_alignments_polish.smk'
include: 'rules/index_short_read_alignments_polish.smk'
include: 'rules/second_polish.smk'
include: 'rules/second_polish_qc.smk'
include: 'rules/rename_contigs.smk'
include: 'rules/scaffolding.smk'
include: 'rules/clean_scaffolds.smk'
include: 'rules/scaffolding_qc.smk'
include: 'rules/assembly_comparison.smk'
include: 'rules/final_assembly_report.smk'


########################
# RULE POST-PROCESSING #
########################

for rule_obj in workflow.rules:
  if rule_obj.name in config['rules']:
    if 'resources' in config['rules'][rule_obj.name]:
      # handle resource assignments
      for key in config['rules'][rule_obj.name]['resources']:
        rule_obj.resources[key] = config['rules'][rule_obj.name]['resources'][key]
    if 'params' in config['rules'][rule_obj.name]:
      # handle parameter assignments
      for key in config['rules'][rule_obj.name]['params']:
        rule_obj.set_params(**config['rules'][rule_obj.name]['params'])

  if 'cores' not in rule_obj.resources:
    # provide a default core resource of 1
    rule_obj.resources['cores'] = 1
  if 'mem_mb' not in rule_obj.resources:
    # provide a default memory resource of 1000 MB
    rule_obj.resources['mem_mb'] = 1000
  # make sure threads are also set properly
  rule_obj.resources['_cores'] = rule_obj.resources['cores']
  # calculate memory per core as a resource
  rule_obj.resources['mem_mb_per_core'] = round(
    rule_obj.resources['mem_mb'] / rule_obj.resources['_cores']
  )
