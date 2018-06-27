def get_final_assembly_report_input(wildcards):
  """ Determine the correct inputs for the final qc report """
  input_files = {}
  input_files['long_reads_stats'] = expand(rules.raw_long_reads_qc.output.summary_table, asm=wildcards.asm, long_readset_id=range(0, len(config['assemblies'][wildcards.asm]['long_read_paths'])))
  input_files['long_reads_length_dist_plots'] = expand(rules.raw_long_reads_qc.output.read_length_dist_plot, asm=wildcards.asm, long_readset_id=range(0, len(config['assemblies'][wildcards.asm]['long_read_paths'])))
  input_files['initial_assembly_qc'] = str(rules.initial_assembly_qc.output).format(asm=wildcards.asm)
  input_files['long_read_alignments_qc'] = expand(rules.long_read_alignments_qc.output, asm=wildcards.asm, long_readset_id=range(0, len(config['assemblies'][wildcards.asm]['long_read_paths'])))
  input_files['first_polish_qc'] = str(rules.first_polish_qc.output).format(asm=wildcards.asm)
  if 'paired_short_read_paths' in config['assemblies'][wildcards.asm]:
    # whole genome short reads provided, go through second polish
    input_files['raw_short_reads_qc'] = [s for e in [expand(output, asm=wildcards.asm, paired_short_readset_id=range(0, len(config['assemblies'][wildcards.asm]['paired_short_read_paths']))) for output in rules.raw_short_reads_qc.output] for s in e]
    input_files['trimmed_short_reads_qc'] = [s for e in [expand(output, asm=wildcards.asm, paired_short_readset_id=range(0, len(config['assemblies'][wildcards.asm]['paired_short_read_paths']))) for output in rules.trimmed_short_reads_qc.output] for s in e]
    input_files['short_read_alignments_qc'] = expand(rules.short_read_alignments_qc.output, asm=wildcards.asm, paired_short_readset_id=range(0, len(config['assemblies'][wildcards.asm]['paired_short_read_paths'])))
    input_files['second_polish_qc'] = str(rules.second_polish_qc.output).format(asm=wildcards.asm)
  if 'reference' in config['assemblies'][wildcards.asm]:
    # reference provided, perform scaffolding
    input_files['scaffolding_qc'] = str(rules.scaffolding_qc.output).format(asm=wildcards.asm)
  else:
    # no reference provided, don't scaffold
    input_files['contigs'] = str(rules.rename_contigs.output.contigs).format(asm=wildcards.asm)
  input_files['assembly_comparison_report'] = str(rules.assembly_comparison.output.report_tsv).format(asm=wildcards.asm)
  return input_files

rule final_assembly_report:
  input:
    unpack(get_final_assembly_report_input)
  output:
    touch("{asm}/reports/final.touch")
