rule first_polish:
  """ correct sequence errors using long read mappings to the de novo assembly """
  input:
    alignments = rules.sort_long_read_alignments_polish.output.sorted_alignments,
    alignment_index = rules.index_long_read_alignments_polish.output.index,
    reference = rules.initial_assembly.output.contigs_fasta,
    reference_index = rules.index_initial_assembly.output.index
  output:
    variants = protected("{asm}/assembly/first_polish/variants.gff"),
    consensus_fasta = protected("{asm}/assembly/first_polish/consensus.fasta"),
    consensus_fastq = protected("{asm}/assembly/first_polish/consensus.fastq")
  benchmark:
    "{asm}/benchmarks/first_polish.txt"
  log:
    stdout = "{asm}/logs/first_polish.out",
    stderr = "{asm}/logs/first_polish.err"
  shell:
    str(Path(config['paths']['smrtcmds_bin']) / 'arrow') +
    " -j{threads} " +
    "{input.alignments} " +
    "-r {input.reference} " +
    "-o {output.variants} " +
    "-o {output.consensus_fasta} " +
    "-o {output.consensus_fastq} " +
    "2> {log.stderr} > {log.stdout}"
