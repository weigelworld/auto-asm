rule second_polish:
  """ correct sequence errors using short read mappings to the polished assembly """
  input:
    short_read_alignments = rules.sort_short_read_alignments_polish.output.sorted_alignments,
    short_read_alignments_index = rules.index_short_read_alignments_polish.output.index,
    assembly = rules.first_polish.output.consensus_fasta
  conda:
    '../envs/second_polish.yaml'
  params:
    out_dir = lambda wildcards, output: str(Path(output['consensus']).parent),
    java_mem_mb = round(config['rules']['second_polish']['resources']['mem_mb'] * 0.85)  # Leave 15% memory as buffer for Java
  output:
    consensus = protected("{asm}/assembly/second_polish/{asm}.fasta"),
    changes = protected("{asm}/assembly/second_polish/{asm}.changes")
  benchmark:
    "{asm}/benchmarks/second_polish.txt"
  log:
    stdout = "{asm}/logs/second_polish.out",
    stderr = "{asm}/logs/second_polish.err"
  shell:
      "pilon " +
      "-Xmx{params.java_mem_mb}M " +
      "--threads {threads} " +
      "--genome {input.assembly} " +
      "--bam {input.short_read_alignments} " +
      "--outdir {params.out_dir} " +
      "--output {wildcards.asm} " +
      "--changes " +
      "2> {log.stderr} > {log.stdout}"
