{% set asm = snakemake.wildcards['asm'] %}
{% set long_readset_id = snakemake.wildcards['long_readset_id']|int %}
Read length distribution for readset
{{ snakemake.wildcards['long_readset_id'] }}
({{ snakemake.config['assemblies'][asm]['long_read_paths'][long_readset_id] }})
of
{% if 'name' in snakemake.config['assemblies'][asm] %}
{{ snakemake.config['assemblies'][asm]['name'] }}
{% else %}
{{ asm }}
{% endif %}
