auto-asm generated a separate folder for each assembly in the specified output
directory, named after the assembly id from the configuration file (e.g. asm1).
Paths specified in the report are children of each assembly directory, unless
specified otherwise.


Raw Data
--------

The raw long and short read data given in the configuration has been symlinked
to ``raw_data``. Statistics such as output, read N50/N90, coverage, length
distributions, etc. can be found in ``reports/raw_data``.


Assemblies
----------

Depending on the provided input data, auto-asm generated at least two, but up to
four, assemblies:

1. Initial *de novo* (always)

    - ``assembly/initial``

2. Long read polished (always)

    - ``assembly/first_polish``

3. Short read polished (if short reads provided)

    - ``assembly/second_polish``

4. Reference-guided Scaffolded (if reference provided)

    - ``assembly/scaffolded``

If short reads were provided, the second polish contig ids were standarized to
"CTG\_XXX" and a mapping file between old and new ids was generated, along with
the new contig FASTA file, in ``assembly/``. Otherwise, the first polish contig
ids were standarized instead.

The assembly comparison report can be found in ``reports/assembly/comparison``.
