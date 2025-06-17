# orthomcl
Python module to find VEuPathDB species orthologues via the OrthoMCL database.

## Installation

The easiest way to install the `orthomcl` python module is via pip:
```bash
pip install git+https://github.com/ulido/orthomcl
```

## API

Querying the database uses a simple key-value interface. This illustrative
example queries OrthoMCL for the Trypanosoma brucei orthologs of the
Leishmania mexicana PF16 gene (ID LmxM.20.1400):
```python-repl
>>> from orthomcl import OrthoMCL
>>> gene_id = "LmxM.20.1400"
>>> results = OrthoMCL[gene_id]
>>> results["tbre"]
[OrthoEntry(group='OG7_0009222', organism='tbre', gene_id='Tb427_010020200')]
```

The database can also be queried by OrthoMCL group name:
```python-repl
>>> results = OrthoMCL.by_group("OG7_0009222")
>>> results["lmex"]
[OrthoEntry(group='OG7_0009222', organism='lmex', gene_id='LmxM.20.1400')]
```

## Command line interface

The OrthoMCL database can be queried via the command line using
```bash
orthomcl -o tbre LmxM.20.1400
```

This looks for orthologues of the /Leishmania mexicana/ gene `LmxM.20.1400` in the /Trypanosoma brucei/ genome. The `-o` flag is optional, if not given, all ortholog gene IDs in all organisms are listed. More than one query gene can be given.