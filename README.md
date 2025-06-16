# orthomcl
Python module to find VEuPathDB species orthologues via the OrthoMCL database.

## Command line interface

The OrthoMCL database can be queried via the command line using
```bash
orthomcl -o tbre LmxM.20.1400
```

This looks for orthologues of the /Leishmania mexicana/ gene `LmxM.20.1400` in the /Trypanosoma brucei/ genome. The `-o` flag is optional, if not given, all ortholog gene IDs in all organisms are listed. More than one query gene can be given.