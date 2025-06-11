MOFChecker is an algorithm designed for MOF structure curation. It includes duplicated structure screening, geometric error checking, and charge error checking in MOF structures. Based on the check result, we also provide a workflow to clean and heal the structures with geometric or charge errors. The previous version can be found in (https://github.com/kjappelbaum/mofchecker) and (https://github.com/Au-4/mofchecker_2.0)

1.Installation

Development version:

```bash
pip install git+https://github.com/sxm13/mofchecker_2.0
```

Note that you need to install [zeopp](https://anaconda.org/conda-forge/zeopp-lsmo) if you want to use the porosity features.

```bash
conda install -c conda-forge zeopp-lsmo
```

2.Manual
### Input file
We support standard CIF format file. To avoid potential issues, we strongly recommend you first read and rewrite your CIF file with ase.
```python
from ase.io import read, write
mof = read(<path_to_cif>)
write(<path_to_new_cif>, mof)
```

### Command line interface

```bash
mofchecker --help # list options
mofchecker structure1.cif structure2.cif  # prints JSON output
mofchecker -d has_metal -d has_atomic_overlaps *.cif  # compute only selected descriptors
```

### In Python

```python
from mofchecker import MOFChecker
mofchecker = MOFChecker.from_cif(<path_to_cif>)
# or: MOFChecker(structure=my_pymatgen_structure)

# Test for OMS
mofchecker.has_oms

# Test for clashing atoms
mofchecker.has_atomic_overlaps

# Run basic checks on a list of cif paths (sample_structures)
results = []

for structure in sample_structures:
    mofchecker = MOFChecker.from_cif(structure)
    results.append(mofchecker.get_mof_descriptors())
```
Checking and cleaning cases can be found in tests folder.
The positive charge of metal site can be obtained from Oximachine (https://github.com/kjappelbaum/oximachinerunner)


### ⚖️ License

The code in this package is licensed under the Academic Non-Commercial Share-Alike(ANCSA) License.

