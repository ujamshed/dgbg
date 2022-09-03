# Docking Grid Box Generator (DGBG)

Python module to automatically generate a grid box around a ligand present in a pdb file.

![GitHub last commit](https://img.shields.io/github/last-commit/ujamshed/dbg)
<!-- last commit -->

![GitHub issues](https://img.shields.io/github/issues-raw/ujamshed/dbg)
<!-- tracks issues in  project and gets updated automatically -->

![GitHub pull requests](https://img.shields.io/github/issues-pr/ujamshed/dbg)
<!--: that tracks pull requests.-->

![GitHub](https://img.shields.io/github/license/ujamshed/dbg)
<!-- license -->

## Preview

![screenshot](result.gif)

## Table of contents

- [Demo-Preview](#demo-preview)
- [Table of contents](#table-of-contents)
- [Installation](#installation)
- [Usage](#usage)
- [Development](#development)
- [Acknowledgements](#acknowledgements)
- [License](#license)


## Installation
[(Back to top)](#table-of-contents)

To use this project, pip install the module on your device using the command below:

```pip install dbg```

## Usage
[(Back to top)](#table-of-contents)

```python
from dbg import dbg

box = dbg.binding_box('ex_file_path.pdb')

ligand_data = box.get_ligand_data('ligand_ID', 'chain_ID')

centroid, box_size = box.create_box(spacing=0.375, padding=np.array([0, 0, 0]))

# Visualize using the following commands in Jupyter Notebook.

result = box.show_result(centroid, box_size)

result

```

## Examples
[(Back to top)](#table-of-contents)

1. General example use case.

2. Adjusting grid box after generation.

3. Integration with DeepChem and its docking functionality.

4. Integration with AutoDock Vina.

5. Multiple identical ligands on the same protein chain.

## Development
[(Back to top)](#table-of-contents)

1. Git Clone Repository
```git clone git@github.com:ujamshed/dbg.git```

2. Create Virtual Environment
```python3  -m venv .venv```

3. Activate Virtual Environment
```source .venv/bin/activate```

4. Install Project Dependencies
```pip install -r requirements.txt```

## Acknowledgements

This software uses the following open source packages:

- [NumPy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [MDTraj](https://www.mdtraj.org/1.9.8.dev0/index.html)
- [nglview](https://github.com/nglviewer/nglview)

## License
[(Back to top)](#table-of-contents)

[MIT](https://opensource.org/licenses/MIT)