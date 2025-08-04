
# EnzyMM - The Enzyme Motif Miner

## ️Overview

Enzyme Motif Miner uses geometric template matching to identify known arrangements of catalytic residues called templates in protein structures. It searches protein structures provided by the user against a database of templates. `EnzyMM` ships with a library of catalytic templates derived from the [Mechanism and Catalytic Site Atlas](https://www.ebi.ac.uk/thornton-srv/m-csa/) (M-CSA) but you can also generate your own. These templates represent consensus arrangements of catalytic sites found in active sites of experimental protein structures.   

As catalytic sites are both highly conserved and absolutely critical for the function of a protein, identifying them offers many biological insights. This method has two key advantages. Firstly, as it doesn't rely on sequence or (global) fold similarity, similar catalytic arrangements can be found accross great evolutionary distances offering insights into the divergence or even convergence of enyzmes. Secondly, as geometric matching is very fast, `EnzyMM` scales along side databases of predicted protein structures. Expect to scan a protein structure in a matter of seconds on consumer laptops.  

As a database driven method, `EnzyMM` is inherently limited by the coverage of residue arrangements in its template library. The provided template library covers nearly the entire M-CSA and thus around 3/4 of enzyme mechanisms classified by the Enzyme Commission to the 3rd level. Catalytic arrangements not found in the PDBe won't be included in the M-CSA. Of course, the user can also provide their own library of templates. While primarily intended for catalytic sites, you are invited to search with your own library of templates.  

For the actual geometric matching `EnzyMM` relies on [PyJess](https://github.com/althonos/pyjess) - a [Cython](https://cython.org/) wrapper of [Jess](https://github.com/iriziotis/jess).

[![Actions](https://img.shields.io/github/actions/workflow/status/RayHackett/enzymm/test.yml?branch=main&style=flat-square&maxAge=300)](https://github.com/RayHackett/Enzymm/actions/workflows/test.yml)
[![License](https://img.shields.io/badge/license-mit.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/RayHackett/enzymm/)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/rayHackett/enzymm/blob/main/CHANGELOG.md)
[![Issues](https://img.shields.io/github/issues/RayHackett/enzymm.svg?style=flat-square&maxAge=600)](https://github.com/RayHackett/enzymm/issues)
[![PyPI](https://img.shields.io/pypi/v/enzymm.svg?style=flat-square&maxAge=3600)](https://pypi.python.org/pypi/gecco-tool)
[![Versions](https://img.shields.io/pypi/pyversions/enzymm.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/enzymm/#files)
[![Wheel](https://img.shields.io/pypi/wheel/enzymm?style=flat-square&maxAge=3600)](https://pypi.org/project/enzymm/#files)


## 🔧 Installing EnzyMM

`EnzyMM` is implemented in [Python](https://www.python.org/), 
and supports [all versions](https://endoflife.date/python) from Python 3.8 on Linux and MacOS. It requires
additional libraries that can be installed directly from
[PyPI](https://pypi.org), the Python Package Index.

Use [`pip`](https://pip.pypa.io/en/stable/) to install `EnzyMM` on your
machine:
```console
$ pip install enzymm
```

This will both install `EnzyMM` and also download a library of catalytic templates together with important metadata. This requires around 16MB of data to be downloaded.
It should also run on windows (though this is not tested for on release).

## 🔎 Running EnzyMM

Once `EnzyMM` is installed, you can run it from the terminal. The user can either provide a path to a single protein structure `-i` or to run multiple queries at once, the path to a text file `-l` which itself contains a list of paths to protein structures.
Optionally, an output directory for pdb structures of the identified matches per query protein can be supplied with the `--pdbs` flag.

```console
$ enzymm -i some_structure.pdb -o results.tsv --pdbs dir_to_save_matches
```

Additional parameters of interest are:

- `--n-jobs` or `-n`, which controls the number of threads used to parallelize the search.
  By default, it will use one thread less than available on your system using
  [`os.cpu_count`](https://docs.python.org/3/library/os.html#os.cpu_count).
- `--unfiltered` or `-u`, which disables filtering of matches by RMSD and residue orientation.
  By default, filtering is enabled.
- `--skip-smaller-hits`, which skips searches with smaller templates on a query
  if a match to a larger template has already been found.
- `--jess` or `-j`, which controls the RMSD threshold and pairwise distance threshold applied. By default sensible thresholds are selected. Refer to the Docs for details
- `--template-dir` or `-t`, though which the user may supply their own template library. By default, a library of catalytic templates derived from the M-CSA is loaded.
- `--conservation-cutoff` or `-c`, which can be set to exclude atoms with B-factors or pLDDT scores below this threshold from matching. This is not set by default.

Further, `EnyzMM` is designed with modularity in mind and comes with a fully usable internal API.
Please refer to the Docs for further reference.

## 🖹 Results

`EnzyMM` will create a single output file:

- `{output}.tsv`: A `.tsv` file containing a summary of all results. One row is printed per match.

For visual exploration of matches, you can optionally save an alignment of the template and the matched query residues to a pdb file which can be viewed with any pdb viewer.
To do so, supply an output directory after the `--pdbs` flag for the `.pdb` files.

This will also create:

- `{pdbs_dir}/{query_identifier}_matches.pdb`: One `.pdb` file per query with a structural alignment between template and query residues. This can be further configured.

Add additional information to each `.pdb` file with the following flags:

- `--transform`, which causes the query to be aligned to the the template instead of vice versa.
- `--include-query`, which also writes the entire query pdb structure to the `.pdb` file

Currently, `--transform` and `--include-query` should not be used together.
Hopefully I'll get around to fixing this soon.


## 🔖 Citations
`EnyzMM` is academic software but relies on many previous approaches.  
  
`EnzyMM` itself can not yet be cited but a preprint is in preparation.
We intend to publish during the summer of 2025.  

We kindly ask you to cite both:  
- PyJess for instance as:
> PyJess, a Python library binding to Jess (Barker *et al.*, 2003).
- Mechanism and Catalytic Site Atlas as:
> Ribeiro AJM et al. (2017), Nucleic Acids Res, 46, D618-D623. Mechanism and Catalytic Site Atlas (M-CSA): a database of enzyme reaction mechanisms and active sites. DOI:10.1093/nar/gkx1012. PMID:29106569.

## 💭 Feedback

### ⚠️ Issue Tracker

Please report any bugs or feature requests though the [GitHub issue tracker](https://github.com/RayHackett/enzymm/issues).
Please also feel free to ask any questions and I will do my best to answer them.  
If reporting a bug, please include as much information as you can about the issue and try to recreate the same bug.
Ideally include a little test example so I can quickly troubleshoot.

### 🏗️ Contributing
Contributions are more than welcome!
Raise an issue or shoot me an email under `r.e.hackett` AT `lumc.nl`  
I'm happy to help.

## ⚖️ License

This software is provided under an [MIT](https://choosealicense.com/licenses/mit/) licence. Though conceived at the [EMBL-EBI](https://www.ebi.ac.uk/) in Hinxton, UK in the [Thornton Group](https://www.ebi.ac.uk/research/thornton/), `EnzyMM` is now developed by Raymund Hackett and the [Zeller Group](https://zellerlab.org/) at the [Leiden University Medical Center](https://www.lumc.nl/en/) in Leiden in the Netherlands with continuing support from the Thornton Group.