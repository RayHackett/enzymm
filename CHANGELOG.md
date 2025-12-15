# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [v0.3.1] - 2025-12-15
[v0.3.1]: https://github.com/rayhackett/enzymm/compare/v0.3.0..v0.3.1

### Fixed
- Added check for molecule id - molecule id will be derived from structure file or failing that the filename
- Jess object will only be created once per template size batch and passed along
- Added a cif file as input to the unittests

## [v0.3.0] - 2025-11-25
[v0.3.0]: https://github.com/rayhackett/enzymm/compare/v0.2.0..v0.3.0

### Added
- Added [Documentation](https://enzymm.readthedocs.io/en/latest/)
- `--include-template` flag added. Set it to include the template pdb in the output
- templates and molecules can now be dumped too (with the ability to transform them) -> bumped requirement to PyJess 0.9.0

### Changed
- **breaking** different structures (molecules, hits, templates) will now be seperated in `MODEL`/`ENDMDL` blocks (compatibility with some molecular viewers).
- **breaking** `--transform` behaviour changed. Setting it will write in the template coordinate system.
- `--transform` and `--include-query` are now compatible. `--transform` will write one output pdb per matched template

### Fixed
- removed two duplicate 3-residue templates; we are down to 7605 templates in total now
- fixed some EC numbers which contained apostrophes
- will now properly recognize CIF structures starting with a comment line as CIF structures
- outsourced dumping molecules, templates and hits to PyJess
- `HEADER` lines in the pdb output were changed to regular `REMARK` lines.
- improved various docstrings
- testing for pdb file output

## [v0.2.0] - 2025-09-17
[v0.2.0]: https://github.com/rayhackett/enzymm/compare/v0.1.7..v0.2.0

### Added
- Support for mmCIF files through [`gemmi`](https://gemmi.readthedocs.io/). integration in PyJess version `>= 0.7.0`
- Support for [Selenocysteine](https://en.wikipedia.org/wiki/Selenocysteine) and [Pyrrolysine](https://en.wikipedia.org/wiki/Pyrrolysine) residues in templates
- Added reference_residue and reference_pdb attributes to AnnotatedResidue Class
- Added iter and len methods to Residue and Template Classes in line with PyJess

### Changed
- **breaking**: Command line argument `-j`/`--jess` replaced with `-p`/`--parameters`
- **breaking**: Command line argument `-n`/`--n-jobs` replaced with `-j`/`--jobs`
- **breaking**: Residue.residue_name attribute replaced with Residue.name
- **breaking**: Residue.residue_number attribute replaced with Residue.number
- **breaking**: Missing input files will now raise errors instead of just a warning
- **breaking**: Bugfix related to max_candidates in PyJess where supplying many templates might supress some matches (requires PyJess version `>= 0.7.0`)
- Will now match sites across chain interfaces even if not specified in the template (requires PyJess version `>= 0.7.0`)

### Fixed
- Should run a lot faster due to many [optimizations](https://pyjess.readthedocs.io/en/latest/guide/optimizations.html) to PyJess version `>= 0.6.0`
- Improved handling of logistic regression models and model ensembles through new classes
- Included PyJess version in the output tsv file
- Included Command line argument (if run from cli) in the output tsv file
- enzymm.template.load_templates() will now load supplied templates by default without further kwargs
- enzymm.template.load_templates() will now by default use only one thread (improves performance)
- Removed useless argument conservation_score of Matcher. Behaviour unchanged.

## [v0.1.7] - 2025-08-21
[v0.1.7]: https://github.com/rayhackett/enzymm/compare/v0.1.6..v0.1.7

### Fixed
- Improved thread-based parallelism over both molecules and batches of templates to avoid waiting for costly molecule/template searches
- Fixed the progress bar
- Pickling of template objects

### Added
- bash entry points for apptainer and docker containers

## [v0.1.6] - 2025-08-21
[v0.1.6]: https://github.com/rayhackett/enzymm/compare/v0.1.5..v0.1.6

### Fixed
- Fixed overwriting the path variable. now appending to it.

## [v0.1.5] - 2025-08-21
[v0.1.5]: https://github.com/rayhackett/enzymm/compare/v0.1.4..v0.1.5

### Fixed
- Installing apptainer manually in apptainer workflow to avoid later issues with apt

## [v0.1.4] - 2025-08-21
[v0.1.4]: https://github.com/rayhackett/enzymm/compare/v0.1.3..v0.1.4

### Fixed
- Error in downloading and unpacking oras

## [v0.1.3] - 2025-08-21
[v0.1.3]: https://github.com/rayhackett/enzymm/compare/v0.1.2..v0.1.3

### Fixed
- Write permissions given to apptainer workflow. Skip pypi upload (silently!) if version tag exists.

## [v0.1.2] - 2025-08-20
[v0.1.2]: https://github.com/rayhackett/enzymm/compare/1094cde..v0.1.2

### Fixed
- Github actions workflow for apptainer should now generate a release artifact and upload via oras.

## [v0.1.1] - 2025-08-20
[v0.1.1]: https://github.com/rayhackett/enzymm/compare/26de8cc..1094cde

### Fixed
- Github actions workflow for apptainer runs directly on unbuntu base without --fakeroot. New tag should satisfy pypi.

## [v0.1.0] - 2025-08-20
[v0.1.0]: https://github.com/rayhackett/enzymm/compare/603a1bd..26de8cc

### Added
- Unittests for Annotated Templates and residues
- Added Apptainer via ORAS built in github actions

### Changed
- Disabling checks on multichain template and query pairs for chain relationships !

### Fixed
- Eliminated the wait times until all molecules had been scanned before smaller templates used for searches
- Removed unnecessary duplicate EMO function tags from residue annotations
- Fixed some inconsistent ptm residue annotations
- Eliminated some unnecessary steps in the unittests which added considerable compute time
- Check in the CLI for pairwise distances for which no prediction models exist if the tag --unfilteredqwas not passed
- Spelling and badges in the README.md

## [v0.0.3] - 2025-08-07
[v0.0.3]: https://github.com/rayhackett/enzymm/compare/6dad6cd..603a1bd

### Added
- Added Information to README.md
- Added Docker Container via Github actions to GHCR

### Changed
- Disabling checks on multichain template and query pairs for chain relationships

### Fixed
- Fixed Github actions url
- cpu counting on linux systems
- Fixed attempt at filtering for pairwise distances without determined logistic models

## [v0.0.2] - 2025-07-20
[v0.0.2]: https://github.com/rayhackett/enzymm/compare/ea71726..6dad6cd

### Fixed
- Fixed Github actions to properly build from source

## [v0.0.1] - 2025-07-20
[v0.0.1]: https://github.com/RayHackett/enzymm/tree/ea7172665215e5073f70b27ce2aa07a49b72eb48

Initial release.