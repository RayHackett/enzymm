# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).
## Unreleased
## [v0.2.0-alpha.1] - 2025-09-02
[v0.2.0-alpha.1]: https://github.com/rayhackett/enzymm/compare/v0.1.7..v0.2.0-alpha.1

### Changed
- Command line argument `-j`/`--jess` replaced with `-p`/`--parameters`
- Command line argument `-n`/`--n-jobs` replaced with `-j`/`--jobs`
- Missing input files will now raise errors instead of just a warning
- Bugfix where supplying many templates might supress some matches (requires PyJess version `>= 0.7.0a1`)
- Will now match sites across chain interfaces even if not specified in the template (requires PyJess version `>= 0.7.0a2`)

### Fixed

- Should run a lot faster due to many [optimizations](https://pyjess.readthedocs.io/en/latest/guide/optimizations.html) to PyJess version `>= 0.6.0`
- Improved handling of logistic regression models and model ensembles through new classes
- Included PyJess version in the output tsv file
- Included Command line argument (if run from cli) in the output tsv file
- enzymm.template.load_templates() will now load supplied templates by default without further kwargs
- enzymm.template.load_templates() will now by default use only one thread (improves performance)

## [v0.1.7] - 2025-08-21
[v0.1.6]: https://github.com/rayhackett/enzymm/compare/v0.1.6..v0.1.7

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