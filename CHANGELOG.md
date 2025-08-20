# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [v0.1.0] - 2025-08-20
[v0.1.0]: https://github.com/rayhackett/enzymm/compare/v0.1.0..0603a1bd

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