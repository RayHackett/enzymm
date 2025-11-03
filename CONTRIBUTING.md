# Contributing to EnzyMM

For bug fixes or new features, please file an issue before submitting a
pull request. If the change isn't trivial, it may be best to wait for
feedback.

## Coding guidelines

### Versions

This project targets Python 3.7 or later.

Python objects should be typed; Since this library also targets older python versions,
please use the typing module as you would prior to python 3.9.

### Format

Code is formated with [Black](https://github.com/psf/black)

## Setting up a local repository

```console
$ git clone --recursive https://github.com/rayhackett/enzymm
```

## Running tests

Tests are written as usual Python unit tests with the `unittest` module of
the standard library. Running them requires the extension to be built
locally:

```console
$ python -m pip install -v -e .
$ python -m unittest enzymm.tests -vv
```