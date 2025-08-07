from __future__ import annotations

import collections
import math
import warnings
import csv
import itertools
import io
import sys
import os
import json
import functools
from multiprocessing.pool import ThreadPool
from typing import List, Tuple, Dict, Optional, IO, ClassVar, Iterable
from functools import cached_property
from dataclasses import dataclass, field
from pathlib import Path

import rich
import pyjess
from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

from enzymm.template import Template, AnnotatedTemplate, Vec3, check_template
from enzymm.utils import chunks, ranked_argsort, DummyPool
from enzymm.utils import PROTEINOGENIC_AMINO_ACIDS, SPECIAL_AMINO_ACIDS

__all__ = [
    "LogisticRegressionModel",
    "Match",
    "Matcher",
]


@dataclass(frozen=True)
class LogisticRegressionModel:
    """
    Class for storing a Logistc Regression Model for predicting if a match is correct.

    f(Xn) = 1/(1+e^-(beta0 + beta1*x1 + beta2*x2 + ...))

    Trained on data for matches with a particular size at a particular pairwise distance

    Attributes:
        coeficients: `list` of `floats` beta coefficients
        intercept: `float` beta0 intercept
        threshold: `float` optimal threshold for this model
    """

    coefficents: List[float]
    intercept: float
    threshold: float


@dataclass
class Match:
    """
    Class for storing annotated PyJess hits.

    This class is a wrapper around `pyjess.Hit` and the original template object
    that was used for the query.

    Attributes:
        hit: `~pyjess.Hit` instance
        complete: `bool` If the query matched all other templates within the same cluster
        pairwise_distance: `float` Pairwise distance at which this match was found
        index: `int` internal index of this match. Default 0

    """

    hit: pyjess.Hit
    complete: bool = field(default=False)
    pairwise_distance: float = field(default=0)
    index: int = field(default=0)
    _logistic_regression_models: ClassVar[
        Dict[str, Dict[str, List[LogisticRegressionModel]]]
    ]

    def dumps(self, header: bool = False) -> str:
        """
        Dump `Match` to a string. Calls `Match.dump()`

        Arguments:
            header: if a header line should be dumped to the string too.
        """
        buffer = io.StringIO()
        self.dump(buffer, header=header)
        return (
            buffer.getvalue()
        )  # returns entire content temporary file object as a string

    def dump2pdb(
        self, file: IO[str], include_query: bool = False, transform: bool = False
    ):
        """
        Dump the 3D coordinates of the `Match` to a '.pdb' file.

        Arguments:
            file: file type object to write to
            include_query: `bool` If the full query molecule should be written too.
            transform: `bool` If the matched atoms should be written to the template reference frame.

        Note:
            The query is never transformed to the target reference frame.

        Note:
            By default, atoms are written in the coordinate reference frame of the query.
        """

        # TODO option to include template atoms too. esp. with --transform
        # this requires a dump method for pyjess.TemplateAtoms
        # if NOT transform, rotate the template atoms into the query reference frame.
        # if transform then in pymol once could fetch the template reference pdb
        # select the template atoms in the reference_pdb (selection called ref_site)
        # then call align ref_site, csa*** name of template
        # TODO once writing one match to a single file is done
        # with --transform and --include query, rotate the query to the template reference too
        # simply by using the roation matrix or by selecting the atoms in the query
        # and aligning to the matched residues object.

        def write_atom_line(atom: pyjess.Atom) -> str:
            one_char_elements = {
                "H",
                "B",
                "C",
                "N",
                "O",
                "F",
                "P",
                "S",
                "K",
                "V",
                "Y",
                "I",
                "W",
                "U",
            }
            if atom.element in one_char_elements:
                return f"ATOM  {atom.serial:>5}  {atom.name:<3s}{atom.altloc if atom.altloc is not None else '':<1}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:>4}{atom.insertion_code:1s}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.temperature_factor:>6.2f}      {atom.segment:<4s}{atom.element:>2s} \n"
            else:
                return f"ATOM  {atom.serial:>5} {atom.name:<4s}{atom.altloc if atom.altloc is not None else '':<1}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:>4}{atom.insertion_code:1s}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.temperature_factor:>6.2f}      {atom.segment:<4s}{atom.element:>2s} \n"

        if include_query:  # write the original query molecule too
            file.write(f"HEADER MOLECULE_ID {self.hit.molecule().id}\n")
            for atom in self.hit.molecule(transform=transform):
                file.write(write_atom_line(atom))
            file.write("END\n\n")

        file.write(
            f"HEADER {self.predicted_correct} MATCH {self.hit.molecule().id} {self.index}\n"
        )

        file.write(
            f'REMARK TEMPLATE_PDB {str(self.hit.template.pdb_id)}_{",".join(set(res.chain_id for res in self.hit.template.residues))}\n'
        )

        # alias for improved readability
        template = self.hit.template
        cluster = template.cluster

        if cluster:
            file.write(
                f"REMARK TEMPLATE CLUSTER {cluster.id}_{str(cluster.member)}_{str(cluster.size)}\n"
            )
        if template.represented_sites:
            file.write(f"REMARK TEMPLATE RESIDUES {template.template_id_string}\n")
        file.write(f"REMARK MOLECULE_ID {str(self.hit.molecule().id)}\n")
        file.write(f"REMARK MATCH INDEX {self.index}\n")

        if transform:
            file.write("REMARK TEMPLATE COORDINATE FRAME\n")
        else:
            file.write("REMARK QUERY COORDINATE FRAME\n")

        for atom in self.hit.atoms(
            transform=transform
        ):  # if transform == True then coordinates are transformed to the template reference frame
            file.write(write_atom_line(atom))
        file.write("END\n\n")

    def dump(
        self, file: IO[str], header: bool = False, predict_correctness: bool = True
    ):
        """
        Dump the information associated with a `Match` to a '.tsv' like line.

        Arguments:
            file: `file-like` object to write to
            header: `bool` If a header line should be written too
            predict_correctness: 'bool' Wether to predict if the match is correct or not

        Note:
            Coordinate information is not written.
        """
        writer = csv.writer(
            file, dialect="excel-tab", delimiter="\t", lineterminator="\n"
        )
        # aliases for improved readability
        template = self.hit.template
        cluster = self.hit.template.cluster

        if header:
            writer.writerow(
                [
                    "query_id",
                    "pairwise_distance",
                    "match_index",
                    "template_pdb_id",
                    "template_pdb_chains",
                    "template_cluster_id",
                    "template_cluster_member",
                    "template_cluster_size",
                    "template_effective_size",
                    "template_dimension",
                    "template_mcsa_id",
                    "template_uniprot_id",
                    "template_ec",
                    "template_cath",
                    "template_multimeric",
                    "query_multimeric",
                    "query_atom_count",
                    "query_residue_count",
                    "rmsd",
                    "log_evalue",
                    "orientation",
                    "preserved_order",
                    "completeness",
                    "predicted_correct",
                    "matched_residues",
                    "number_of_mutated_residues",
                    "number_of_side_chain_residues_(template,reference)",
                    "number_of_metal_ligands_(template,reference)",
                    "number_of_ptm_residues_(template, reference)",
                    "total_reference_residues",
                ]
            )

        content = [
            str(self.hit.molecule().id),
            str(self.pairwise_distance),
            str(self.index),
            str(template.pdb_id if template.pdb_id else ""),
            (",".join(set(res.chain_id for res in template.residues))),
            str(cluster.id if cluster else ""),
            str(self.hit.template.cluster.member if self.hit.template.cluster else ""),
            str(cluster.size if cluster else ""),
            str(template.effective_size),
            str(template.dimension),
            str(template.mcsa_id if template.mcsa_id else ""),
            str(template.uniprot_id if template.uniprot_id else ""),
            ",".join(template.ec if template.ec is not None else ""),
            ",".join(template.cath if template.cath else ""),
            str(template.multimeric),
            str(self.multimeric),
            str(self.query_atom_count),
            str(self.query_residue_count),
            str(round(self.hit.rmsd, 5)),
            str(round(self.hit.log_evalue, 5)),
            str(round(self.orientation, 5)),
            str(self.preserved_resid_order),
            str(self.complete),
            str(self.predicted_correct) if predict_correctness else "",
            (",".join("_".join(t) for t in self.matched_residues)),
        ]

        # check if the template was annotated with M-CSA information
        if isinstance(template, AnnotatedTemplate):
            content.extend(
                [
                    str(template.number_of_mutated_residues),
                    ",".join(str(i) for i in template.number_of_side_chain_residues),
                    ",".join(str(i) for i in template.number_of_metal_ligands),
                    ",".join(str(i) for i in template.number_of_ptm_residues),
                    str(template.total_reference_residues),
                ]
            )
        else:
            content.extend(["", "", "", "", "", ""])

        writer.writerow(content)

    def get_identifying_attributes(self) -> Tuple[int, int, int]:
        """
        `tuple` of (M-CSA id, cluster id and template dimension).`
        """
        # return the tuple (hit.template.m-csa, hit.template.cluster.id, hit.template.dimension)
        template = self.hit.template
        return (
            template.mcsa_id,
            template.cluster.id,
            template.dimension,
        )

    @property
    def predicted_correct(self) -> bool:
        """`bool`: If the match is predicted as correct based on logistic regression models."""
        if (
            str(self.hit.template.effective_size)
            not in self._logistic_regression_models
        ):  # Since no models for 5+ residue matches is provided, these are predicted as true
            return True
        else:
            try:
                predictions = []
                for model in self._logistic_regression_models[
                    str(self.hit.template.effective_size)
                ][str(self.pairwise_distance)]:
                    value = 1 / (
                        1
                        + math.e
                        ** -(
                            model.intercept
                            + model.coefficents[0] * self.hit.rmsd
                            + model.coefficents[1] * self.orientation
                        )
                    )
                    predictions.append(value >= model.threshold)

                # majority decision from all models
                return bool(
                    sum(predictions)
                    >= round(
                        (
                            len(
                                self._logistic_regression_models[
                                    str(self.hit.template.effective_size)
                                ][str(self.pairwise_distance)]
                            )
                            / 2
                        ),
                        0,
                    )
                )

            except KeyError as exc:
                raise KeyError(
                    f"Missing appropriate model parameters to predict correctness. Encountered either unexpected dictionary structure or no models for the pairwise distance {self.pairwise_distance} were provided"
                ) from exc
            except IndexError as exc:
                raise IndexError(
                    "Missing coefficients for both RMSD and Residue Orientation. Expecting models with 2 coeficients."
                ) from exc

    @cached_property
    def atom_triplets(self) -> List[Tuple[pyjess.Atom, pyjess.Atom, pyjess.Atom]]:
        """
        `list`: of `~pyjess.Atom` triplets belonging to the same matched query residue.
        """
        # list with matched residues
        # # Hit.atoms is a list of matched atoms with all info on residue numbers and residue chain ids and atom types, this should conserve order if Hit.atoms is a list!!!
        atom_triplets = []
        for atom_triplet in chunks(
            self.hit.atoms(transform=True), 3
        ):  # yield chunks of 3 atoms each, transform true because for angle calculation atoms need to be in template reference frame
            if len(atom_triplet) != 3:
                raise ValueError(
                    f"Failed to construct residues. Got only {len(atom_triplet)} ATOM lines"
                )
            # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
            unique_residues = {
                (atom.residue_name, atom.chain_id, atom.residue_number)
                for atom in atom_triplet
            }
            if len(unique_residues) != 1:
                raise ValueError(
                    f"Mixed up atom triplets {unique_residues}. The atoms come from different residues!"
                )
            atom_triplets.append(atom_triplet)
        return atom_triplets

    @property
    def matched_residues(self) -> List[Tuple[str, str, str]]:
        """
        `list`:  with information on all matched query residues.
        Elements have are `tuple` (`~pyjess.Atom.residue_name`, `~pyjess.Atom.chain_id`, `~pyjessAtom.residue_number`)
        """
        return [
            (
                atom_triplet[0].residue_name,
                atom_triplet[0].chain_id,
                str(atom_triplet[0].residue_number),
            )
            for atom_triplet in self.atom_triplets
        ]

    @property
    def multimeric(self) -> bool:
        """
        `bool`: If the matched atoms in the query stem from multiple protein chains
        """
        # note that these are pyjess atom objects!
        return not all(
            atom.chain_id == self.hit.atoms()[0].chain_id for atom in self.hit.atoms()
        )

    @property
    def preserved_resid_order(self) -> bool:
        """
        `bool`: If the residues in the template and in the matched query structure have the same relative order.

        Note:
            This is a good filtering parameter but excludes hits on examples of convergent evolution or circular permutations

        Note:
            Will always return `False` if either template or query is multimeric
        """
        if self.hit.template.multimeric or self.multimeric:
            return False
        else:
            # Now extract relative atom order in hit
            return (
                ranked_argsort(
                    [
                        atom_triplet[0].residue_number
                        for atom_triplet in self.atom_triplets
                    ]
                )
                == self.hit.template.relative_order
            )

    @cached_property
    def match_vector_list(cls) -> List[Vec3]:
        """
        `list` of `Vec3`: of orientation vectors for each matched residue in the query
        """
        # !!! atom coordinates must be in template coordinate system!
        vector_list = []
        for residue_index, residue in enumerate(cls.hit.template.residues):
            first_atom_index, second_atom_index = residue.orientation_vector_indices
            if (
                second_atom_index == 9
            ):  # Calculate orientation vector going from middle_atom to mitpoint between side1 and side2
                middle_atom = cls.atom_triplets[residue_index][first_atom_index]
                side1, side2 = [
                    atom
                    for atom in cls.atom_triplets[residue_index]
                    if atom != middle_atom
                ]
                midpoint = (Vec3.from_xyz(side1) + Vec3.from_xyz(side2)) / 2
                vector_list.append(midpoint - Vec3.from_xyz(middle_atom))
            else:
                # Calculate orientation vector going from first_atom to second_atom_index
                first_atom = cls.atom_triplets[residue_index][first_atom_index]
                second_atom = cls.atom_triplets[residue_index][second_atom_index]
                vector_list.append(
                    Vec3.from_xyz(second_atom) - Vec3.from_xyz(first_atom)
                )
        return vector_list

    @property
    def template_vector_list(self) -> List[Vec3]:
        """
        `list` of `Vec3`: of orientation vectors for each residue in the template
        """
        return [res.orientation_vector for res in self.hit.template.residues]

    @property
    def orientation(self) -> float:  # average angle
        """
        `float`: The arithmetic mean of per-residue orientation angles
        for matched pairs of template and query residues
        """
        if len(self.template_vector_list) != len(self.match_vector_list):
            raise ValueError(
                "Vector lists for Template and matching Query structure had different lengths."
            )

        # now calculate the angle between the vector of the template and the query per residue
        angle_list = []
        for i in range(len(self.template_vector_list)):
            angle_list.append(
                self.template_vector_list[i].angle_to(self.match_vector_list[i])
            )
        return sum(angle_list) / len(angle_list)

    @property
    def query_atom_count(self) -> int:
        """
        `int`: The number of atoms in the query molecule
        """
        return len(self.hit.molecule())

    @property
    def query_residue_count(self) -> int:
        """
        `int`: The number of residues in the query molecule
        """
        all_residue_numbers = set()
        for atom in self.hit.molecule():
            if atom.residue_name in PROTEINOGENIC_AMINO_ACIDS + SPECIAL_AMINO_ACIDS:
                all_residue_numbers.add(atom.residue_number)
        return len(all_residue_numbers)


# Load the logistic regression models from the json into a dictionary
with resource_files(__package__).joinpath("data", "logistic_regression_models.json").open() as f:  # type: ignore
    logistic_regression_models: Dict[str, Dict[str, List[LogisticRegressionModel]]] = {}
    model_dict = json.load(f)

    for template_size, pairwise_dict in model_dict["match_size"].items():
        logistic_regression_models[template_size] = {}
        for pairwide_distance, model_list in pairwise_dict["pairwise_distance"].items():
            list_of_log_regression_models = []
            for param_dict in model_list["model_list"]:
                list_of_log_regression_models.append(
                    LogisticRegressionModel(
                        coefficents=param_dict["coef"],
                        intercept=param_dict["intercept"],
                        threshold=param_dict["threshold"],
                    )
                )
            logistic_regression_models[template_size][
                pairwide_distance
            ] = list_of_log_regression_models

    Match._logistic_regression_models = logistic_regression_models


def load_molecules(
    molecule_paths: List[Path], conservation_cutoff: float = 0, warn: bool = False
) -> List[pyjess.Molecule]:
    """Load query molecules"""
    molecules = []
    stem_counter: Dict[str, int] = collections.defaultdict(int)
    for molecule_path in molecule_paths:
        stem = Path(molecule_path).stem
        stem_counter[stem] += 1
        if stem_counter[stem] > 1:
            # In case the same stem occurs multiple times, create a unique ID using the stem and a running number starting from 2
            unique_id = f"{stem}_{stem_counter[stem]}"
        else:
            unique_id = stem

        molecule = pyjess.Molecule.load(
            str(molecule_path), id=unique_id
        )  # by default it will stop at ENDMDL
        if conservation_cutoff:
            molecule.conserved(conservation_cutoff)
            # molecule = molecule.conserved(args.conservation_cutoff)
            # conserved is a method called on a molecule object that returns a filtered molecule
            # atoms with a temperature-factor BELOW the conservation cutoff will be excluded
        if molecule:
            molecules.append(
                molecule
            )  # load a molecule and filter it by conservation_cutoff
        elif warn:
            warnings.warn(f"received an empty molecule from {molecule_path}")

    if not molecules and warn:
        warnings.warn("received no molecules from input")

    return molecules


class Matcher:
    """
    Class from which a query `~pyjess.Molecule` is matched to a `list` of `Template`.
    """

    _DEFAULT_JESS_PARAMS = {
        3: {"rmsd": 2, "distance": 0.9, "max_dynamic_distance": 0.9},
        4: {"rmsd": 2, "distance": 1.7, "max_dynamic_distance": 1.7},
        5: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
        6: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
        7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
        8: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
    }

    def __init__(
        self,
        templates: List[Template],
        jess_params: Optional[Dict[int, Dict[str, float]]] = None,
        conservation_cutoff: int = 0,
        warn: bool = False,
        verbose: bool = False,
        skip_smaller_hits: bool = False,
        match_small_templates: bool = False,
        cpus: int = (
            len(os.sched_getaffinity(0))
            if sys.platform == "linux"
            else os.cpu_count() or 1
        ),
        filter_matches: bool = True,
        console: rich.console.Console | None = None,
    ):
        """
        Initialize a `Matcher` instance

        Arguments:
            templates: `list` of `Template` to match
            jess_params: `dict` Dictionary of PyJess parameters to apply. Will superseed defaults.
            conservation_cutoff: `float` Atoms below this cutoff will not be matched. Default 0.
            warn: `bool` If warnings about issues during matching should be printed. Default `False`
            verbose: `bool` If progress statements on matching should be printed. Default `False`
            skip_smaller_hits: `bool` Continue searching the query against smaller templates, after a match against any larger one was found. Default `False`
            match_small_templates: `bool` If matches for Templates with fewer than 3 side-chain residues should be reported. Default `False`
            cpus: `int` The number of cpus for multithreading. If 0 (default), use all. If <0 leave this number of threads free.
            filter_matches: `bool` If matches should be filtered by wether they are predicted to be correct. Default `True`

        Note:
            Default jess parameters depend on the size of the template::

                _DEFAULT_JESS_PARAMS = {
                    3: {"rmsd": 2, "distance": 0.9, "max_dynamic_distance": 0.9},
                    4: {"rmsd": 2, "distance": 1.7, "max_dynamic_distance": 1.7},
                    5: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                    6: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                    7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                    8: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                }

        """

        self.templates = templates
        self.cpus = cpus
        self.conservation_cutoff = conservation_cutoff
        self.warn = warn
        self.verbose = verbose
        self.skip_smaller_hits = skip_smaller_hits
        self.match_small_templates = match_small_templates
        self.filter_matches = filter_matches
        self.jess_params = (
            self._DEFAULT_JESS_PARAMS if jess_params is None else jess_params
        )
        self.console = rich.console.Console(quiet=True) if console is None else console

        # # optional code to find duplicates
        # def find_duplicates(objects):
        #     d = collections.defaultdict(list)
        #     for obj in objects:
        #         d[obj].append(obj)

        #     for k,v in d.items():
        #         if len(v) > 1:
        #             for dup in v:
        #                 print(dup.__dict__)
        #             break

        unique_templates = set(self.templates)
        if len(unique_templates) < len(self.templates):
            # find_duplicates(self.templates)
            raise ValueError("Duplicate templates were found.")

        if self.cpus <= 0:
            os_cpu_count = (
                len(os.sched_getaffinity(0))
                if sys.platform == "linux"
                else os.cpu_count()
            )
            if os_cpu_count is not None:
                self.cpus = max(1, os_cpu_count + self.cpus)
            else:
                self.cpus = 1

        self.verbose_print(f"PyJess Version: {pyjess.__version__}")
        self.verbose_print(f"Running on {self.cpus} Thread(s)")
        self.verbose_print(f"Warnings are set to {self.warn}")
        self.verbose_print(
            f"Skip_smaller_hits search is set to {self.skip_smaller_hits}"
        )

        if self.conservation_cutoff:
            self.verbose_print(f"Conservation Cutoff set to {self.conservation_cutoff}")

        # check each template and if it passes add it to the dictionary of templates
        self.templates_by_effective_size: Dict[int, List[Template]] = (
            collections.defaultdict(list)
        )  # Dictionary of List of AnnoatedTemplate objects grouped by Template.effective_size as keys
        for template in templates:
            if check_template(
                template, warn=self.warn
            ):  # returns True if the Template passed all checks or if warn is set to False
                self.templates_by_effective_size[template.effective_size].append(
                    template
                )

        if self.verbose:
            template_number_dict: Dict[int, int] = {}
            for size, template_list in self.templates_by_effective_size.items():
                if not self.match_small_templates and size < 3:
                    continue
                template_number_dict[size] = len(template_list)
            print(
                f"Templates by effective size: {collections.OrderedDict(sorted(template_number_dict.items()))}"
            )

        self.template_effective_sizes = list(self.templates_by_effective_size.keys())
        self.template_effective_sizes.sort(
            reverse=True
        )  # get a list of template_sizes in decending order

        if self.warn:
            smaller_sizes = [i for i in self.template_effective_sizes if i < 3]
            if smaller_sizes:
                small_templates = []
                for i in smaller_sizes:
                    small_templates.extend(self.templates_by_effective_size[i])

                if self.match_small_templates:
                    warnings.warn(
                        f"{len(small_templates)} Templates with an effective size smaller than 3 defined sidechain residues were supplied.\nFor small templates Jess parameters for templates of 3 residues will be used."
                    )
                else:
                    warnings.warn(
                        f"{len(small_templates)} Templates with an effective size smaller than 3 defined sidechain residues were supplied.\nThese will be excluded since these templates are too general."
                    )

                self.verbose_print(
                    "The templates with the following ids are too small:"
                )
                self.verbose_print([st.id for st in small_templates])

    def verbose_print(self, *args):
        """
        Make a print statement only in verbose mode
        """
        if self.verbose:
            print(*args)

    def _get_jess_parameters(self, template_size: int) -> Tuple[float, float, float]:
        if template_size < 3:
            parameter_size = 3
        elif template_size > 8:
            parameter_size = 8
        else:
            parameter_size = template_size

        rmsd = self.jess_params[parameter_size]["rmsd"]
        distance = self.jess_params[parameter_size]["distance"]
        max_dynamic_distance = self.jess_params[parameter_size]["max_dynamic_distance"]

        return rmsd, distance, max_dynamic_distance

    @staticmethod
    def _check_completeness(matches: List[Match]) -> List[Match]:
        # only after all templates of a certain size have been scanned could we compute the complete tag
        # This requries cluster and mcsa.id to be set! Otherwise I assume there is no cluster and therefore the match is complete by default!
        groupable_matches = []
        lone_matches = []

        for match in matches:
            if (
                match.hit.template.mcsa_id is not None
                and match.hit.template.cluster is not None
            ):
                groupable_matches.append(match)
            else:
                lone_matches.append(match)

        grouped_matches = [
            list(g)
            for _, g in itertools.groupby(
                sorted(groupable_matches, key=Match.get_identifying_attributes),
                Match.get_identifying_attributes,
            )
        ]

        for cluster_matches in grouped_matches:
            # For each query check if all Templates assigned to the same cluster targeted that structure
            #
            # TODO report statistics on this: This percentage of queries had a complete active site as reported by the complete tag
            # Filter this by template clusters with >1 member of course or report seperately by the number of clustermembers
            # or say like: This template cluster was always complete while this template cluster was only complete X times out of Y Queries matched to one member
            #
            # check if all the cluster members up to and including cluster_size are present in the group,
            indexed_possible_cluster_members = list(range(cluster_matches[0].hit.template.cluster.size))  # type: ignore
            possible_cluster_members = [x + 1 for x in indexed_possible_cluster_members]

            found_cluster_members = [match.hit.template.cluster.member for match in cluster_matches]  # type: ignore
            found_cluster_members.sort()

            if found_cluster_members == possible_cluster_members:
                for match in cluster_matches:
                    match.complete = True

        for match in lone_matches:
            match.complete = True

        return matches

    @staticmethod
    def _run_jess(
        molecule: pyjess.Molecule,
        templates: List[Template],
        rmsd_threshold: float = 2.0,
        distance_cutoff: float = 1.5,
        max_dynamic_distance: float = 1.5,
        max_candidates: int = 10000,
    ) -> List[Match]:
        """`list` of `Match`: Match the `list` of `Template` to one `~pyjess.Molecule`"""

        # killswitch is controlled by max_candidates. Internal default is currently 1000
        # killswitch serves to limit the iterations in cases where the template would be too general,
        # and the program would run in an almost endless loop

        jess = pyjess.Jess(
            templates
        )  # Create a Jess instance and use it to query a molecule (a PDB structure) against the stored templates:
        query = jess.query(
            molecule=molecule,
            rmsd_threshold=rmsd_threshold,
            distance_cutoff=distance_cutoff,
            max_dynamic_distance=max_dynamic_distance,
            max_candidates=max_candidates,
            best_match=True,
            ignore_chain=True,
        )  # query is pyjess.Query object which is an iterator over pyjess.Hits

        # ignore_chain=True disables checks for chain relationship
        # i.e. if two atoms are on a different chain in the template
        # they must also be on different chains in the target

        # best_match=True reports only the single best match between template and target
        # For this to make sense consider that:
        # A template is not encoded as coordinates, rather as a set of constraints.
        # For example, it would not contain the exact positions of THR and ASN atoms,
        # but instructions like
        # "Cα of ASN should be X angstrom away from the Cα of THR plus the allowed distance."

        # Multiple solutions = Matches to a template, satisfying all constraints may therefore exist
        # Jess produces matches to templates by looking for any combination of atoms,
        # residue_types, elements etc. and ANY positions which satisfy the constraints in the template

        # thus the solutions that Jess finds are NOT identical to the template at all
        # rather they are all possible solutions to the set constraints.
        # Solutions may completely differ from the template geometry
        # or atom composition if allowed by the set constraints.
        # by setting best_match=True we turn on filtering by rmsd to return only the best match
        # for every molecule template pair. Currently best_match is not exposed to the user.
        # This should be the only use case (I think)

        # Thus we hope to return the one solution to the constraints
        # which most closely resembles the original template - this is not guaranteed of course

        matches: List[Match] = []
        for hit in query:  # hit is pyjess.Hit
            matches.append(Match(hit=hit, pairwise_distance=distance_cutoff))

        return matches

    def _single_query_run(
        self,
        molecule: pyjess.Molecule,
        templates: List[Template],
        rmsd_threshold: float,
        distance_cutoff: float,
        max_dynamic_distance: float,
        max_candidates: int = 10000,
    ) -> List[Match]:

        matches = self._run_jess(
            molecule=molecule,
            templates=templates,
            rmsd_threshold=rmsd_threshold,
            distance_cutoff=distance_cutoff,
            max_dynamic_distance=max_dynamic_distance,
            max_candidates=max_candidates,
        )
        self._check_completeness(matches)

        return matches

    def _filter_matches(
        self,
        all_matches: Iterable[Tuple[pyjess.Molecule, List[Match]]],
        processed_molecules: Dict[pyjess.Molecule, List[Match]],
    ):

        # counter for all matches found for a given tempate size
        total_matches = 0

        # keep only matches predicted as correct
        if self.filter_matches:
            for molecule, matches in all_matches:
                for match in matches:
                    if match.predicted_correct:
                        total_matches += 1
                        processed_molecules[molecule].append(match)

        # all matches are added to processed_molecules
        else:
            for molecule, matches in all_matches:
                if matches:
                    processed_molecules[molecule].extend(matches)
                total_matches += len(matches)

        self.verbose_print(f"{total_matches} matches found!")
        self.verbose_print(f"{len(processed_molecules)} target structures processed!")

        return processed_molecules

    def run(
        self, molecules: List[pyjess.Molecule]
    ) -> Dict[pyjess.Molecule, List[Match]]:
        """
        Run the matcher against a `list` of query `~pyjess.Molecule` to search.

        Arguments:
            molecules: `list` of `~pyjess.Molecule` to search

        Note:
            Will first search against larger and then progressively smaller templates

        Note:
            Templates with fewer than 3 residues are skipped

        Returns:
            `dict` of `~pyjess.Molecule` --> `list` of `Match`: Dictionary of query molecules as keys and all found matches as values.
        """

        processed_molecules: Dict[pyjess.Molecule, List[Match]] = (
            collections.defaultdict(list)
        )

        pool: DummyPool | ThreadPool = (
            DummyPool() if self.cpus == 1 else ThreadPool(self.cpus)
        )
        with pool, Progress(
            SpinnerColumn(),
            *Progress.get_default_columns(),
            TimeElapsedColumn(),
            "Structures {task.completed}/{task.total}",
            console=self.console,
        ) as progress:

            for template_size in self.template_effective_sizes:

                if template_size < 3 and not self.match_small_templates:
                    # skip the rest of the for-loop
                    continue

                templates = self.templates_by_effective_size[template_size]

                rmsd, distance, max_dynamic_distance = self._get_jess_parameters(
                    template_size
                )

                self.verbose_print(
                    f"Now matching query structure(s) to template of size {template_size}"
                )
                self.verbose_print(
                    f"jess parameters are: {rmsd} {distance} {max_dynamic_distance}"
                )

                ################# Running Jess ###################################

                if (
                    self.skip_smaller_hits
                ):  # only search the molecule if no larger hits have been found for it already
                    query_molecules = [
                        molecule
                        for molecule in molecules
                        if molecule not in processed_molecules
                    ]
                else:
                    query_molecules = molecules

                _single_query_run = functools.partial(
                    self._single_query_run,
                    templates=templates,
                    rmsd_threshold=rmsd,
                    distance_cutoff=distance,
                    max_dynamic_distance=max_dynamic_distance,
                )

                task_id = progress.add_task(
                    description=f"[green]Searching with {template_size}-residue templates ...",
                    total=len(query_molecules),
                )

                all_matches = []
                # NOTE we cannot use imap_unordered here ... mb work on this
                for mol, result in zip(
                    query_molecules, pool.imap(_single_query_run, query_molecules)  # type: ignore
                ):
                    all_matches.append((mol, result))
                    progress.advance(task_id)

                processed_molecules = self._filter_matches(
                    all_matches=all_matches,
                    processed_molecules=processed_molecules,
                )

        return processed_molecules

    def run_single(self, molecule: pyjess.Molecule) -> List[Match]:
        """
        Run the matcher against a single query `~pyjess.Molecule`.

        Argument:
            molecule: `~pyjess.Molecule` to search

        Returns;
            `list`: of `Match` found for the query `~pyjess.Molecule`
        """
        return self.run([molecule])[molecule]
