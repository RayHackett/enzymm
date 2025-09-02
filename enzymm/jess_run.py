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
from multiprocessing.pool import ThreadPool
import functools
from typing import (
    List,
    Tuple,
    Dict,
    Optional,
    IO,
    Iterable,
    Sequence,
    Callable,
    ClassVar,
)
from functools import cached_property
from dataclasses import dataclass, field
from pathlib import Path

import rich
import pyjess
from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn, TaskID, Task
from readerwriterlock import rwlock

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

from enzymm import __version__
from enzymm.template import Template, AnnotatedTemplate, Vec3, check_template
from enzymm.utils import chunks, ranked_argsort, DummyPool
from enzymm.utils import PROTEINOGENIC_AMINO_ACIDS, SPECIAL_AMINO_ACIDS

__all__ = [
    "LogisticRegressionModel",
    "ModelEnsemble",
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
        coef: `list` of `floats` beta coefficients
        intercept: `float` beta0 intercept
        threshold: `float` optimal threshold for this model
    """

    coef: List[float]
    intercept: float
    threshold: float

    @cached_property
    def logit_threshold(self) -> float:
        """Precomputed logit value of the threshold"""
        return math.log(self.threshold / (1 - self.threshold))

    def __call__(
        self,
        rmsd: float,
        orientation: float,
    ) -> bool:
        """
        Make a prediction with the Logistc Regression Model based on RMSD and residue orientation.

        Attributes:
            rmsd `float`: RMSD value of the match
            orientation `float`: residue orientation of the match

        Returns:
            `bool`: If the match is predicted as correct or not
        """

        predicted_logit = (
            self.intercept + self.coef[0] * rmsd + self.coef[1] * orientation
        )
        # We compare logits
        return predicted_logit >= self.logit_threshold


@dataclass(frozen=True)
class ModelEnsemble:
    """
    Ensemble of Models which each produce a binary prediction. The ensemble takes a majority vote.

    Attributes:
        ensemble: `Dict[int, Dict[float, List[Callable[..., bool]]]] Dictonary of template_effective_size of Dictonaries of pairwise_distance of a List of callable models.
        min_true_template_size: `int` Minimum effective size of a template to be considered always correct.
        minimum_effective_size: `int` Smallest template_effective_size for which there are models. Smaller template will be treated as if they had 3 residues.

    Note:
        The `ensemble` dictionary should cover at least 3 and 4 residue matches
        at pairwise distances in the "usual" range - about 0.7 to 2.0A!
        The call method will raise an error otherwise!

    """

    ensemble: Dict[int, Dict[float, List[Callable[..., bool]]]]
    min_true_template_size: int
    minimum_effective_size: int

    @classmethod
    def from_json(
        cls,
        json_file: IO[str],
        model_cls: Callable,
    ) -> ModelEnsemble:
        """Build an ensemble model directly from an open JSON file"""

        model_dict = json.load(json_file)

        ensemble: Dict[int, Dict[float, List[Callable[..., bool]]]] = {}
        for template_size, pairwise_dict in model_dict["match_size"].items():
            ensemble[int(template_size)] = {}
            for pairwise_distance, model_list in pairwise_dict[
                "pairwise_distance"
            ].items():
                ensemble[int(template_size)][float(pairwise_distance)] = [
                    model_cls(**param_dict)  # just unpacks the key, value pairs
                    for param_dict in model_list["model_list"]
                ]

        min(ensemble.keys())

        return cls(
            ensemble=ensemble,
            minimum_effective_size=min(ensemble.keys()),
            min_true_template_size=max(ensemble.keys()) + 1,
        )

    def number_of_models(
        self,
        *,
        template_effective_size: int,
        pairwise_distance: float,
    ) -> int:
        """
        Number of models for a given template_effective_size and pairwise_distance

        Attributes:
            template_effective_size: `int` Number of side chain residues in the template
            pairwise_distance: `float` Pairwise distance of the match

        Returns:
            `int`: Number of models
        """
        return len(self.ensemble[template_effective_size][pairwise_distance])

    def __call__(
        self,
        *,
        template_effective_size: int,
        pairwise_distance: float,
        model_kwargs: Dict[str, float],
    ) -> bool:
        """
        Make an ensemble prediction at a given template_effective_size and pairwise_distance

        Attributes:
            template_effective_size: `int` Number of side chain residues in the template
            pairwise_distance: `float` Pairwise distance of the match
            model_paramters: `float` Named floats to pass parameters to the individual models

        Returns:
            `bool`: Wether the match is predicted correct or false by the ensemble model
        """

        if template_effective_size >= self.min_true_template_size:
            # Matches with 5+ residues are considered true
            return True
        else:

            # treat templates with smaller effective sizes as if they had 3 residues
            if template_effective_size < self.minimum_effective_size:
                template_effective_size = 3

            predictions = []
            for model in self.ensemble[template_effective_size][pairwise_distance]:
                predictions.append(model(**model_kwargs))

            num_models = self.number_of_models(
                template_effective_size=template_effective_size,
                pairwise_distance=pairwise_distance,
            )

            # majority decision from all models
            return sum(predictions) >= (num_models + 1) // 2


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
    complete: bool = False
    pairwise_distance: float = 0
    index: int = 0
    ensemble_model: ClassVar[ModelEnsemble]

    def __reduce_ex__(self, protocol):
        return (
            functools.partial(
                Match,
                hit=self.hit,
                complete=self.complete,
                pairwise_distance=self.pairwise_distance,
                index=self.index,
                # we skip the _logistic_regression_models since it is a ClassVar
            ),
            (),
        )

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
            file.write(
                f"# Enzymm Version {__version__} running PyJess Version {pyjess.__version__}\n"
            )
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
        """`bool`: If the match is predicted as correct based the ensemble model"""
        return self.ensemble_model(
            pairwise_distance=self.pairwise_distance,
            template_effective_size=self.hit.template.effective_size,
            model_kwargs={"rmsd": self.hit.rmsd, "orientation": self.orientation},
        )

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


with resource_files(__package__).joinpath("data", "logistic_regression_models.json").open() as f:  # type: ignore
    Match.ensemble_model = ModelEnsemble.from_json(
        f,
        model_cls=LogisticRegressionModel,
    )


def load_molecules(
    molecule_paths: List[Path], conservation_cutoff: float = 0
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
        else:
            raise ValueError(
                f"Received an empty molecule from {molecule_path}. Is this file in PDB format?"
            ) from None

    if not molecules:
        raise FileNotFoundError("Received no molecules from -i or -l input!") from None

    return molecules


@dataclass
class QueryMolecule:
    molecule: pyjess.Molecule
    lock: rwlock.RWLockRead = field(default_factory=rwlock.RWLockRead)
    # the lock is similar to: threading.Lock = threading.Lock()
    hit_found: bool = False
    hit_size: int = 0


class StructuresColumn(rich.progress.ProgressColumn):
    def render(self, task: Task):
        job_batches = task.fields.get("job_batches", 1)
        return f"Structures {task.completed // job_batches}/{task.total // job_batches}"


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
        templates: Sequence[Template],
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
        )  # Dictionary of List of Template objects grouped by Template.effective_size as keys

        for template in templates:
            # skip smaller templates if match_small_templates is not set!
            if not self.match_small_templates and template.effective_size < 3:
                continue
            if check_template(
                template, warn=self.warn
            ):  # returns True if the Template passed all checks or if warn is set to False
                self.templates_by_effective_size[template.effective_size].append(
                    template
                )

        if self.verbose:
            template_number_dict: Dict[int, int] = {}
            for size, template_list in self.templates_by_effective_size.items():
                template_number_dict[size] = len(template_list)
            print(
                f"Templates by effective size: {collections.OrderedDict(sorted(template_number_dict.items()))}"
            )

        self.template_effective_sizes = list(self.templates_by_effective_size.keys())
        self.template_effective_sizes.sort(
            reverse=True
        )  # get a list of template_sizes in decending order

        # print a warning if match_small_templates was set.
        if self.warn and self.match_small_templates:
            smaller_sizes = [i for i in self.template_effective_sizes if i < 3]
            if smaller_sizes:
                small_templates = []
                for i in smaller_sizes:
                    small_templates.extend(self.templates_by_effective_size[i])

                warnings.warn(
                    f"{len(small_templates)} Templates with an effective size smaller than 3 defined sidechain residues were supplied.\nFor small templates Jess parameters for templates of 3 residues will be used."
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
        templates: Iterable[Template],
        rmsd_threshold: float = 2.0,
        distance_cutoff: float = 1.5,
        max_dynamic_distance: float = 1.5,
        max_candidates: Optional[int] = None,
    ) -> List[Match]:
        """`list` of `Match`: Match the `list` of `Template` to one `~pyjess.Molecule`"""

        # killswitch is controlled by max_candidates. Internal default is None
        # Which disabled the killswitch. by setting it to an integer like 1000 or 10000
        # Only that many matches (total matches not best matches!)
        # for a given query-template pair are returned
        # killswitch limits the iterations when the template would be too general,
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
            ignore_chain="residues",
        )  # query is pyjess.Query object which is an iterator over pyjess.Hits

        # ignore_chain=None will behave like the previous False - will observe relative chain assignments in the template.
        # ignore_chain="atoms" works like the previous True - might result in matches in which residues can be split across chains.
        # ignore_chain="residues" will check for chain membership only between residues. Atoms within a residue always belong to the same chain.

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

    def _filter_molecule_matches(
        self,
        all_matches: List[Match],
    ):

        # keep only matches predicted as correct
        if self.filter_matches:
            filtered_matches = []
            for match in all_matches:
                if match.predicted_correct:
                    filtered_matches.append(match)
            return filtered_matches

        # return unchanged
        else:
            return all_matches

    def _single_query_run(
        self,
        molecule: pyjess.Molecule,
        templates: Iterable[Template],
        rmsd_threshold: float,
        distance_cutoff: float,
        max_dynamic_distance: float,
        max_candidates: Optional[int] = None,
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

        return self._filter_molecule_matches(matches)

    def run(
        self, molecules: List[pyjess.Molecule]
    ) -> Dict[pyjess.Molecule, List[Match]]:
        """
        Run the matcher against a `list` of query `~pyjess.Molecule` to search.

        Arguments:
            molecules: `list` of `~pyjess.Molecule` to search

        Returns:
            `dict` of `~pyjess.Molecule` --> `list` of `Match`: Dictionary of query molecules as keys and all found matches as values.
        """
        processed_molecules: Dict[pyjess.Molecule, List[Match]] = (
            collections.defaultdict(list)
        )

        query_molecules = [QueryMolecule(molecule=mol) for mol in molecules]

        # batches of templates with only one size of templates within each batch
        # batches must be pure in terms of template size
        ordered_template_batches: List[Tuple[Tuple[pyjess.Template, ...], int]] = []
        for template_size in self.template_effective_sizes:
            for batch in chunks(
                iterable=self.templates_by_effective_size[template_size], n=100
            ):
                ordered_template_batches.append((batch, template_size))

        job_batches = []
        for template_batch, template_size in ordered_template_batches:

            rmsd, distance, max_dynamic_distance = self._get_jess_parameters(
                template_size
            )

            the_function = functools.partial(
                self._single_query_run,
                templates=template_batch,
                rmsd_threshold=rmsd,
                distance_cutoff=distance,
                max_dynamic_distance=max_dynamic_distance,
            )

            job_batches.append(the_function)

        def process(
            batch_partial: functools.partial,
            molecule: QueryMolecule,
            progress: Progress,
            task: TaskID,
        ):

            template_size = batch_partial.keywords["templates"][0].effective_size

            if self.skip_smaller_hits:
                with molecule.lock.gen_rlock():
                    if molecule.hit_found and molecule.hit_size > template_size:
                        progress.advance(task_id=task)
                        return []

            # call the partial function with the molecule
            results = list(batch_partial(molecule=molecule.molecule))

            if results and self.skip_smaller_hits:
                with molecule.lock.gen_wlock():
                    molecule.hit_found = True
                    molecule.hit_size = template_size

            progress.advance(task_id=task)

            return results

        pool: DummyPool | ThreadPool = (
            DummyPool() if self.cpus == 1 else ThreadPool(self.cpus)
        )

        with pool, Progress(
            SpinnerColumn(),
            *Progress.get_default_columns(),
            TimeElapsedColumn(),
            StructuresColumn(),
            console=self.console,
        ) as progress:

            task = progress.add_task(
                "Searching Structures",
                total=len(job_batches) * len(molecules),
                job_batches=len(job_batches),
            )
            _process = functools.partial(process, progress=progress, task=task)

            results = list(
                pool.starmap(_process, itertools.product(job_batches, query_molecules))
            )

        # this is to reverse the itertools product
        for (_, qmol), matches in zip(
            itertools.product(job_batches, query_molecules),
            results,
        ):
            if matches:
                # # NOTE
                # # due to parallelism some smaller results might have get computed
                # # and will not be skipped. These run at no extra time cost
                # # to clean those up too, uncomment this
                # if self.skip_smaller_hits:
                #     matches = [
                #         match for match in matches
                #         if match.hit.template.effective_size == qmol.hit_size
                #     ]
                processed_molecules[qmol.molecule].extend(matches)

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
