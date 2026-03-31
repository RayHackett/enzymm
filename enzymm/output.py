from __future__ import annotations

import csv
import io
from itertools import chain
from typing import (
    overload,
    Iterable,
    Any,
    TextIO,
    List,
    Literal,
    Type,
    Dict,
    TYPE_CHECKING,
    Tuple,
)

import pyjess

if TYPE_CHECKING:
    import polars
    from enzymm.jess_run import Match

from enzymm import __version__
from enzymm.template import AnnotatedTemplate, AnnotatedResidue


TableKind = Literal["simple", "full", "residue"]


class _BaseTable:

    def __init__(
        self,
        matches: List[Match],
    ):
        self.matches = matches

    @classmethod
    def columns(cls) -> List[str]:
        """
        Return List of column names.
        """
        raise NotImplementedError

    @classmethod
    def from_match(cls, match: Match):
        """
        Generate the table from a single Match
        """
        raise NotImplementedError

    @classmethod
    def from_matches(cls, matches: Iterable[Match]) -> _BaseTable:
        """
        Generate the table from a list of Match
        """
        return cls(list(matches))

    @classmethod
    def _polars_schema(cls) -> Dict[str, Any]:
        """
        Return Polars dtypes keyed by column name.
        """
        raise NotImplementedError

    def dumps(self, header: bool = False) -> str:
        """
        Dump `Match` to a string. Calls `Match.dump()`

        Arguments:
            header: if a header line should be dumped to the string too.
        """
        buffer = io.StringIO()
        self.write_tsv(buffer, header=header)
        return (
            buffer.getvalue()
        )  # returns entire content temporary file object as a string

    def write_tsv(self, file: TextIO, header: bool = True) -> None:
        """
        Dump the Table to a '.tsv' like file.

        Arguments:
            file: `file-like` object to write to
            header: `bool` If a header line should be written too

        Note:
            Coordinate information is not written.
        """
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")

        if header:
            file.write(
                f"# Enzymm Version {__version__} running PyJess Version {pyjess.__version__}\n"
            )
            writer.writerow(self.columns())

        rows = (
            row.row_str()
            for row in chain.from_iterable(self.from_match(m) for m in self.matches)
        )
        writer.writerows(rows)

    def to_polars(self) -> polars.DataFrame:
        try:
            import polars
        except ImportError:
            raise RuntimeError("polars is not installed")

        rows = (
            row.row()
            for row in chain.from_iterable(self.from_match(m) for m in self.matches)
        )

        df = polars.DataFrame(
            rows,
            schema=self._polars_schema(),
            orient="row",
        )

        return df.select(self._polars_schema().keys())


class FullMatchTable(_BaseTable):

    @classmethod
    def columns(cls) -> List[str]:
        return [
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
            "number_of_ptm_residues_(template,reference)",
            "total_reference_residues",
        ]

    @classmethod
    def _polars_schema(cls):
        try:
            import polars
        except ImportError:
            raise RuntimeError("polars is not installed")

        return {
            "query_id": polars.Utf8,
            "pairwise_distance": polars.Float64,
            "match_index": polars.UInt32,
            "template_pdb_id": polars.Utf8,
            "template_pdb_chains": polars.List(polars.Utf8),
            "template_cluster_id": polars.UInt8,
            "template_cluster_member": polars.UInt8,
            "template_cluster_size": polars.UInt8,
            "template_effective_size": polars.UInt8,
            "template_dimension": polars.UInt8,
            "template_mcsa_id": polars.UInt16,
            "template_uniprot_id": polars.Utf8,
            "template_ec": polars.List(polars.Utf8),
            "template_cath": polars.List(polars.Utf8),
            "template_multimeric": polars.Boolean,
            "query_multimeric": polars.Boolean,
            "query_atom_count": polars.UInt32,
            "query_residue_count": polars.UInt32,
            "rmsd": polars.Float64,
            "log_evalue": polars.Float64,
            "orientation": polars.Float64,
            "preserved_order": polars.Boolean,
            "completeness": polars.Boolean,
            "predicted_correct": polars.Boolean,
            "matched_residues": polars.List(
                polars.Struct(
                    [
                        polars.Field("res", polars.Utf8),
                        polars.Field("chain", polars.Utf8),
                        polars.Field("num", polars.UInt32),
                    ]
                )
            ),
            "number_of_mutated_residues": polars.UInt8,
            "number_of_side_chain_residues": polars.Array(polars.UInt8, 2),
            "number_of_metal_ligands": polars.Array(polars.UInt8, 2),
            "number_of_ptm_residues": polars.Array(polars.UInt8, 2),
            "total_reference_residues": polars.UInt8,
        }

    class Row:
        def __init__(self, match: Match):
            self.match = match

        @staticmethod
        def to_res_struct_list(
            data: List[Tuple[str, str, str]],
        ) -> List[Dict[str, str | int]]:
            return [{"res": a, "chain": b, "num": int(c)} for a, b, c in data]

        def row(self) -> List[Any]:
            t = self.match.hit.template()
            c = t.cluster

            base = [
                self.match.hit.molecule().id,
                self.match.pairwise_distance,
                self.match.index,
                t.pdb_id,
                list({r.chain_id for r in t.residues}),
                c.id,
                c.member,
                c.size,
                t.effective_size,
                t.dimension,
                t.mcsa_id,
                t.uniprot_id,
                t.ec or [],
                t.cath or [],
                t.multimeric,
                self.match.multimeric,
                self.match.query_atom_count,
                self.match.query_residue_count,
                self.match.hit.rmsd,
                self.match.hit.log_evalue,
                self.match.orientation,
                self.match.preserved_resid_order,
                self.match.complete,
                self.match.predicted_correct,
                self.to_res_struct_list(self.match.matched_residues),
            ]

            if isinstance(t, AnnotatedTemplate):
                base.extend(
                    [
                        t.number_of_mutated_residues,
                        t.number_of_side_chain_residues,
                        t.number_of_metal_ligands,
                        t.number_of_ptm_residues,
                        t.total_reference_residues,
                    ]
                )
            else:
                base.extend([None, None, None, None, None])

            return base

        def row_str(self) -> List[str]:
            t = self.match.hit.template()
            c = t.cluster

            base = [
                str(self.match.hit.molecule().id),
                str(self.match.pairwise_distance),
                str(self.match.index),
                str(t.pdb_id or ""),
                ",".join({r.chain_id for r in t.residues}),
                str(c.id if c else ""),
                str(c.member if c else ""),
                str(c.size if c else ""),
                str(t.effective_size),
                str(t.dimension),
                str(t.mcsa_id or ""),
                str(t.uniprot_id or ""),
                ",".join(t.ec or []),
                ",".join(t.cath or []),
                str(t.multimeric),
                str(self.match.multimeric),
                str(self.match.query_atom_count),
                str(self.match.query_residue_count),
                round(self.match.hit.rmsd, 5),
                round(self.match.hit.log_evalue, 5),
                round(self.match.orientation, 5),
                str(self.match.preserved_resid_order),
                str(self.match.complete),
                str(self.match.predicted_correct or ""),
                ",".join("_".join(x) for x in self.match.matched_residues),
            ]

            if isinstance(t, AnnotatedTemplate):
                base.extend(
                    [
                        str(t.number_of_mutated_residues),
                        ",".join(map(str, t.number_of_side_chain_residues)),
                        ",".join(map(str, t.number_of_metal_ligands)),
                        ",".join(map(str, t.number_of_ptm_residues)),
                        str(t.total_reference_residues),
                    ]
                )
            else:
                base.extend(["", "", "", "", ""])

            return base

    @classmethod
    def from_match(cls, match: Match) -> Iterable[FullMatchTable.Row]:
        yield cls.Row(match)


class SimpleMatchTable(_BaseTable):

    @classmethod
    def columns(cls) -> List[str]:
        return [
            "query_id",
            "pdb_id",
            "pdb_chains",
            "effective_size",
            "size_dimension",
            "mcsa_id",
            "uniprot_id",
            "ec",
            "cath",
            "rmsd",
            "orientation_deg",
            "preserved_order",
            "metal_ligands_(template,reference)",
            "total_reference_residues",
            "matched_residues",
            "predicted_correct",
        ]

    @classmethod
    def _polars_schema(cls):
        try:
            import polars
        except ImportError:
            raise RuntimeError("polars is not installed")

        return {
            "query_id": polars.Utf8,
            "pdb_id": polars.Utf8,
            "pdb_chains": polars.List(polars.Utf8),
            "effective_size": polars.UInt8,
            "size_dimension": polars.UInt8,
            "mcsa_id": polars.UInt16,
            "uniprot_id": polars.Utf8,
            "ec": polars.List(polars.Utf8),
            "cath": polars.List(polars.Utf8),
            "rmsd": polars.Float64,
            "orientation_deg": polars.Float64,
            "preserved_order": polars.Boolean,
            "metal_ligands": polars.Array(polars.UInt8, 2),
            "total_reference_residues": polars.UInt8,
            "matched_residues": polars.List(
                polars.Struct(
                    [
                        polars.Field("res", polars.Utf8),
                        polars.Field("chain", polars.Utf8),
                        polars.Field("num", polars.UInt32),
                    ]
                )
            ),
            "predicted_correct": polars.Boolean,
        }

    class Row:
        def __init__(self, match: Match):
            self.match = match

        @staticmethod
        def to_res_struct_list(
            data: List[Tuple[str, str, str]]
        ) -> List[Dict[str, str | int]]:
            return [{"res": a, "chain": b, "num": int(c)} for a, b, c in data]

        def row(self) -> List[Any]:
            t = self.match.hit.template()
            row = [
                str(self.match.hit.molecule().id),
                t.pdb_id,
                list({r.chain_id for r in t.residues}),
                t.effective_size,
                t.dimension,
                t.mcsa_id,
                t.uniprot_id,
                t.ec or [],
                t.cath or [],
                self.match.hit.rmsd,
                self.match.orientation,
                self.match.preserved_resid_order,
            ]

            if isinstance(t, AnnotatedTemplate):
                row.extend(
                    [
                        t.number_of_metal_ligands,
                        t.total_reference_residues,
                    ]
                )
            else:
                row.extend(
                    [
                        None,
                        None,
                        self.to_res_struct_list(self.match.matched_residues),
                        self.match.predicted_correct,
                    ]
                )
            return row

        def row_str(self) -> List[str]:
            t = self.match.hit.template()

            row = [
                str(self.match.hit.molecule().id),
                str(t.pdb_id or ""),
                ",".join({r.chain_id for r in t.residues}),
                str(t.effective_size),
                str(t.dimension),
                str(t.mcsa_id or ""),
                str(t.uniprot_id or ""),
                ",".join(t.ec or []),
                ",".join(t.cath or []),
                round(self.match.hit.rmsd, 5),
                round(self.match.orientation, 5),
                str(self.match.preserved_resid_order),
            ]
            if isinstance(t, AnnotatedTemplate):
                row.extend(
                    [
                        (
                            str(t.number_of_metal_ligands[0])
                            if t.number_of_metal_ligands
                            else ""
                        ),
                        str(t.total_reference_residues),
                    ]
                )
            else:
                row.extend(
                    [
                        "",
                        "",
                        ", ".join("_".join(x) for x in self.match.matched_residues),
                        str(self.match.predicted_correct or ""),
                    ]
                )
            return row

    @classmethod
    def from_match(cls, match: Match) -> Iterable[SimpleMatchTable.Row]:
        yield cls.Row(match)


class MatchResidueTable(_BaseTable):

    @classmethod
    def columns(cls) -> List[str]:
        return [
            "query_id",
            "match_index",
            "query_residue",
            "template_pdb_id",
            "template_residue",
            "reference_pdb_id",
            "reference_residue",
            "roles",
        ]

    def _polars_schema(cls):
        try:
            import polars
        except ImportError:
            raise RuntimeError("polars is not installed")

        return {
            "query_id": polars.Utf8,
            "match_index": polars.UInt32,
            "query_residue": polars.Struct(
                [
                    polars.Field("res", polars.Utf8),
                    polars.Field("chain", polars.Utf8),
                    polars.Field("num", polars.UInt32),
                ]
            ),
            "template_pdb_id": polars.Utf8,
            "template_residue": polars.Struct(
                [
                    polars.Field("res", polars.Utf8),
                    polars.Field("chain", polars.Utf8),
                    polars.Field("num", polars.UInt32),
                ]
            ),
            "reference_pdb_id": polars.Utf8,
            "referece_residue": polars.Struct(
                [
                    polars.Field("res", polars.Utf8),
                    polars.Field("chain", polars.Utf8),
                    polars.Field("num", polars.UInt32),
                ]
            ),
            "roles": polars.List(polars.Utf8),
        }

    class Row:
        def __init__(
            self,
            q_res: Tuple[str, str, str],
            t_res: AnnotatedResidue,
            match: Match,
        ):
            self.q_res = q_res
            self.t_res = t_res
            self.match = match

        @staticmethod
        def to_res_struct(data: Tuple[str, str, str | int]) -> Dict[str, str | int]:
            return {"res": data[0], "chain": data[1], "num": int(data[2])}

        def row(self) -> List[Any]:
            row = [
                str(self.match.hit.molecule().id),
                self.match.index,
                self.to_res_struct((self.q_res[0], self.q_res[1], self.q_res[2])),
                self.match.hit.template().pdb_id,
                self.to_res_struct(
                    (self.t_res.name, self.t_res.chain_id, self.t_res.number)
                ),
            ]

            if isinstance(self.match.hit.template(), AnnotatedTemplate):
                row.extend(
                    [
                        self.t_res.reference_pdb.pdb_id,
                        self.to_res_struct(
                            (
                                self.t_res.reference_residue.name,
                                self.t_res.reference_pdb.chain_id,
                                self.t_res.reference_residue.auth_number,
                            )
                        ),
                        self.t_res.roles_summary,
                    ]
                )
            else:
                row.extend([None, None, None])
            return row

        def row_str(self) -> List[str]:
            row = [
                str(self.match.hit.molecule().id),
                str(self.match.index),
                "_".join((self.q_res[0], self.q_res[1], self.q_res[2])),
                self.match.hit.template().pdb_id,
                "_".join(
                    (self.t_res.name, self.t_res.chain_id, str(self.t_res.number))
                ),
                self.t_res.reference_pdb.pdb_id,
                "_".join(
                    (
                        self.t_res.reference_residue.name,
                        self.t_res.reference_pdb.chain_id,
                        str(self.t_res.reference_residue.auth_number),
                    )
                ),
                ",".join(self.t_res.roles_summary),
            ]
            if isinstance(self.match.hit.template(), AnnotatedTemplate):
                row.extend(
                    [
                        self.t_res.reference_pdb.pdb_id,
                        "_".join(
                            (
                                self.t_res.reference_residue.name,
                                self.t_res.reference_pdb.chain_id,
                                str(self.t_res.reference_residue.auth_number),
                            )
                        ),
                        ",".join(self.t_res.roles_summary),
                    ]
                )
            else:
                row.extend(["", "", ""])
            return row

    @classmethod
    def from_match(cls, match: Match) -> Iterable[MatchResidueTable.Row]:
        t = match.hit.template()
        for q_res, t_res in zip(
            match.matched_residues,
            t.residues,
            strict=True,
        ):  # ty:ignore[no-matching-overload]
            yield cls.Row(q_res, t_res, match)


class Tables:
    # a registry dict for my different table types
    _map: dict[TableKind, Type[_BaseTable]] = {
        "simple": SimpleMatchTable,
        "full": FullMatchTable,
        "residue": MatchResidueTable,
    }

    @overload
    @classmethod
    def create(
        cls, kind: Literal["full"], matches: Iterable[Match]
    ) -> FullMatchTable: ...

    @overload
    @classmethod
    def create(
        cls, kind: Literal["simple"], matches: Iterable[Match]
    ) -> SimpleMatchTable: ...

    @overload
    @classmethod
    def create(
        cls, kind: Literal["residue"], matches: Iterable[Match]
    ) -> MatchResidueTable: ...

    @classmethod
    def create(cls, kind: TableKind, matches: Iterable[Match]) -> _BaseTable:
        return cls._map[kind].from_matches(matches)
