import io
import os
import unittest
import math
import textwrap

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

import pyjess
import numpy

import enzymm
from enzymm import __version__
from . import test_data
from enzymm import template, jess_run
from enzymm.output import Tables


class TestMatch(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(
            resource_files(enzymm).joinpath(
                "jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb",
            ),
            "r",
        ) as f:  # ty:ignore[no-matching-overload]
            cls.template_text1 = f.read()
            cls.template1 = template.AnnotatedTemplate.loads(
                cls.template_text1, warn=False
            )
            jess_1 = pyjess.Jess([cls.template1])

        with resource_files(test_data).joinpath("1AMY.pdb").open() as f:
            cls.molecule = pyjess.Molecule.load(f)  # ty:ignore[invalid-argument-type]

        query = jess_1.query(cls.molecule, 2, 1.5, 1.5, best_match=True)
        best_hits = list(query)

        cls.match1 = jess_run.Match(
            hit=best_hits[0], pairwise_distance=1.5, complete=True, index=0
        )

        with open(
            resource_files(enzymm).joinpath(
                "jess_templates_20230210/3_residues/results/csa3d_0415/csa3d_0415.cluster_1_1_2.1be0_A124-A175-A125-A289-A260.template.pdb",
            ),
            "r",
        ) as f:  # ty:ignore[no-matching-overload]
            template_text2 = f.read()
            cls.template2 = template.AnnotatedTemplate.loads(template_text2, warn=False)
            jess_2 = pyjess.Jess([cls.template2])

        query = jess_2.query(cls.molecule, 2.0, 1.0, 1.0, best_match=True)
        best_hits = list(query)

        cls.match2 = jess_run.Match(hit=best_hits[0], pairwise_distance=1.0)

        cls.maxDiff = None

    def test_match(self):
        self.assertEqual(self.match1.hit.molecule().id, "1AMY")
        self.assertEqual(self.match1.index, 0)
        self.assertEqual(self.match1.hit.template().pdb_id, self.template1.pdb_id)
        self.assertEqual(
            self.match1.hit.template().effective_size, self.template1.effective_size
        )
        self.assertEqual(self.match1.hit.template().dimension, self.template1.dimension)
        self.assertEqual(self.match1.hit.template().mcsa_id, self.template1.mcsa_id)
        self.assertEqual(
            self.match1.hit.template().uniprot_id, self.template1.uniprot_id
        )
        self.assertEqual(self.match1.hit.template().ec, self.template1.ec)
        self.assertEqual(self.match1.hit.template().cath, self.template1.cath)
        self.assertEqual(
            self.match1.hit.template().multimeric, self.template1.multimeric
        )
        self.assertIsNotNone(self.match1.multimeric)
        self.assertFalse(self.match1.multimeric)
        self.assertEqual(self.match1.query_atom_count, 3339)
        self.assertEqual(self.match1.query_residue_count, 403)
        self.assertAlmostEqual(self.match1.hit.rmsd, 0.32093143)
        self.assertAlmostEqual(self.match1.hit.log_evalue, -3.08424478)
        self.assertAlmostEqual(self.match1.orientation, 0.15327054322)
        self.assertIsNotNone(self.match1.predicted_correct)
        self.assertTrue(self.match1.predicted_correct)
        self.assertEqual(
            self.match1.template_vector_list,
            [res.orientation_vector for res in self.template1.residues],
        )
        expected = [
            template.Vec3(
                x=0.2290067979141952, y=-0.3853409610281773, z=0.377114677867322
            ),
            template.Vec3(
                x=0.4249816660862038, y=-0.21966898402981627, z=-0.3540863184957992
            ),
            template.Vec3(
                x=0.45459385444007694, y=-0.34869961601989985, z=0.10687378206512577
            ),
            template.Vec3(
                x=-0.8733960645698886, y=0.2563504028143271, z=-0.9840695023070225
            ),
            template.Vec3(
                x=-0.510183600042339, y=-0.1958417994791759, z=0.18963368325429997
            ),
        ]

        actual = self.match1.match_vector_list
        self.assertEqual(len(actual), len(expected))

        for a, e in zip(actual, expected):
            self.assertTrue(math.isclose(a.x, e.x, rel_tol=1e-9, abs_tol=1e-9))
            self.assertTrue(math.isclose(a.y, e.y, rel_tol=1e-9, abs_tol=1e-9))
            self.assertTrue(math.isclose(a.z, e.z, rel_tol=1e-9, abs_tol=1e-9))

        self.assertTrue(self.match1.preserved_resid_order)
        self.assertTrue(self.match1.complete)
        self.assertEqual(
            self.match1.matched_residues,
            [
                ("GLU", "A", "204"),
                ("ASP", "A", "87"),
                ("ASP", "A", "179"),
                ("HIS", "A", "288"),
                ("ASP", "A", "289"),
            ],
        )
        self.assertEqual(self.match2.hit.template().pdb_id, self.template2.pdb_id)
        self.assertEqual(
            self.match2.hit.template().effective_size, self.template2.effective_size
        )
        self.assertEqual(self.match2.hit.template().dimension, self.template2.dimension)
        self.assertEqual(self.match2.hit.template().mcsa_id, self.template2.mcsa_id)
        self.assertEqual(
            self.match2.hit.template().uniprot_id, self.template2.uniprot_id
        )
        self.assertEqual(self.match2.hit.template().ec, self.template2.ec)
        self.assertEqual(self.match2.hit.template().cath, self.template2.cath)
        self.assertEqual(
            self.match2.hit.template().multimeric, self.template2.multimeric
        )
        self.assertIsNotNone(self.match2.multimeric)
        self.assertFalse(self.match2.multimeric)
        self.assertAlmostEqual(self.match2.hit.rmsd, 1.7353479120)
        self.assertAlmostEqual(self.match2.hit.log_evalue, 1.8517964201)
        self.assertAlmostEqual(self.match2.orientation, 1.6503123465442575)
        self.assertIsNotNone(self.match2.predicted_correct)
        self.assertFalse(self.match2.predicted_correct)
        self.assertEqual(
            self.match2.template_vector_list,
            [res.orientation_vector for res in self.template2.residues],
        )
        self.assertFalse(self.match2.preserved_resid_order)
        self.assertFalse(self.match2.complete)
        self.assertEqual(
            self.match2.matched_residues,
            [("TRP", "A", "38"), ("HIS", "A", "288"), ("ASP", "A", "289")],
        )

    def test_tfm(self):
        # test the transformation matrix
        expected_tfm = numpy.array(
            [
                [0.90385903, 0.39426385, 0.16611706, 18.85708649],
                [0.3477188, -0.45075881, -0.82213632, 90.1209215],
                [-0.2492599, 0.80085736, -0.54451537, -21.38900405],
                [
                    0,
                    0,
                    0,
                    1,
                ],
            ]
        )
        got_tfm = numpy.asarray(self.match1.hit.transformation)
        self.assertTrue(numpy.allclose(expected_tfm, got_tfm))
        self.assertTrue(
            numpy.allclose(
                numpy.asarray(self.match1.hit.inverse_transformation),
                numpy.linalg.inv(expected_tfm),
            )
        )

        tfm_mol = self.molecule.transform(got_tfm)  # ty:ignore[invalid-argument-type]
        self.assertAlmostEqual(tfm_mol[0].x, 46.59270033179604)
        self.assertAlmostEqual(tfm_mol[0].y, 55.9905432060407)
        self.assertAlmostEqual(tfm_mol[0].z, 6.5389174388503974)

    def test_match_tsv_write(self):
        buffer = io.StringIO()
        tbl = Tables.create(kind="full", matches=[self.match1])

        tbl.write_tsv(file=buffer, header=False)
        with resource_files(test_data).joinpath("results.tsv").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_match_dumps(self):
        expected = textwrap.dedent(
            f"""
            # Enzymm Version {__version__} running PyJess Version {pyjess.__version__}
            query_id	pairwise_distance	match_index	template_pdb_id	template_pdb_chains	template_cluster_id	template_cluster_member	template_cluster_size	template_effective_size	template_dimension	template_mcsa_id	template_uniprot_id	template_ec	template_cath	template_multimeric	query_multimeric	query_atom_count	query_residue_count	rmsd	log_evalue	orientation	preserved_order	completeness	predicted_correct	matched_residues	number_of_mutated_residues	number_of_side_chain_residues_(template,reference)	number_of_metal_ligands_(template,reference)	number_of_ptm_residues_(template,reference)	total_reference_residues
            1AMY	1.5	0	1uh3	A	1	1	1	5	5	285	Q60053	3.2.1.10,3.2.1.135	2.60.40.10,2.60.40.1180,3.20.20.80	False	False	3339	403	0.32093	-3.08424	0.15327	True	True	True	GLU_A_204,ASP_A_87,ASP_A_179,HIS_A_288,ASP_A_289	0	5,5	0,0	0,0	5
            """
        )
        self.maxDiff = None
        self.assertMultiLineEqual(
            self.match1.dumps(header=True).strip(), expected.strip()
        )

    def test_match_dump2pdb(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=False)
        with resource_files(test_data).joinpath("1AMY_matches.pdb").open() as f:
            self.assertMultiLineEqual(buffer.getvalue().strip(), f.read().strip())

    def test_match_dump2pdb_transformed(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=True)
        with resource_files(test_data).joinpath(
            "1AMY_matches_template.pdb"
        ).open() as f:
            self.assertMultiLineEqual(buffer.getvalue(), f.read())

    def test_query_dump(self):
        buffer = io.StringIO()
        self.match1.dump_query(buffer, transform=False)

        expected = textwrap.dedent(
            """
            REMARK MOLECULE_ID 1AMY
            REMARK QUERY COORDINATE FRAME
            ATOM      1  N   GLN A   1       6.240  48.686  17.460  1.00 27.79           N
            ATOM      2  CA  GLN A   1       5.440  49.851  17.773  1.00 16.62           C
            ATOM      3  C   GLN A   1       6.628  50.721  18.086  1.00 14.24           C
            ATOM      4  O   GLN A   1       7.313  50.396  19.052  1.00 11.61           O
            """
        )
        buffer.seek(0)
        first5 = "".join(buffer.getvalue().splitlines(keepends=True)[:6])
        self.assertMultiLineEqual(first5.strip(), expected.strip())

    def test_template_dump(self):
        buffer = io.StringIO()
        self.match1.dump_template(buffer, transform=True)

        expected = textwrap.dedent(
            """
            REMARK TEMPLATE COORDINATE FRAME
            REMARK TEMPLATE
            REMARK CLUSTER 1_1_1
            REMARK REPRESENTING 93 CATALYTIC SITES
            REMARK ID 1uh3_A396-A262-A356-A471-A472
            REMARK MCSA_ID 285
            REMARK PDB_ID 1uh3
            REMARK UNIPROT_ID Q60053
            REMARK EC 3.2.1.10,3.2.1.135
            REMARK CATH 2.60.40.10,2.60.40.1180,3.20.20.80
            REMARK ENZYME alpha-amylase I(E.C.3.2.1.1)
            REMARK EXPERIMENTAL_METHOD X-ray diffraction
            REMARK RESOLUTION 2.6
            REMARK ORGANISM_NAME Thermoactinomyces vulgaris
            REMARK ORGANISM_ID 2026
            REMARK PDB_ID 5-residues-1uh3_A396-A262-A356-A471-A472_cluster_1-1-1
            ATOM      3  CD  GLU A 396      53.233  36.186  11.266 E     3.20
            ATOM      3  OE1 GLU A 396      52.490  35.259  11.650 E     3.20
            ATOM      3  OE2 GLU A 396      54.381  36.361  11.722 E     3.20
            ATOM      3  CG  ASP A 262      49.175  39.646  17.664 DE    1.90
            ATOM      3  OD1 ASP A 262      50.360  40.036  17.776 DE    1.90
            ATOM      3  OD2 ASP A 262      48.818  38.741  16.875 DE    1.90
            ATOM      3  CG  ASP A 356      50.337  34.680  15.420 DE    2.27
            ATOM      3  OD1 ASP A 356      50.120  33.502  15.786 DE    2.27
            ATOM      3  OD2 ASP A 356      51.480  35.175  15.320 DE    2.27
            ATOM      0  CG  HIS A 471      56.763  40.562  17.238 H     3.77
            ATOM      8  ND1 HIS A 471      56.123  41.172  16.180 H     3.77
            ATOM      8  CD2 HIS A 471      56.129  39.378  17.430 H     3.77
            ATOM      3  CG  ASP A 472      58.152  37.388  14.843 D     1.23
            ATOM      3  OD1 ASP A 472      58.311  37.254  16.076 D     1.23
            ATOM      3  OD2 ASP A 472      57.076  37.105  14.265 D     1.23
            """
        )

        self.assertMultiLineEqual(buffer.getvalue().strip(), expected.strip())

    # TODO looks like pyjess.Hit has no equality comparison
    # def test_pickling(self):
    #     match1_picle = pickle.loads(pickle.dumps(self.match1))
    #     print(self.match1.hit.atoms())
    #     print(match1_picle.hit.atoms())
    #     self.assertEqual(self.match1.__dict__, match1_picle.__dict__)


class TestMatcher(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        template_files = [
            "jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/5_residues/results/csa3d_0045/csa3d_0045.cluster_1_1_1.2cxg_A227-A229-A257-A327-A328.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_2.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_2_2.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0421/csa3d_0421.cluster_1_1_2.1bf2_A229-A232-A230-A259-A375-A435-A510-A128.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0896/csa3d_0896.cluster_2_1_2.2qy1_A135-A179-A137-A230-A175.template.pdb",
        ]

        template_texts = []
        for path in template_files:
            with open(
                resource_files(enzymm).joinpath(path), "r"
            ) as f:  # ty:ignore[no-matching-overload]
                template_texts.append(f.read())

        cls.template_list = []
        for txt in template_texts:
            cls.template_list.append(template.AnnotatedTemplate.loads(txt))

        res5_templates = list(
            template.load_templates(
                template_dir=resource_files(enzymm).joinpath(
                    "jess_templates_20230210/5_residues/results/csa3d_0285/"
                ),  # ty:ignore[invalid-argument-type]
            )
        )
        res4_templates = list(
            template.load_templates(
                template_dir=resource_files(enzymm).joinpath(
                    "jess_templates_20230210/4_residues/results/csa3d_0285"
                ),  # ty:ignore[invalid-argument-type]
            )
        )
        res3_templates = list(
            template.load_templates(
                template_dir=resource_files(enzymm).joinpath(
                    "jess_templates_20230210/3_residues/results/csa3d_0344"
                ),  # ty:ignore[invalid-argument-type]
            )
        )

        res5_and_res4_templates = res5_templates + res4_templates
        with_smaller_templates = res5_and_res4_templates + res3_templates

        cls.template_matcher = jess_run.Matcher(
            templates=res5_and_res4_templates, cpus=2, warn=True
        )
        cls.template_matcher2 = jess_run.Matcher(
            templates=res5_and_res4_templates, skip_smaller_hits=True, cpus=1
        )
        cls.template_matcher3 = jess_run.Matcher(
            templates=res5_and_res4_templates,
            skip_smaller_hits=True,
        )
        with cls.assertWarns(cls, Warning):
            cls.template_matcher4 = jess_run.Matcher(
                templates=with_smaller_templates,
                match_small_templates=True,
                warn=True,
                cpus=-1,
            )

        cls.unfiltered_matcher = jess_run.Matcher(
            templates=cls.template_list,
            jess_params={
                3: {"rmsd": 2, "distance": 1.2, "max_dynamic_distance": 1.2},
                4: {"rmsd": 2, "distance": 1.7, "max_dynamic_distance": 1.7},
                5: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                6: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                8: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
            },
            filter_matches=False,
        )

        cls.filtered_matcher = jess_run.Matcher(
            templates=cls.template_list,
            jess_params={
                3: {"rmsd": 2, "distance": 1.2, "max_dynamic_distance": 1.2},
                4: {"rmsd": 2, "distance": 1.7, "max_dynamic_distance": 1.7},
                5: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                6: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                8: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
            },
            filter_matches=True,
        )

        cls.filtered_matcher_different_cutoffs = jess_run.Matcher(
            templates=cls.template_list,
            jess_params={
                3: {"rmsd": 2, "distance": 0.6, "max_dynamic_distance": 0.6},
                4: {"rmsd": 2, "distance": 2.5, "max_dynamic_distance": 2.5},
                5: {"rmsd": 2, "distance": 0.6, "max_dynamic_distance": 0.6},
                6: {"rmsd": 2, "distance": 0.6, "max_dynamic_distance": 0.6},
                7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                8: {"rmsd": 2, "distance": 2.5, "max_dynamic_distance": 2.5},
            },
            filter_matches=True,
            skip_smaller_hits=False,
        )

        cls.unfiltered_matcher_different_cutoffs = jess_run.Matcher(
            templates=cls.template_list,
            jess_params={
                3: {"rmsd": 2, "distance": 0.6, "max_dynamic_distance": 0.6},
                4: {"rmsd": 2, "distance": 2.5, "max_dynamic_distance": 2.5},
                5: {"rmsd": 2, "distance": 0.6, "max_dynamic_distance": 0.6},
                6: {"rmsd": 2, "distance": 0.6, "max_dynamic_distance": 0.6},
                7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
                8: {"rmsd": 2, "distance": 2.5, "max_dynamic_distance": 2.5},
            },
            filter_matches=False,
            skip_smaller_hits=False,
        )

        with resource_files(test_data).joinpath("1AMY.pdb").open() as f:
            cls.molecule = pyjess.Molecule.load(f)  # ty:ignore[invalid-argument-type]

        with resource_files(test_data).joinpath(
            "AF-P0DUB6-F1-model_v6.cif"
        ).open() as f:
            cls.molecule2 = pyjess.Molecule.load(f)  # ty:ignore[invalid-argument-type]
            cls.molecule3 = cls.molecule2.conserved(80)

    def test_init(self):
        jess_params = {
            3: {"rmsd": 2, "distance": 0.9, "max_dynamic_distance": 0.9},
            4: {"rmsd": 2, "distance": 1.7, "max_dynamic_distance": 1.7},
            5: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
            6: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
            7: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
            8: {"rmsd": 2, "distance": 2.0, "max_dynamic_distance": 2.0},
        }

        # Testing jess param defaults are stored correctly
        self.assertEqual(self.template_matcher.jess_params, jess_params)
        # Testing if template effective size is detected correctly
        self.assertEqual(self.template_matcher.template_effective_sizes, [5, 4])
        # Test if cpu counts work correctly
        self.assertEqual(self.template_matcher.cpus, 2)
        self.assertEqual(self.template_matcher2.cpus, 1)
        self.assertEqual(self.template_matcher3.cpus, (os.cpu_count() or 1))
        self.assertEqual(self.template_matcher4.cpus, (os.cpu_count() or 1) - 1)

        # check for duplicated templates in Matcher
        with self.assertRaises(ValueError):
            jess_run.Matcher(templates=self.template_list + self.template_list)

    def test_Matcher_run(self):
        processed_molecules_1 = self.template_matcher.run(
            molecules=[self.molecule, self.molecule2]
        )
        processed_molecules_2 = self.template_matcher2.run(
            molecules=[self.molecule, self.molecule3]
        )
        processed_molecules_3 = self.template_matcher3.run(
            molecules=[self.molecule, self.molecule3]
        )

        self.assertEqual(len(processed_molecules_1[self.molecule]), 2)
        self.assertEqual(len(processed_molecules_1[self.molecule2]), 2)
        #  processed_molecules_2 has skip smaller hits and 1 thread
        self.assertEqual(len(processed_molecules_2[self.molecule]), 1)
        self.assertEqual(len(processed_molecules_2[self.molecule3]), 1)

        # # NOTE
        # # processed_molecules_3 has skip smaller hits and all cpu threads
        # # skip small hits might fail with more than one thread
        # # because the result has already been computed due to parallelism
        # self.assertEqual(len(processed_molecules_3[self.molecule]), 1)
        # self.assertEqual(len(processed_molecules_3[self.molecule3]), 1)

        # self.assertEqual(
        #     [i.query_residue_count for i in processed_molecules_1[self.molecule]], [403, 403]
        # )

        self.assertEqual(
            [i.query_residue_count for i in processed_molecules_1[self.molecule2]],
            [511, 511],
        )

        self.assertEqual(
            [i.query_residue_count for i in processed_molecules_2[self.molecule3]],
            [495],
        )  # conservation cutoff applied

        processed_molecules_4 = self.template_matcher4.run(molecules=[self.molecule])
        self.assertEqual(
            len(processed_molecules_4[self.molecule]), 2
        )  # match-small-templates

    def test_Matcher_single_run(self):
        # we only have to check that filtering works,
        unfiltered_matches = self.unfiltered_matcher.run_single(molecule=self.molecule)

        filtered_matches = self.filtered_matcher.run_single(molecule=self.molecule)

        unfiltered_pdb_ids = []
        for match in unfiltered_matches:
            unfiltered_pdb_ids.append(match.hit.template().pdb_id)

        filtered_pdb_ids = []
        for match in filtered_matches:
            filtered_pdb_ids.append(match.hit.template().pdb_id)

        self.assertEqual(len(filtered_matches), 5)
        self.assertEqual(
            sorted(filtered_pdb_ids), ["1bf2", "1uh3", "1uh3", "1uh3", "2cxg"]
        )
        self.assertEqual(len(unfiltered_matches), 6)
        self.assertEqual(
            sorted(unfiltered_pdb_ids), ["1bf2", "1uh3", "1uh3", "1uh3", "2cxg", "2qy1"]
        )

        for match in unfiltered_matches:
            if match.index == 1:
                self.assertTrue(match.complete)
                self.assertEqual(match.hit.template().pdb_id, "2cxg")
            elif match.index == 2:
                self.assertTrue(match.complete)
                self.assertEqual(match.hit.template().pdb_id, "1uh3")
            elif match.index == 3:
                self.assertFalse(match.complete),
                self.assertEqual(match.hit.template().pdb_id, "2qy1")
            elif match.index == 4:
                self.assertFalse(match.complete)
                self.assertEqual(match.hit.template().pdb_id, "1bf2")
            elif match.index == 5:
                self.assertTrue(match.complete)
                self.assertEqual(match.hit.template().pdb_id, "1uh3")
            elif match.index == 6:
                self.assertTrue(match.complete)
                self.assertEqual(match.hit.template().pdb_id, "1uh3")
            else:
                ValueError("Shouldn't get here in the tests")

    def test_other_distances(self):
        # Test how pairwise distances outside the range of 0.7 to 2.0 are handeled
        unfiltered_matches = self.unfiltered_matcher_different_cutoffs.run_single(
            molecule=self.molecule
        )
        # here we expect returned matches all of which cannot be predicted on.
        # so they are none
        self.assertIsNone(unfiltered_matches[0].predicted_correct)

        filtered_matches = self.filtered_matcher_different_cutoffs.run_single(
            molecule=self.molecule
        )
        # filtering does not remove predicted_correct=None matches
        self.assertIsNone(filtered_matches[0].predicted_correct)
