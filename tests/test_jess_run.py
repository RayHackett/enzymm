import io
import os
import unittest
import importlib.resources
from importlib.resources import files
from . import test_data

import pyjess

import enzymm
from enzymm import template, jess_run


class TestMatch(unittest.TestCase):

    maxDiff = None

    @classmethod
    def setUpClass(cls):
        with open(
            importlib.resources.files(enzymm).joinpath(
                "jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb",
            ),
            "r",
        ) as f:
            template_text1 = f.read()
            cls.template1 = template.AnnotatedTemplate.loads(template_text1, warn=False)
            jess_1 = pyjess.Jess([cls.template1])

        with files(test_data).joinpath("1AMY.pdb").open() as f:
            molecule = pyjess.Molecule.load(f)

        query = jess_1.query(
            molecule, 2, 1.5, 1.5, max_candidates=10000, best_match=True
        )
        best_hits = list(query)

        cls.match1 = jess_run.Match(
            hit=best_hits[0], pairwise_distance=1.5, complete=True, index=0
        )

        with open(
            importlib.resources.files(enzymm).joinpath(
                "jess_templates_20230210/3_residues/results/csa3d_0415/csa3d_0415.cluster_1_1_2.1be0_A124-A175-A125-A289-A260.template.pdb",
            ),
            "r",
        ) as f:
            template_text2 = f.read()
            cls.template2 = template.AnnotatedTemplate.loads(template_text2, warn=False)
            jess_2 = pyjess.Jess([cls.template2])

        query = jess_2.query(molecule, 2, 1, 1, max_candidates=10000, best_match=True)
        best_hits = list(query)

        cls.match2 = jess_run.Match(hit=best_hits[0])

    def test_match(self):
        self.assertEqual(self.match1.hit.molecule().id, "1AMY")
        self.assertEqual(self.match1.index, 0)
        self.assertEqual(self.match1.hit.template.pdb_id, self.template1.pdb_id)
        self.assertEqual(
            self.match1.hit.template.effective_size, self.template1.effective_size
        )
        self.assertEqual(self.match1.hit.template.dimension, self.template1.dimension)
        self.assertEqual(self.match1.hit.template.mcsa_id, self.template1.mcsa_id)
        self.assertEqual(self.match1.hit.template.uniprot_id, self.template1.uniprot_id)
        self.assertEqual(self.match1.hit.template.ec, self.template1.ec)
        self.assertEqual(self.match1.hit.template.cath, self.template1.cath)
        self.assertEqual(self.match1.hit.template.multimeric, self.template1.multimeric)
        self.assertEqual(self.match1.multimeric, False)
        self.assertEqual(self.match1.query_atom_count, 3339)
        self.assertEqual(self.match1.query_residue_count, 403)
        self.assertAlmostEqual(self.match1.hit.rmsd, 0.32093143)
        self.assertAlmostEqual(self.match1.hit.log_evalue, -3.08424478)
        self.assertAlmostEqual(self.match1.orientation, 0.15327054322)
        self.assertEqual(
            self.match1.template_vector_list,
            [res.orientation_vector for res in self.template1.residues],
        )
        self.assertEqual(
            self.match1.match_vector_list,
            [
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
            ],
        )
        self.assertEqual(self.match1.preserved_resid_order, True)
        self.assertEqual(self.match1.complete, True)
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
        self.assertEqual(self.match2.hit.template.pdb_id, self.template2.pdb_id)
        self.assertEqual(
            self.match2.hit.template.effective_size, self.template2.effective_size
        )
        self.assertEqual(self.match2.hit.template.dimension, self.template2.dimension)
        self.assertEqual(self.match2.hit.template.mcsa_id, self.template2.mcsa_id)
        self.assertEqual(self.match2.hit.template.uniprot_id, self.template2.uniprot_id)
        self.assertEqual(self.match2.hit.template.ec, self.template2.ec)
        self.assertEqual(self.match2.hit.template.cath, self.template2.cath)
        self.assertEqual(self.match2.hit.template.multimeric, self.template2.multimeric)
        self.assertEqual(self.match2.multimeric, False)
        self.assertAlmostEqual(self.match2.hit.rmsd, 1.7353479120)
        self.assertAlmostEqual(self.match2.hit.log_evalue, 1.8517964201)
        self.assertAlmostEqual(self.match2.orientation, 1.6503123465442575)
        self.assertEqual(
            self.match2.template_vector_list,
            [res.orientation_vector for res in self.template2.residues],
        )
        self.assertEqual(self.match2.preserved_resid_order, False)
        self.assertEqual(self.match2.complete, False)
        self.assertEqual(
            self.match2.matched_residues,
            [("TRP", "A", "38"), ("HIS", "A", "288"), ("ASP", "A", "289")],
        )

    def test_match_dump(self):
        buffer = io.StringIO()
        self.match1.dump(buffer, header=True)
        with files(test_data).joinpath("results.tsv").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_match_dumps(self):
        dumps_string = "query_id	pairwise_distance	match_index	template_pdb_id	template_pdb_chains	template_cluster_id	template_cluster_member	template_cluster_size	template_effective_size	template_dimension	template_mcsa_id	template_uniprot_id	template_ec	template_cath	template_multimeric	query_multimeric	query_atom_count	query_residue_count	rmsd	log_evalue	orientation	preserved_order	completeness	predicted_correct	matched_residues	number_of_mutated_residues	number_of_side_chain_residues_(template,reference)	number_of_metal_ligands_(template,reference)	number_of_ptm_residues_(template, reference)	total_reference_residues\n1AMY	1.5	0	1uh3	A	1	1	1	5	5	285	Q60053	3.2.1.10,3.2.1.135	2.60.40.10,2.60.40.1180,3.20.20.80	False	False	3339	403	0.32093143180639955	-3.084244780540347	0.1532705432273403	True	True	True	GLU_A_204,ASP_A_87,ASP_A_179,HIS_A_288,ASP_A_289	0	5,5	0,0	0,0	5\n"
        self.assertEqual(self.match1.dumps(header=True), dumps_string)

    def test_match_dump2pdb(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=False, include_query=False)
        with files(test_data).joinpath("1AMY_matches_no_query.pdb").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_match_dump2pdb_with_query(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=False, include_query=True)
        with files(test_data).joinpath("1AMY_matches_query_included.pdb").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_match_dump2pdb_transformed(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=True, include_query=False)
        with files(test_data).joinpath("1AMY_matches_template.pdb").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())


class TestMatcher(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        template_files = [
            "jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/5_residues/results/csa3d_0045/csa3d_0045.cluster_1_1_1.2cxg_A227-A229-A257-A327-A328.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_2.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_2_2.1uh3_A396-A262-A356-A471-A472.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0421/csa3d_0421.cluster_1_1_3.1bf2_A229-A232-A230-A259-A375-A435-A510-A128.template.pdb",
            "jess_templates_20230210/3_residues/results/csa3d_0896/csa3d_0896.cluster_2_1_2.2qy1_A135-A179-A137-A230-A175.template.pdb",
        ]

        template_texts = []
        for path in template_files:
            with open(importlib.resources.files(enzymm).joinpath(path), "r") as f:
                template_texts.append(f.read())

        cls.template_list = []
        for txt in template_texts:
            cls.template_list.append(template.AnnotatedTemplate.loads(txt))

        res5_templates = list(
            template.load_templates(
                template_dir=importlib.resources.files(enzymm).joinpath(
                    "jess_templates_20230210/5_residues/results/csa3d_0285/"
                ),
            )
        )
        res4_templates = list(
            template.load_templates(
                template_dir=importlib.resources.files(enzymm).joinpath(
                    "jess_templates_20230210/4_residues/results/csa3d_0285"
                ),
            )
        )
        res3_templates = list(
            template.load_templates(
                template_dir=importlib.resources.files(enzymm).joinpath(
                    "jess_templates_20230210/3_residues/results/csa3d_0344"
                ),
            )
        )

        res5_and_res4_templates = res5_templates + res4_templates
        with_smaller_templates = res5_and_res4_templates + res3_templates

        cls.template_matcher = jess_run.Matcher(
            templates=res5_and_res4_templates, cpus=2, warn=True
        )
        cls.template_matcher2 = jess_run.Matcher(
            templates=res5_and_res4_templates, skip_smaller_hits=True
        )
        with cls.assertWarns(cls, Warning):
            cls.template_matcher3 = jess_run.Matcher(
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

        with files(test_data).joinpath("1AMY.pdb").open() as f:
            cls.molecule = pyjess.Molecule.load(f)

        with files(test_data).joinpath("AF-P0DUB6-F1-model_v4.pdb").open() as f:
            cls.molecule2 = pyjess.Molecule.load(f)
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
        self.assertEqual(self.template_matcher2.cpus, (os.cpu_count() or 1))
        self.assertEqual(self.template_matcher3.cpus, (os.cpu_count() or 1) - 1)

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

        self.assertEqual(list(processed_molecules_1.keys())[0], self.molecule)
        self.assertEqual(list(processed_molecules_1.keys())[1], self.molecule2)
        self.assertEqual(len(processed_molecules_1[self.molecule]), 2)
        self.assertEqual(len(processed_molecules_1[self.molecule2]), 2)
        #  processed_molecules_2 has skip smaller hits
        self.assertEqual(len(processed_molecules_2[self.molecule]), 1)
        self.assertEqual(len(processed_molecules_2[self.molecule3]), 1)

        # self.assertEqual(
        #     [i.query_residue_count for i in processed_molecules_1[self.molecule]], [403, 403]
        # )

        self.assertEqual(
            [i.query_residue_count for i in processed_molecules_1[self.molecule2]],
            [511, 511],
        )

        self.assertEqual(
            [i.query_residue_count for i in processed_molecules_2[self.molecule3]],
            [494],
        )  # conservation cutoff applied

        processed_molecules_3 = self.template_matcher3.run(molecules=[self.molecule])
        self.assertEqual(
            len(processed_molecules_3[self.molecule]), 3
        )  # match-small-templates

    def test_Matcher_single_run(self):
        # we only have to check that filtering works,
        unfiltered_matches = self.unfiltered_matcher.run_single(molecule=self.molecule)

        filtered_matches = self.filtered_matcher.run_single(molecule=self.molecule)

        unfiltered_pdb_ids = []
        for match in unfiltered_matches:
            unfiltered_pdb_ids.append(match.hit.template.pdb_id)

        filtered_pdb_ids = []
        for match in filtered_matches:
            filtered_pdb_ids.append(match.hit.template.pdb_id)

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
                self.assertEqual(match.complete, True)
                self.assertEqual(match.hit.template.pdb_id, "1uh3")
            elif match.index == 2:
                self.assertEqual(match.complete, True)
                self.assertEqual(match.hit.template.pdb_id, "2cxg")
            elif match.index == 3:
                self.assertEqual(match.complete, False),
                self.assertEqual(match.hit.template.pdb_id, "2qy1")
            elif match.index == 4:
                self.assertEqual(match.complete, True)
                self.assertEqual(match.hit.template.pdb_id, "1uh3")
            elif match.index == 5:
                self.assertEqual(match.complete, True)
                self.assertEqual(match.hit.template.pdb_id, "1uh3")
            elif match.index == 6:
                self.assertEqual(match.complete, False)
                self.assertEqual(match.hit.template.pdb_id, "1bf2")
            else:
                ValueError("Shouldn't get here in the tests")
