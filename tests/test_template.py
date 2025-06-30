import unittest
import matcher
from matcher import template
from pathlib import Path
import math

import pyjess

import importlib.resources
from importlib.resources import files
from . import test_data


class TestIntegration(unittest.TestCase):

    def test_load_templates(self):
        with self.assertRaises(NotADirectoryError):
            list(template.load_templates(Path("/some/bogus/folder")))

        with self.assertRaises(NotADirectoryError):
            list(
                template.load_templates(
                    files(test_data).joinpath("bad_templates/pdb_id_none.pdb")
                )
            )

        with self.assertRaises(ValueError):
            list(template.load_templates(Path(files(test_data))))

        # End to End test loading all supplied templates
        self.assertEqual(
            len(list(template.load_templates(template_dir=None))), 7607
        )  # check if all the expected templates were found

    def test_get_template_paths(self):
        data_path = importlib.resources.files(matcher).joinpath("data/")
        found_paths = sorted(template._get_paths_by_extension(data_path, ".json"))
        expected_paths = sorted(
            [
                Path(data_path, "pdb_sifts.json"),
                Path(data_path, "MCSA_EC_mapping.json"),
                Path(data_path, "MCSA_CATH_mapping.json"),
                Path(data_path, "logistic_regression_models.json"),
                Path(data_path, "catalytic_residue_homologs_information.json"),
            ]
        )
        self.assertEqual(found_paths, expected_paths)

    def test__get_paths_by_extension(self):
        with self.assertRaises(FileNotFoundError):
            template._get_paths_by_extension(Path(files(test_data)), ".xyz")


class TestVec3(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vec1 = template.Vec3(1.0, 2.0, 3.0)
        cls.vec2 = template.Vec3(-0.275, 2.8, 0.837)
        cls.vec3 = template.Vec3(0, 0, 0)
        cls.vec4 = template.Vec3(1.00000000001, 0, 0)
        cls.vec5 = template.Vec3(-1, 0, 0)
        cls.vec6 = template.Vec3(
            -0.43667809853452577, 0.6652199071133453, 0.6056357927338702
        )
        cls.vec7 = template.Vec3(
            -0.4366780985345203, 0.6652199071133459, 0.6056357927338736
        )
        cls.notavec = (0, 3, 5)

    def test_attributes(self):
        self.assertEqual(self.vec1.x, 1.0)
        self.assertEqual(self.vec1.y, 2.0)
        self.assertEqual(self.vec1.z, 3.0)
        self.assertEqual(self.vec2.x, -0.275)
        self.assertEqual(self.vec2.y, 2.8)
        self.assertEqual(self.vec2.z, 0.837)

    def test_norm(self):
        self.assertEqual(self.vec1.norm, math.sqrt(14))
        self.assertEqual(self.vec3.norm, math.sqrt(0))

    def test_normalize(self):
        self.assertEqual(self.vec1.normalize().x, 1 / math.sqrt(14))
        self.assertEqual(self.vec1.normalize().y, 2 / math.sqrt(14))
        self.assertEqual(self.vec1.normalize().z, 3 / math.sqrt(14))
        self.assertEqual(self.vec3.normalize(), self.vec3)

    def test_add(self):
        result = self.vec1 + self.vec2
        expected = template.Vec3(x=0.725, y=4.8, z=3.837)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        result = self.vec1 + 5
        expected = template.Vec3(x=6, y=7, z=8)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        self.assertEqual(self.vec1 + self.vec3, self.vec1)

        with self.assertRaises(TypeError):
            self.vec1 + self.notavec

    def test_sub(self):
        result = self.vec1 - self.vec2
        expected = template.Vec3(x=1.275, y=-0.8, z=2.163)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        result = self.vec1 - 5
        expected = template.Vec3(x=-4, y=-3, z=-2)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        self.assertEqual(self.vec1 - self.vec3, self.vec1)

        with self.assertRaises(TypeError):
            self.vec1 - self.notavec

    def test_truediv(self):
        result = self.vec1 / self.vec2
        expected = template.Vec3(x=-3.6363636, y=0.7142857, z=3.58422939)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        result = self.vec1 / 2
        expected = template.Vec3(x=0.5, y=1, z=1.5)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        with self.assertRaises(ZeroDivisionError):
            self.vec1 / self.vec3

        with self.assertRaises(TypeError):
            self.vec1 / self.notavec

    def test_matmul(self):
        self.assertEqual(self.vec1 @ self.vec3, 0)
        self.assertAlmostEqual((self.vec1 @ self.vec2), 7.836)
        self.assertAlmostEqual(
            (self.vec1.normalize() @ self.vec2.normalize()), 0.713465
        )

        with self.assertRaises(TypeError):
            self.vec1 @ self.notavec

    def test_angle_to(self):
        self.assertAlmostEqual(self.vec1.angle_to(self.vec2), math.acos(0.713465))
        self.assertAlmostEqual(self.vec4.angle_to(self.vec4), 0)
        self.assertAlmostEqual(self.vec4.angle_to(self.vec5), math.pi)
        self.assertAlmostEqual(self.vec6.angle_to(self.vec7), 0)


class TestCluster(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cluster1 = template.Cluster(3, 1, 2)

    def test_cluster_initialization(self):
        with self.assertRaises(ValueError):
            template.Cluster(1, 2, 1)

    def test_attributes(self):
        self.assertEqual(self.cluster1.id, 3)
        self.assertEqual(self.cluster1.size, 2)
        self.assertEqual(self.cluster1.member, 1)


class TestTemplate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(
            importlib.resources.files(matcher).joinpath(
                "jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb",
            ),
            "r",
        ) as f:
            cls.template_text1 = f.read()
        with open(
            importlib.resources.files(matcher).joinpath(
                "jess_templates_20230210/4_residues/results/csa3d_0011/csa3d_0011.cluster_1_1_3.1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179.template.pdb",
            ),
            "r",
        ) as f:
            cls.template_text2 = f.read()

        with files(test_data).joinpath("bad_templates/hetatm.pdb").open() as f:
            cls.template_hetatm = f.read()
        with files(test_data).joinpath(
            "bad_templates/malformed_residues_1.pdb"
        ).open() as f:
            cls.template_malformed_residue_1 = f.read()
        with files(test_data).joinpath(
            "bad_templates/malformed_residues_2.pdb"
        ).open() as f:
            cls.template_malformed_residue_2 = f.read()
        with files(test_data).joinpath(
            "bad_templates/malformed_residues_3.pdb"
        ).open() as f:
            cls.template_malformed_residue_3 = f.read()
        with files(test_data).joinpath(
            "bad_templates/missing_cluster_annotation.pdb"
        ).open() as f:
            cls.template_missing_cluster_annotation = f.read()
        with files(test_data).joinpath("bad_templates/new_remark.pdb").open() as f:
            cls.template_new_remark = f.read()
        with files(test_data).joinpath("bad_templates/pdb_id_none.pdb").open() as f:
            cls.template_pdb_id_none = f.read()
        with files(test_data).joinpath("bad_templates/not_a_template.pdb").open() as f:
            cls.not_a_template = f.read()
        with files(test_data).joinpath("bad_templates/no_remark.pdb").open() as f:
            cls.template_no_remark = f.read()

        magic_atom = pyjess.TemplateAtom(
            chain_id="A",
            residue_number=51,
            residue_names=["ANY"],
            atom_names=["C"],
            distance_weight=1.5,
            match_mode=1,
            x=1.0,
            y=0,
            z=-99.5,
        )

        magic_atom2 = pyjess.TemplateAtom(
            chain_id="A",
            residue_number=51,
            residue_names=["ANY"],
            atom_names=["A"],
            distance_weight=1.5,
            match_mode=1,
            x=1.5,
            y=-1,
            z=15,
        )

        magic_atom3 = pyjess.TemplateAtom(
            chain_id="A",
            residue_number=51,
            residue_names=["ANY"],
            atom_names=["O"],
            distance_weight=1.5,
            match_mode=1,
            x=0,
            y=-2.3,
            z=-2.233,
        )

        cls.artifical_template1 = template.Template(
            id="some_name",
            # pdb_id="Magic_ID",
            template_id_string="A completely fake template",
            cluster=template.Cluster(id=1, member=2, size=2),
            # uniprot_id="NotAUniProtID",
            organism=None,
            organism_id="a string??",
            resolution=99,
            experimental_method="pure creativity",
            enzyme_discription="brain chemistry goes brr",
            cath=["1.1.1.1"],
            ec=["9.9.9.9"],
            residues=template.Residue.construct_residues_from_atoms(
                atoms=[magic_atom, magic_atom2, magic_atom3]
            ),
        )

        cls.artifical_template2 = template.Template(
            id="some_name",
            # pdb_id="Magic_ID_and_different_di",
            template_id_string="A completely fake template",
            cluster=template.Cluster(id=1, member=2, size=2),
            # uniprot_id="A different uniprot",
            organism=None,
            organism_id="a string??",
            resolution=99,
            experimental_method="pure creativity",
            enzyme_discription="brain chemistry goes brr",
            cath=["1.1.1.1"],
            ec=["9.9.9.9"],
            residues=template.Residue.construct_residues_from_atoms(
                atoms=[
                    magic_atom,
                    magic_atom2,
                    magic_atom3,
                    magic_atom,
                    magic_atom2,
                    magic_atom3,
                ]
            ),
        )

        cls.minimal_template = template.Template(
            residues=template.Residue.construct_residues_from_atoms(
                atoms=[magic_atom, magic_atom2, magic_atom3]
            )
        )

        cls.template1 = template.Template.loads(cls.template_text1, warn=False)

        cls.template1_with_id = template.Template.loads(
            cls.template_text1, id="another_id", warn=False
        )
        cls.template2 = template.Template.loads(
            cls.template_text2, id="hello_world", warn=True
        )

    def test_good_loads(self):

        self.assertEqual(self.template1.pdb_id, "1b74")
        self.assertEqual(
            self.template1.template_id_string, "1b74_A147-AA180-AA70-AA178-AA8-AA7"
        )
        self.assertEqual(self.template1.mcsa_id, 1)
        self.assertEqual(self.template1.cluster.id, 1)
        self.assertEqual(self.template1.cluster.member, 1)
        self.assertEqual(self.template1.cluster.size, 1)
        self.assertEqual(self.template1.uniprot_id, "P56868")
        self.assertEqual(self.template1.organism, "Aquifex pyrophilus")
        self.assertEqual(self.template1.organism_id, "2714")
        self.assertEqual(self.template1.resolution, 2.3)
        self.assertEqual(self.template1.experimental_method, "X-ray diffraction")
        self.assertEqual(
            self.template1.ec, tuple(["5.1.1.3"])
        )  # could be more info added through sifts
        self.assertEqual(self.template1.represented_sites, 2)
        self.assertEqual(
            self.template1.enzyme_discription, "GLUTAMATE RACEMASE (E.C.5.1.1.3)"
        )
        self.assertEqual(self.template1.effective_size, 6)
        self.assertEqual(self.template1.dimension, 6)
        self.assertEqual(self.template1.multimeric, True)
        self.assertEqual(self.template1.relative_order, [0])
        self.assertEqual(
            self.template1.cath, tuple(["3.40.50.1860"])
        )  # could be more info added through sifts
        self.assertEqual(len(self.template1.residues), 6)

        self.assertEqual(self.template2.pdb_id, "1qum")
        self.assertEqual(self.template2.id, "hello_world")
        self.assertEqual(
            self.template2.template_id_string,
            "1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179",
        )
        self.assertEqual(self.template2.mcsa_id, 11)
        self.assertEqual(self.template2.cluster.id, 1)
        self.assertEqual(self.template2.cluster.member, 1)
        self.assertEqual(self.template2.cluster.size, 3)
        self.assertEqual(self.template2.uniprot_id, None)
        self.assertEqual(self.template2.organism, "Escherichia coli")
        self.assertEqual(self.template2.organism_id, "562")
        self.assertEqual(self.template2.resolution, 1.55)
        self.assertEqual(self.template2.experimental_method, "X-ray diffraction")
        self.assertEqual(
            self.template2.ec, tuple(["3.1.21.2"])
        )  # could be more info added through sifts
        self.assertEqual(self.template2.represented_sites, 1)
        self.assertEqual(
            self.template2.enzyme_discription, "ENDONUCLEASE IV (E.C.3.1.21.2)/DNA"
        )
        self.assertEqual(self.template2.dimension, 4)
        self.assertEqual(self.template2.effective_size, 4)
        self.assertEqual(self.template2.multimeric, False)
        self.assertEqual(self.template2.relative_order, [2, 1, 4, 3])
        self.assertEqual(
            self.template2.cath, tuple(["3.20.20.150"])
        )  # could be more info added through sifts
        self.assertEqual(len(self.template2.residues), 4)

    def test_copy(self):
        template1_copy = self.template1.copy()
        self.assertEqual(self.template1, self.template1)
        self.assertEqual(
            list(self.template1), list(template1_copy)
        )  # testing atoms copy
        self.assertEqual(self.template1.id, template1_copy.id)
        self.assertEqual(self.template1.pdb_id, template1_copy.pdb_id)
        self.assertEqual(self.template1.mcsa_id, template1_copy.mcsa_id)
        self.assertEqual(
            self.template1.template_id_string, template1_copy.template_id_string
        )
        self.assertEqual(self.template1.cluster, template1_copy.cluster)
        self.assertEqual(self.template1.effective_size, template1_copy.effective_size)
        self.assertEqual(self.template1.uniprot_id, template1_copy.uniprot_id)
        self.assertEqual(self.template1.organism, template1_copy.organism)
        self.assertEqual(self.template1.organism_id, template1_copy.organism_id)
        self.assertEqual(self.template1.resolution, template1_copy.resolution)
        self.assertEqual(
            self.template1.experimental_method, template1_copy.experimental_method
        )
        self.assertEqual(
            self.template1.enzyme_discription, template1_copy.enzyme_discription
        )
        self.assertEqual(
            self.template1.represented_sites, template1_copy.represented_sites
        )
        self.assertEqual(self.template1.ec, template1_copy.ec)
        self.assertEqual(self.template1.cath, template1_copy.cath)
        self.assertEqual(self.template1, template1_copy)
        # Test with lots of None's
        self.assertEqual(self.minimal_template, self.minimal_template.copy())

    def test_artifical_template_copy(self):
        self.assertEqual(self.artifical_template1, self.artifical_template1.copy())

    def test_template_non_equality(self):
        # Test for template equality
        self.assertEqual(hash(self.template1), hash(self.template1.copy()))
        # Changed only the internal id
        self.assertNotEqual(self.template1, self.template1_with_id)
        # completely different templates
        self.assertNotEqual(self.template1, self.template2)

    def test_artifical_template_non_equality(self):
        # changed only the uniprot id
        self.assertNotEqual(
            hash(self.artifical_template1), hash(self.artifical_template2)
        )
        self.assertFalse(self.artifical_template1 == self.artifical_template2)
        self.assertTrue(self.artifical_template1 != self.artifical_template2)
        self.assertNotEqual(self.artifical_template1, self.artifical_template2)
        self.assertNotEqual(
            hash(self.artifical_template1), hash(self.artifical_template2)
        )

    def test_bad_loads(self):
        with self.assertRaises(ValueError):
            template.Template.loads(self.template_hetatm, warn=True)
        with self.assertRaises(IndexError):
            template.Template.loads(self.template_missing_cluster_annotation, warn=True)
        with self.assertRaises(ValueError):
            template.Template.loads(self.template_malformed_residue_1, warn=True)
        with self.assertRaises(ValueError):
            template.Template.loads(self.template_malformed_residue_2, warn=True)
        with self.assertRaises(ValueError):
            template.Template.loads(self.template_malformed_residue_3, warn=True)
        with self.assertRaises(ValueError):
            template.Template.loads(self.not_a_template, warn=True)

        # These do not raise Warnings or errors
        template.Template.loads(self.template_no_remark, warn=True)
        template.Template.loads(self.template_new_remark, warn=True)
        template.Template.loads(self.template_pdb_id_none, warn=True)

    # def test_dumps(self): # TODO
    #     # tests both dump and dumps methods
    #     template1 = template.Template.loads(
    #         self.template_text1, warn=False
    #     )

    #     template1_dumpstring = template1.dumps()

    #     self.assertEqual(self.template_text1, template1_dumpstring)

    def test_annotation_parsing(self):
        with self.assertRaises(ValueError):
            template.Template._parse_pdb_id(tokens=[0, 1, "abcde"], metadata={})
        with self.assertRaises(ValueError):
            template.Template._parse_uniprot_id(tokens=[0, 1, "abcde"], metadata={})
        with self.assertRaises(ValueError):
            template.Template._parse_mcsa_id(tokens=[0, 1, "abcde"], metadata={})
        with self.assertRaises(ValueError):
            template.Template._parse_cluster(tokens=[0, 1, "abcde"], metadata={})
        with self.assertRaises(ValueError):
            template.Template._parse_cluster(tokens=[0, 1, "ab_cd_e"], metadata={})
        with self.assertRaises(ValueError):
            template.Template._parse_resolution(tokens=[0, 1, "abcde"], metadata={})
        with self.assertRaises(ValueError):
            template.Template._parse_ec(tokens=[0, 1, "1.1.1"], metadata={"ec": []})
        with self.assertRaises(ValueError):
            template.Template._parse_ec(tokens=[0, 1, "8.1.1.1"], metadata={"ec": []})
        with self.assertWarns(Warning):
            template.Template._parse_ec(tokens=[0, 1, "1.1.1.n1"], metadata={"ec": []})
        with self.assertRaises(ValueError):
            template.Template._parse_cath(
                tokens=[0, 1, "1.1.1.n1"], metadata={"cath": []}
            )
        template.Template._parse_cath(tokens=[0, 1, "1.1.1.1"], metadata={"cath": []})
        template.Template._parse_cath(
            tokens=[0, 1, "1.1.1800.1"], metadata={"cath": []}
        )
        with self.assertRaises(ValueError):
            template.Template._parse_represented_sites(
                tokens=[0, 1, "abcde"], metadata={}
            )


class TestAnnotatedTemplate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(
            importlib.resources.files(matcher).joinpath(
                "jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb",
            ),
            "r",
        ) as f:
            cls.template_text1 = f.read()
        with open(
            importlib.resources.files(matcher).joinpath(
                "jess_templates_20230210/4_residues/results/csa3d_0011/csa3d_0011.cluster_1_1_3.1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179.template.pdb",
            ),
            "r",
        ) as f:
            cls.template_text2 = f.read()

        with files(test_data).joinpath("bad_templates/hetatm.pdb").open() as f:
            cls.template_hetatm = f.read()
        with files(test_data).joinpath(
            "bad_templates/malformed_residues_1.pdb"
        ).open() as f:
            cls.template_malformed_residue_1 = f.read()
        with files(test_data).joinpath(
            "bad_templates/malformed_residues_2.pdb"
        ).open() as f:
            cls.template_malformed_residue_2 = f.read()
        with files(test_data).joinpath(
            "bad_templates/malformed_residues_3.pdb"
        ).open() as f:
            cls.template_malformed_residue_3 = f.read()
        with files(test_data).joinpath(
            "bad_templates/missing_cluster_annotation.pdb"
        ).open() as f:
            cls.template_missing_cluster_annotation = f.read()
        with files(test_data).joinpath("bad_templates/new_remark.pdb").open() as f:
            cls.template_new_remark = f.read()
        with files(test_data).joinpath("bad_templates/pdb_id_none.pdb").open() as f:
            cls.template_pdb_id_none = f.read()
        with files(test_data).joinpath("bad_templates/not_a_template.pdb").open() as f:
            cls.not_a_template = f.read()
        with files(test_data).joinpath("bad_templates/no_remark.pdb").open() as f:
            cls.template_no_remark = f.read()

        cls.template1 = template.AnnotatedTemplate.loads(cls.template_text1, warn=False)

        cls.template1_with_id = template.AnnotatedTemplate.loads(
            cls.template_text1, id="another_id", warn=False
        )
        cls.template2 = template.AnnotatedTemplate.loads(
            cls.template_text2, id="hello_world", warn=True
        )

    # TODO test the anntoations added via annotated load fn
    def test_good_loads(self):
        pass

    def test_copy(self):
        template1_copy = self.template1.copy()
        self.assertEqual(self.template1, self.template1)
        self.assertEqual(
            list(self.template1), list(template1_copy)
        )  # testing atoms copy
        self.assertEqual(self.template1, template1_copy)

    def test_template_non_equality(self):
        # Test for template equality
        self.assertEqual(hash(self.template1), hash(self.template1.copy()))
        # Changed only the internal id
        self.assertNotEqual(self.template1, self.template1_with_id)
        # completely different templates
        self.assertNotEqual(self.template1, self.template2)

    def test_bad_loads(self):
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate.loads(self.template_hetatm, warn=True)
        with self.assertRaises(IndexError):
            template.AnnotatedTemplate.loads(
                self.template_missing_cluster_annotation, warn=True
            )
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate.loads(
                self.template_malformed_residue_1, warn=True
            )
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate.loads(
                self.template_malformed_residue_2, warn=True
            )
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate.loads(
                self.template_malformed_residue_3, warn=True
            )
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate.loads(self.not_a_template, warn=True)

        # These do not raise Warnings or errors
        template.AnnotatedTemplate.loads(self.template_no_remark, warn=True)
        template.AnnotatedTemplate.loads(self.template_new_remark, warn=True)
        template.AnnotatedTemplate.loads(self.template_pdb_id_none, warn=True)


class TestResidue(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(
            importlib.resources.files(matcher).joinpath(
                "jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb",
            ),
            "r",
        ) as f:
            template1 = template.Template.loads(
                f.read(),
            )
        cls.residue1 = template1.residues[0]
        cls.residue2 = template1.residues[1]

        with files(test_data).joinpath("bad_templates/unk_residue.pdb").open() as f:
            cls.template_unkown_residue_text = f.read()
        with files(test_data).joinpath(
            "bad_templates/bad_atom_names_1.pdb"
        ).open() as f:
            cls.bad_atom_names_1 = f.read()
        with files(test_data).joinpath(
            "bad_templates/bad_atom_names_2.pdb"
        ).open() as f:
            cls.bad_atom_names_2 = f.read()
        with files(test_data).joinpath(
            "bad_templates/bad_atom_names_3.pdb"
        ).open() as f:
            cls.bad_atom_names_3 = f.read()

    def test_unkown_resids(self):
        with self.assertRaises(KeyError):
            template.Template.loads(self.template_unkown_residue_text, warn=True)

    def test_bad_atoms_in_orientation(self):
        with self.assertRaises(ValueError):
            template.Template.loads(self.bad_atom_names_1, warn=True)
        with self.assertRaises(ValueError):
            template.Template.loads(self.bad_atom_names_2, warn=True)
        with self.assertRaises(ValueError):
            template.Template.loads(self.bad_atom_names_3, warn=True)

    def test_attributes(self):
        self.assertEqual(self.residue1.residue_name, "GLU")
        self.assertEqual(self.residue1.allowed_residues, "E")
        self.assertEqual(self.residue1.match_mode, 3)
        self.assertEqual(self.residue1.backbone, False)
        self.assertEqual(self.residue1.residue_number, 147)
        self.assertEqual(self.residue1.chain_id, "A")
        result = (
            self.residue1.orientation_vector
        )  # Glu orientation is from C to midpoint between Os
        expected = template.Vec3(-0.125, -0.347499999, 0.45599999)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)
        self.assertEqual(
            self.residue1.orientation_vector_indices, (0, 9)
        )  # from O to midpoint

        result = (
            self.residue2.orientation_vector
        )  # His orientation is from from CG to ND1
        expected = template.Vec3(1.07900000, -0.63800000, 0.570000000)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)
        self.assertEqual(
            self.residue2.orientation_vector_indices, (0, 1)
        )  # from atom 0 to atom 1


# TODO test annotated Residue class


class TestTemplate_Checking(unittest.TestCase):
    # TODO improve this function and write more tests!
    @classmethod
    def setUpClass(cls):
        with open(
            importlib.resources.files(matcher).joinpath(
                "jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb",
            ),
            "r",
        ) as f:
            cls.template_text1 = f.read()

    def test_check_template(self):
        template1 = template.Template.loads(
            self.template_text1,
        )
        self.assertEqual(template.check_template(template=template1, warn=False), True)
        self.assertEqual(template.check_template(template=template1, warn=True), True)
