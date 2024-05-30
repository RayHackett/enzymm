import unittest

from matcher import template
from pathlib import Path
import math
import io
import warnings

from importlib.resources import files
from . import test_data

matcher = Path(__package__).parent / 'matcher/'

class TestIntegration(unittest.TestCase):

    def test_load_templates(self):
        with self.assertRaises(FileNotFoundError):
            templates = list(template.load_templates(Path("/some/bogus/folder")))
        
        with self.assertRaises(NotADirectoryError):
            templates = list(template.load_templates(Path(files(test_data).joinpath("bad_templates/pdb_id_none.pdb"))))

        with self.assertRaises(ValueError):
            templates = list(template.load_templates(Path(files(test_data))))

        # End to End test loading all supplied templates
        self.assertEqual(len(list(template.load_templates(warn=False))), 7611) # check if all the expected templates were found
    
    def test_get_template_paths(self):
        data_path = Path(matcher / 'data/')
        # check if all the required annotation files are there - This will run the entire script
        self.assertEqual(template._get_paths_by_extension(data_path, '.json'), [Path(data_path, 'pdb_sifts.json'), Path(data_path, 'MCSA_EC_mapping.json'), Path(data_path, 'MCSA_CATH_mapping.json')])

    def test__get_paths_by_extension(self):
        with self.assertRaises(FileNotFoundError):
            template._get_paths_by_extension(Path(files(test_data)), '.xyz')

class TestVec3(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vec1 = template.Vec3(1.0, 2.0, 3.0)
        cls.vec2 = template.Vec3(-0.275, 2.8, 0.837)
        cls.vec3 = template.Vec3(0, 0 ,0)
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
        self.assertEqual(self.vec1.normalize().x, 1/math.sqrt(14))
        self.assertEqual(self.vec1.normalize().y, 2/math.sqrt(14))
        self.assertEqual(self.vec1.normalize().z, 3/math.sqrt(14))
        self.assertEqual(self.vec3.normalize(), self.vec3)
    
    def test_add(self):
        result = self.vec1 + self.vec2
        expected = template.Vec3(x=0.725, y=4.8, z=3.837)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z)

        result = self.vec1 + 5
        expected = template.Vec3(x=6,y=7, z=8)
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
        self.assertAlmostEqual((self.vec1.normalize() @ self.vec2.normalize()), 0.713465)

        with self.assertRaises(TypeError):
            self.vec1 @ self.notavec

    def test_angle_to(self):
        self.assertAlmostEqual(self.vec1.angle_to(self.vec2), math.acos(0.713465))

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

class TestAnnotatedTemplate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            cls.template_text1 = f.read()
        with open(Path(matcher, 'jess_templates_20230210/4_residues/results/csa3d_0011/csa3d_0011.cluster_1_1_3.1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179.template.pdb'), 'r') as f:
            cls.template_text2 = f.read()
        with files(test_data).joinpath("bad_templates/empty_remark.pdb").open() as f:
            cls.template_empty_remark = f.read()
        with files(test_data).joinpath("bad_templates/hetatm.pdb").open() as f:
            cls.template_hetatm = f.read()
        with files(test_data).joinpath("bad_templates/malformed_residues_1.pdb").open() as f:
            cls.template_malformed_residue_1 = f.read()
        with files(test_data).joinpath("bad_templates/malformed_residues_2.pdb").open() as f:
            cls.template_malformed_residue_2 = f.read()
        with files(test_data).joinpath("bad_templates/malformed_residues_3.pdb").open() as f:
            cls.template_malformed_residue_3 = f.read()
        with files(test_data).joinpath("bad_templates/missing_cluster_annotation.pdb").open() as f:
            cls.template_missing_cluster_annotation = f.read()
        with files(test_data).joinpath("bad_templates/new_remark.pdb").open() as f:
            cls.template_new_remark = f.read()
        with files(test_data).joinpath("bad_templates/pdb_id_none.pdb").open() as f:
            cls.template_pdb_id_none = f.read()
        with files(test_data).joinpath("bad_templates/not_a_template.pdb").open() as f:
            cls.not_a_template = f.read()
        with files(test_data).joinpath("bad_templates/no_remark.pdb").open() as f:
            cls.template_no_remark = f.read()
    
    def test_good_loads(self):
        template1 = template.AnnotatedTemplate.annotated_loads(self.template_text1, str(0), warn=True)
        template2 = template.AnnotatedTemplate.annotated_loads(self.template_text2, str(1), warn=True)
        self.assertEqual(template1.id, '0')
        self.assertEqual(template1.pdb_id, '1b74')
        self.assertEqual(template1.template_id_string, '1b74_A147-AA180-AA70-AA178-AA8-AA7')
        self.assertEqual(template1.mcsa_id, 1)
        self.assertEqual(template1.cluster.id, 1)
        self.assertEqual(template1.cluster.member, 1)
        self.assertEqual(template1.cluster.size, 1)
        self.assertEqual(template1.uniprot_id, 'P56868')
        self.assertEqual(template1.organism, 'Aquifex pyrophilus')
        self.assertEqual(template1.organism_id, '2714')
        self.assertEqual(template1.resolution, 2.3)
        self.assertEqual(template1.experimental_method, 'X-ray diffraction')
        self.assertEqual(template1.ec, ['5.1.1.3']) # could be more info added through sifts
        self.assertEqual(template1.represented_sites, 2)
        self.assertEqual(template1.enzyme_discription, 'GLUTAMATE RACEMASE (E.C.5.1.1.3)')
        self.assertEqual(template1.effective_size, 6)
        self.assertEqual(template1.dimension, 6)
        self.assertEqual(template1.multimeric, True)
        self.assertEqual(template1.relative_order, [0])
        self.assertEqual(template1.cath, ['3.40.50.1860']) # could be more info added through sifts
        self.assertEqual(len(template1.residues), 6)

        self.assertEqual(template2.id, '1')
        self.assertEqual(template2.pdb_id, '1qum')
        self.assertEqual(template2.template_id_string, '1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179')
        self.assertEqual(template2.mcsa_id, 11)
        self.assertEqual(template2.cluster.id, 1)
        self.assertEqual(template2.cluster.member, 1)
        self.assertEqual(template2.cluster.size, 3)
        self.assertEqual(template2.uniprot_id, None)
        self.assertEqual(template2.organism, 'Escherichia coli')
        self.assertEqual(template2.organism_id, '562')
        self.assertEqual(template2.resolution, 1.55)
        self.assertEqual(template2.experimental_method, 'X-ray diffraction')
        self.assertEqual(template2.ec, ['3.1.21.2']) # could be more info added through sifts
        self.assertEqual(template2.represented_sites, 1)
        self.assertEqual(template2.enzyme_discription, 'ENDONUCLEASE IV (E.C.3.1.21.2)/DNA')
        self.assertEqual(template2.dimension, 4)
        self.assertEqual(template2.effective_size, 4)
        self.assertEqual(template2.multimeric, False)
        self.assertEqual(template2.relative_order, [2,1,4,3])
        self.assertEqual(template2.cath, ['3.20.20.150']) # could be more info added through sifts
        self.assertEqual(len(template2.residues), 4)

    def test_bad_loads(self):
        with self.assertWarns(Warning):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_empty_remark, str(0), warn=True)
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_hetatm, str(0), warn=True)
        with self.assertRaises(IndexError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_missing_cluster_annotation, str(0), warn=True)
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_malformed_residue_1, str(0), warn=True)
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_malformed_residue_2, str(0), warn=True)
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_malformed_residue_3, str(0), warn=True)

        # These do not raise Warnings or errors
        template1 = template.AnnotatedTemplate.annotated_loads(self.template_no_remark, str(0), warn=True)
        template1 = template.AnnotatedTemplate.annotated_loads(self.template_new_remark, str(0), warn=True)
        template1 = template.AnnotatedTemplate.annotated_loads(self.template_pdb_id_none, str(0), warn=True)
        template1 = template.AnnotatedTemplate.annotated_loads(self.not_a_template, str(0), warn=True)

    def test_dump(self):
        pass # TODO when atoms have dump method

    def test_dumps(self):
        pass # TODO when atoms have dump method

    def test_annotation_parsing(self):
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_pdb_id(tokens=[0, 1, 'abcde'], metadata={})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_uniprot_id(tokens=[0, 1, 'abcde'], metadata={})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_mcsa_id(tokens=[0, 1, 'abcde'], metadata={})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_cluster(tokens=[0, 1, 'abcde'], metadata={})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_cluster(tokens=[0, 1, 'ab_cd_e'], metadata={})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_resolution(tokens=[0, 1, 'abcde'], metadata={})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_ec(tokens=[0, 1, '1.1.1'], metadata={'ec': []})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_ec(tokens=[0, 1, '8.1.1.1'], metadata={'ec': []})
        with self.assertWarns(Warning):
            template.AnnotatedTemplate._parse_ec(tokens=[0, 1, '1.1.1.n1'], metadata={'ec': []})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_cath(tokens=[0, 1, '1.1.1.n1'], metadata={'cath': []})
        template.AnnotatedTemplate._parse_cath(tokens=[0, 1, '1.1.1.1'], metadata={'cath': []})
        template.AnnotatedTemplate._parse_cath(tokens=[0, 1, '1.1.1800.1'], metadata={'cath': []})
        with self.assertRaises(ValueError):
            template.AnnotatedTemplate._parse_represented_sites(tokens=[0, 1, 'abcde'], metadata={})

class TestResidue(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            template1 = template.AnnotatedTemplate.annotated_loads(f.read(), str(0), warn=False)
        cls.residue1 = template1.residues[0]
        cls.residue2 = template1.residues[1]

        with files(test_data).joinpath("bad_templates/unk_residue.pdb").open() as f:
            cls.template_unkown_residue_text = f.read()
        with files(test_data).joinpath("bad_templates/bad_atom_names_1.pdb").open() as f:
            cls.bad_atom_names_1 = f.read()
        with files(test_data).joinpath("bad_templates/bad_atom_names_2.pdb").open() as f:
            cls.bad_atom_names_2 = f.read()
        with files(test_data).joinpath("bad_templates/bad_atom_names_3.pdb").open() as f:
            cls.bad_atom_names_3 = f.read()

    def test_unkown_resids(self):
        with self.assertRaises(KeyError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.template_unkown_residue_text, str(0), warn=True)

    def test_bad_atoms_in_orientation(self):
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.bad_atom_names_1, str(0), warn=True)
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.bad_atom_names_2, str(0), warn=True)
        with self.assertRaises(ValueError):
            template1 = template.AnnotatedTemplate.annotated_loads(self.bad_atom_names_3, str(0), warn=True)

    def test_attributes(self):
        self.assertEqual(self.residue1.residue_name, 'GLU')
        self.assertEqual(self.residue1.allowed_residues, 'E')
        self.assertEqual(self.residue1.match_mode, 3)
        self.assertEqual(self.residue1.backbone, False)
        self.assertEqual(self.residue1.residue_number, 147)
        self.assertEqual(self.residue1.chain_id, 'A')
        result = self.residue1.orientation_vector # Glu orientation is from C to midpoint between Os
        expected = template.Vec3(-0.125,-0.347499999,0.45599999)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z) 
        self.assertEqual(self.residue1.orientation_vector_indices, (0,9)) # from O to midpoint

        result = self.residue2.orientation_vector # His orientation is from from CG to ND1
        expected = template.Vec3(1.07900000,-0.63800000,0.570000000)
        self.assertAlmostEqual(result.x, expected.x)
        self.assertAlmostEqual(result.y, expected.y)
        self.assertAlmostEqual(result.z, expected.z) 
        self.assertEqual(self.residue2.orientation_vector_indices, (0,1)) # from atom 0 to atom 1

class TestTemplate_Checking(unittest.TestCase):
    # TODO improve this function and write more tests!
    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            cls.template_text1 = f.read()
    def test_check_template(self):
        template1 = template.AnnotatedTemplate.annotated_loads(self.template_text1, str(0), warn=False)
        self.assertEqual(template.check_template(template=template1, warn=False), True)
        self.assertEqual(template.check_template(template=template1, warn=True), True)
