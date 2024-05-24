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

        ## End to End test # TODO
        # self.assertEqual(len(list(template.load_templates())), 7611) # check if all the expected templates were found
    
    def test_get_template_paths(self):
        data_path = Path(matcher / 'data/')
        # check if all the required annotation files are there - This will run the entire script
        self.assertEqual(template._get_paths_by_extension(data_path, '.json'), [Path(data_path, 'pdb_sifts.json'), Path(data_path, 'MCSA_EC_mapping.json'), Path(data_path, 'MCSA_CATH_mapping.json')])

class TestVec3(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vec1 = template.Vec3(1.0, 2.0, 3.0)
        cls.vec2 = template.Vec3(-0.275, 2.8, 0.837)
        cls.vec3 = template.Vec3(0, 0 ,0)

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

    def test_matmul(self):
        self.assertEqual(self.vec1 @ self.vec3, 0)
        self.assertAlmostEqual((self.vec1 @ self.vec2), 7.836)
        self.assertAlmostEqual((self.vec1.normalize() @ self.vec2.normalize()), 0.713465)

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

class TestTemplate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            cls.template_text1 = f.read()
        with open(Path(matcher, 'jess_templates_20230210/4_residues/results/csa3d_0011/csa3d_0011.cluster_1_1_3.1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179.template.pdb'), 'r') as f:
            cls.template_text2 = f.read()
    
    def test_loads(self):
        template1 = template.AnnotatedTemplate.loads(self.template_text1, warn=False)
        template2 = template.AnnotatedTemplate.loads(self.template_text2, warn=False)
        self.assertEqual(template1.pdb_id, '1b74')
        self.assertEqual(template1.id_residue_string, '1b74_A147-AA180-AA70-AA178-AA8-AA7')
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

        self.assertEqual(template2.pdb_id, '1qum')
        self.assertEqual(template2.id_residue_string, '1qum_D145-D109-D37-D72-D69-D229-D182-D231-D261-D216-D179')
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

class TestResidue(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            template1 = template.AnnotatedTemplate.loads(f.read(), warn=False)
        cls.residue1 = template1.residues[0]
        cls.residue2 = template1.residues[1]

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

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            cls.template_text1 = f.read()
    def test_check_template(self):
        template1 = template.AnnotatedTemplate.loads(self.template_text1, warn=False)
        template1_path = Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb')
        self.assertEqual(template.check_template(template_tuple=(template1, template1_path), warn=False), True)
        self.assertEqual(template.check_template(template_tuple=(template1, template1_path), warn=True), True)
