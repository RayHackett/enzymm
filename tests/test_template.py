import unittest

from matcher import template
from pathlib import Path
import math
import io

matcher = Path(__package__).parent / 'matcher/'

# how to test warnings -> turn them to errors and test those
'''
with assertRaises(UserWarning):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # change to error and then check for assertRaises some error
        run_whatever_test_code()
'''

class TestIntegration(unittest.TestCase):

    def test_load_templates(self):
        with self.assertRaises(FileNotFoundError):
            templates = list(template.load_templates(Path("/some/bogus/folder")))

        ## End to End test # TODO
        # self.assertEqual(len(list(template.load_templates())), 7611) # check if all the expected templates were found
    
    def test_get_template_paths(self):
        data_path = Path(matcher / 'data/')
        # check if all the required annotation files are there - This will run the entire script
        self.assertEqual(template.get_template_paths(data_path, '.json'), [Path(data_path, 'pdb_sifts.json'), Path(data_path, 'MCSA_EC_mapping.json'), Path(data_path, 'MCSA_CATH_mapping.json')])


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
        with self.assertRaises(ZeroDivisionError):
            self.vec3.normalize()

    def test_matmul(self):
        self.assertEqual(self.vec1 @ self.vec3, 0)
        self.assertAlmostEqual((self.vec1 @ self.vec2), 7.836)
        self.assertAlmostEqual((self.vec1.normalize() @ self.vec2.normalize()), 0.713465)

    def test_angle_to(self):
        self.assertAlmostEqual(self.vec1.angle_to(self.vec2), math.acos(0.713465))

class TestAtom(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.atomstring1 = 'ATOM      3  CD ZGLU A 156      20.734  32.611  20.206 E     1.86 '
        cls.atomstring2 = 'ATOM      3  OE1ZGLU A 156      21.042  31.877  19.239 E     1.86 '
        cls.atomstring3 = 'ATOM      3  OE2ZGLU A 156      21.258  32.515  21.340 E     1.86 '

    def test_attributes(self):
        atom1 = template.Atom.loads(self.atomstring1)
        atom2 = template.Atom.loads(self.atomstring2)
        atom3 = template.Atom.loads(self.atomstring3)
        self.assertEqual(atom3.x, 21.258)
        self.assertEqual(atom3.y, 32.515)
        self.assertEqual(atom3.z, 21.340)
        self.assertEqual(atom3.specificity, 3)
        self.assertEqual(atom3.name, 'OE2')
        self.assertEqual(atom3.altloc, 'Z')
        self.assertEqual(atom3.residue_name, 'GLU')
        self.assertEqual(atom3.allowed_residues, 'E')
        self.assertEqual(atom3.specificity, 3)
        self.assertEqual(atom3.backbone, False)
        self.assertEqual(atom3.residue_number, 156)
        self.assertEqual(atom3.chain_id, 'A')
        self.assertEqual(atom3.flexibility, 1.86)
        
    def test_diff__(self):
        atom1 = template.Atom(x=111, y=222, z=333, name='bla', altloc='Z', residue_name='GLU', allowed_residues='E', specificity=3, backbone=False, residue_number=1, chain_id='A', flexibility=1.86)
        atom2 = template.Atom(x=111, y=222, z=333, name='bla', altloc=None, residue_name='GLU', allowed_residues='E', specificity=3, backbone=False, residue_number=1, chain_id='A', flexibility=None)
        tmp = atom1-atom2
        self.assertAlmostEqual(tmp.x, 0, places=3)
        self.assertAlmostEqual(tmp.y, 0, places=3)
        self.assertAlmostEqual(tmp.z, 0, places=3)


class TestTemplate(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            cls.template_text = f.read()
    
    def test_loads(self):
        template1 = template.Template.loads(self.template_text)
        self.assertEqual(template1.pdb_id, '1b74')
        self.assertEqual(template1.mcsa_id, 1)
        self.assertEqual(template1.cluster_id, 1)
        self.assertEqual(template1.cluster_member, 1)
        self.assertEqual(template1.cluster_size, 1)
        self.assertEqual(template1.uniprot_id, 'P56868')
        self.assertEqual(template1.organism, 'Aquifex pyrophilus')
        self.assertEqual(template1.organism_id, 2714)
        self.assertEqual(template1.resolution, 2.3)
        self.assertEqual(template1.experimental_method, 'X-ray diffraction')
        self.assertEqual(template1.ec, {'5.1.1.3'}) # could be more info added through sifts
        self.assertEqual(template1.represented_sites, 2)
        self.assertEqual(template1.enzyme_discription, 'GLUTAMATE RACEMASE (E.C.5.1.1.3)')
        self.assertEqual(template1.size, 6)
        self.assertEqual(template1.true_size, 6)
        self.assertEqual(template1.multimeric, True)
        self.assertEqual(template1.relative_order, [0])
        self.assertEqual(template1.cath, {'3.40.50.1860'}) # could be more info added through sifts
        self.assertEqual(len(template1.residues), 6)

        # add test for non-multimeric template

class TestResidue(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        atom1 = template.Atom(x=100, y=100, z=100, name='CD', altloc='Z', residue_name='GLU', allowed_residues='E', specificity=3, backbone=False, residue_number=1, chain_id='A', flexibility=0)
        atom2 = template.Atom(x=150, y=100, z=0, name='bla', altloc='Z', residue_name='GLU', allowed_residues='E', specificity=3, backbone=False, residue_number=1, chain_id='A', flexibility=0)
        atom3 = template.Atom(x=50, y=100, z=0, name='bla', altloc='Z', residue_name='GLU', allowed_residues='E', specificity=3, backbone=False, residue_number=1, chain_id='A', flexibility=0)
        atom4 = template.Atom(x=0, y=0, z=0, name='CA', altloc='Z', residue_name='VAL', allowed_residues='AV', specificity=3, backbone=False, residue_number=2, chain_id='A', flexibility=0)
        atom5 = template.Atom(x=100, y=0, z=0, name='CB', altloc='Z', residue_name='VAL', allowed_residues='AV', specificity=3, backbone=False, residue_number=2, chain_id='A', flexibility=0)
        atom6 = template.Atom(x=0, y=100, z=0, name='bla', altloc='Z', residue_name='VAL', allowed_residues='AV', specificity=3, backbone=False, residue_number=2, chain_id='A', flexibility=0)

        cls.residue1 = template.Residue((atom1, atom2, atom3))
        cls.residue2 = template.Residue((atom4, atom5, atom6))

    def test_attributes(self):
        self.assertEqual(self.residue1.residue_name, 'GLU')
        self.assertEqual(self.residue1.allowed_residues, 'E')
        self.assertEqual(self.residue1.specificity, 3)
        self.assertEqual(self.residue1.backbone, False)
        self.assertEqual(self.residue1.residue_number, 1)
        self.assertEqual(self.residue1.chain_id, 'A')
        self.assertEqual(self.residue1.orientation_vector, template.Vec3(0,0,-100)) # Glu orientation is from C to midpoint between Os
        self.assertEqual(self.residue2.orientation_vector, template.Vec3(100,0,0)) # Val orientation is from from CA to CB

