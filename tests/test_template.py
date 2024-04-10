import unittest

from matcher import template
from pathlib import Path

matcher = Path(__package__).parent / 'matcher/'

class TestIntegration(unittest.TestCase):

    def test_load_templates(self):
        with self.assertRaises(FileNotFoundError):
            templates = list(template.load_templates("/some/bogus/folder"))

        ## End to End test # TODO
        # self.assertEqual(len(list(template.load_templates())), 7611) # check if all the expected templates were found
    
    def test_get_template_paths(self):
        data_path = Path(matcher / 'data/')
        # check if all the required annotation files are there - This will run the entire script
        self.assertEqual(template.get_template_paths(data_path, '.json'), [Path(data_path, 'pdb_sifts.json'), Path(data_path, 'MCSA_CATH_mapping.json'), Path(data_path, 'MCSA_EC_mapping.json')])


class TestVec3(unittest.TestCase):
    def setUp(self):
        self.vec1 = template.Vec3(1.0, 2.0, 3.0)
        self.vec2 = template.Vec3(-0.275, 2.8, 0.837)
        self.vec3 = template.Vec3(0, 0 ,0)

    def test_attributes(self):
        self.assertEqual(self.vec1.x, 1.0)
        self.assertEqual(self.vec1.y, 2.0)
        self.assertEqual(self.vec1.z, 3.0)
        self.assertEqual(self.vec2.x, -0.275)
        self.assertEqual(self.vec2.y, 2.8)
        self.assertEqual(self.vec2.z, 0.837)

class TestAtom(unittest.TestCase):
    def setUp(self):
        atomstring1 = 'ATOM      3  CD ZGLU A 156      20.734  32.611  20.206 E     1.86 '
        atomstring2 = 'ATOM      3  OE1ZGLU A 156      21.042  31.877  19.239 E     1.86 '
        atomstring3 = 'ATOM      3  OE2ZGLU A 156      21.258  32.515  21.340 E     1.86 '
        self.atom1 = template.Atom.loads(atomstring1)
        self.atom2 = template.Atom.loads(atomstring2)
        self.atom3 = template.Atom.loads(atomstring3)

    def test_attributes(self):
        self.assertEqual(self.atom3.x, 21.258)
        self.assertEqual(self.atom3.y, 32.515)
        self.assertEqual(self.atom3.z, 21.340)
        self.assertEqual(self.atom3.name, 'OE2')
        self.assertEqual(self.atom3.residue_name, 'GLU')
        self.assertEqual(self.atom3.allowed_residues, 'E')
        self.assertEqual(self.atom3.specificity, 3)
        self.assertEqual(self.atom3.backbone, False)
        self.assertEqual(self.atom3.residue_number, 156)
        self.assertEqual(self.atom3.chain_id, 'A')

        #self.assertEqual(self.atom1-self.atom2, template.Vec3(-0.308, 0.734, 0.967)) #TODO numerical error

class TestTemplate(unittest.TestCase):
    def setUp(self):
        with open(Path(matcher, 'jess_templates_20230210/6_residues/results/csa3d_0001/csa3d_0001.cluster_1_1_1.1b74_A147-AA180-AA70-AA178-AA8-AA7.template.pdb'), 'r') as f:
            self.template1 = template.Template.load(f.read())
    
    def test_attributes(self):
        self.assertEqual(self.template1.pdb_id, '1b74')
        self.assertEqual(self.template1.mcsa_id, 1)
        self.assertEqual(self.template1.cluster, '1_1_1')
        self.assertEqual(self.template1.uniprot_id, 'P56868')
        self.assertEqual(self.template1.organism, 'Aquifex pyrophilus')
        self.assertEqual(self.template1.organism_id, '2714')
        self.assertEqual(self.template1.resolution, 2.3)
        self.assertEqual(self.template1.experimental_method, 'X-ray diffraction')
        self.assertEqual(self.template1.ec, ['5.1.1.3']) # could be more info added through sifts
        self.assertEqual(self.template1.size, 6)
        self.assertEqual(self.template1.true_size, 6)
        self.assertEqual(self.template1.multimeric, True)
        self.assertEqual(self.template1.relative_order, [0])
        self.assertEqual(self.template1.cath, ['3.40.50.1860']) # could be more info added through sifts
        self.assertEqual(len(self.template1.residues), 6)

        # add test for non-multimeric template

class TestResidue(unittest.TestCase):
    def setUp(self):
        self.residue1 = TestTemplate.template1.residues[0]
        self.residue2 = TestTemplate.template1.residues[1]

    def test_attributes(self):
        self.assertEqual(self.residue1.residue_name, 'GLU')
        self.assertEqual(self.residue1.allowed_residues, 'E')
        self.assertEqual(self.residue1.specificity, 3)
        self.assertEqual(self.residue1.backbone, False)
        self.assertEqual(self.residue1.residue_number, 147)
        self.assertEqual(self.residue1.chain_id, 'A')
        self.assertEqual(self.residue1.orientation_vector, template.Vec3(x,y,z)) # Glu orientation is from C to midpoint between Os

