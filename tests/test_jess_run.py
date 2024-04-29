import unittest
from importlib.resources import files
from . import test_data

from pathlib import Path
import io

from matcher import template
import pyjess # type: ignore
from matcher.jess import filter_hits # type: ignore


matcher = Path(__package__).parent / 'matcher/'

class TestMatcher(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb'), 'r') as f:
            template_text1 = f.read()
            cls.template1 = template.Template.loads(template_text1, warn=False)
            pyjess_template1 = pyjess.Template.load(f, id=cls.template1.id)
            jess_1 = pyjess.Jess([pyjess_template1])

        with files(test_data).joinpath("1AMY.pdb").open() as f:
            molecule = pyjess.Molecule.load(f)

        query = jess_1.query(molecule=molecule, rmsd=2, distance_cutoff=2, max_dynamic_distance=2, max_candidates=10000)
        hits = list(query)
        print(hits)

        best_hits = filter_hits.filter_hits(hits)
        
        cls.match = Match(hit=best_hits[0], template=template1)

    def test_match(self):
        self.assertEqual(self.match.hit.molecule_id, '1amy')
        self.assertEqual(self.match.tempalte.id, self.template1.id)
        self.assertEqual(self.match.template.size, self.template1.size)
        self.assertEqual(self.match.template.true_size, self.template1.true_size)
        self.assertEqual(self.match.template.mcsa_id, self.template1.mcsa_id)
        self.assertEqual(self.match.template.uniprot_id, self.template1.uniprot_id)
        self.assertEqual(self.match.template.pdb_id, self.template1.pdb_id)
        self.assertEqual(self.match.template.ec, self.template1.ec)
        self.assertEqual(self.match.template.cath, self.template1.cath)
        self.assertEqual(self.match.template.multimeric, self.template1.multimeric)
        self.assertEqual(self.match.multimeric, False)
        self.assertAlmostEqual(self.match.hit.rmsd, 0.32093143)
        self.assertAlmostEqual(self.match.hit.log_evalue, -3.08424478)
        self.assertEqual(self.match.orientation, 0.188114)
        self.assertEqual(self.match.template_vector_list, [res.orientation_vector for res in self.template1.residues])
        self.assertEqual(self.match.match_vector_list, [Vec3()])
        self.assertEqual(self.match.preserved_resid_order, True)
        self.assertEqual(self.match.complete, True)
        self.assertEqual(self.match.matched_residues, ['GLU_A_204','ASP_A_87','ASP_A_179','HIS_A_288','ASP_A_289'])

    def test_match_dumps(self):
        dumps_string = 'query_id	template_id	template_size	template_true_size	template_mcsa_id	template_uniprot_id	template_pdb_id	template_ec	template_cath	template_multimeric	query_multimeric	rmsd	log_evalue	orientation	preserved_order	completeness	matched_residues\n1AMY	285_5_1uh3_1_1_1	5	5	285	Q60053	1uh3	3.2.1.10,3.2.1.13,3.2.1.135	3.20.20.80	False	False	0.32093143180639955	-3.084244780540347'
        self.assertEqual(self.match.dumps(), dumps_string)

    def test_match_dump2pdb(self):
        buffer = io.StringIO()
        self.match.dump2pdb(buffer)
        with files(test_data).joinpath("1AMY_matches.pdb").open() as f:
            self.assertEqual(buffer.getvalue(), f.getvalue())
