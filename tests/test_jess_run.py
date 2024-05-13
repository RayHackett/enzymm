import unittest
from importlib.resources import files
from . import test_data

from pathlib import Path
import io

from matcher import template, jess_run
import pyjess # type: ignore
from matcher.jess import filter_hits # type: ignore


matcher = Path(__package__).parent / 'matcher/'

class TestMatcher(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb'), 'r') as f:
            template_text1 = f.read()
            cls.template1 = template.Template.loads(template_text1, warn=False)
            pyjess_template1 = pyjess.Template.loads(template_text1, id=cls.template1.id)
            jess_1 = pyjess.Jess([pyjess_template1])

        with files(test_data).joinpath("1AMY.pdb").open() as f:
            molecule = pyjess.Molecule.load(f)

        query = jess_1.query(molecule, 2, 1.5, 1.5, max_candidates=10000)
        hits = list(query)
        best_hits = filter_hits.filter_hits(hits)
        
        cls.match1 = jess_run.Match(hit=best_hits[0], template=cls.template1)
        cls.match1.complete = True  # this is false because we did not run the completeness check here. we manually set it to true

        with open(Path(matcher, 'jess_templates_20230210/3_residues/results/csa3d_0415/csa3d_0415.cluster_1_1_2.1be0_A124-A175-A125-A289-A260.template.pdb'), 'r') as f:
            template_text2 = f.read()
            cls.template2 = template.Template.loads(template_text2, warn=False)
            pyjess_template2 = pyjess.Template.loads(template_text2, id=cls.template2.id)
            jess_2 = pyjess.Jess([pyjess_template2])

        query = jess_2.query(molecule, 2, 1, 1, max_candidates=10000)
        hits = list(query)
        best_hits = filter_hits.filter_hits(hits)
        
        cls.match2 = jess_run.Match(hit=best_hits[0], template=cls.template2)

    def test_match(self):
        self.assertEqual(self.match1.hit.molecule.id, '1AMY')
        self.assertEqual(self.match1.template.id, self.template1.id)
        self.assertEqual(self.match1.template.size, self.template1.size)
        self.assertEqual(self.match1.template.true_size, self.template1.true_size)
        self.assertEqual(self.match1.template.mcsa_id, self.template1.mcsa_id)
        self.assertEqual(self.match1.template.uniprot_id, self.template1.uniprot_id)
        self.assertEqual(self.match1.template.pdb_id, self.template1.pdb_id)
        self.assertEqual(self.match1.template.ec, self.template1.ec)
        self.assertEqual(self.match1.template.cath, self.template1.cath)
        self.assertEqual(self.match1.template.multimeric, self.template1.multimeric)
        self.assertEqual(self.match1.multimeric, False)
        self.assertEqual(self.match1.query_atom_count, 3339)
        self.assertEqual(self.match1.query_residue_count, 558)
        self.assertAlmostEqual(self.match1.hit.rmsd, 0.32093143)
        self.assertAlmostEqual(self.match1.hit.log_evalue, -3.08424478)
        self.assertAlmostEqual(self.match1.orientation, 0.18811427657)
        self.assertEqual(self.match1.template_vector_list, [res.orientation_vector for res in self.template1.residues])
        self.assertEqual(self.match1.match_vector_list, [template.Vec3(x=0.2290067979141952, y=-0.3853409610281773, z=0.377114677867322), template.Vec3(x=0.4249816660862038, y=-0.21966898402981627, z=-0.3540863184957992), template.Vec3(x=0.45459385444007694, y=-0.34869961601989985, z=0.10687378206512577), template.Vec3(x=-0.8733960645698886, y=0.0, z=-0.9840695023070225), template.Vec3(x=-0.510183600042339, y=-0.1958417994791759, z=0.18963368325429997)])
        self.assertEqual(self.match1.preserved_resid_order, True)
        self.assertEqual(self.match1.complete, True)
        self.assertEqual(self.match1.matched_residues, [('GLU', 'A', '204'),('ASP', 'A', '87'),('ASP', 'A', '179'),('HIS', 'A', '288'),('ASP', 'A', '289')])

        self.assertEqual(self.match2.template.id, self.template2.id)
        self.assertEqual(self.match2.template.size, self.template2.size)
        self.assertEqual(self.match2.template.true_size, self.template2.true_size)
        self.assertEqual(self.match2.template.mcsa_id, self.template2.mcsa_id)
        self.assertEqual(self.match2.template.uniprot_id, self.template2.uniprot_id)
        self.assertEqual(self.match2.template.pdb_id, self.template2.pdb_id)
        self.assertEqual(self.match2.template.ec, self.template2.ec)
        self.assertEqual(self.match2.template.cath, self.template2.cath)
        self.assertEqual(self.match2.template.multimeric, self.template2.multimeric)
        self.assertEqual(self.match2.multimeric, False)
        self.assertAlmostEqual(self.match2.hit.rmsd, 1.7353479120)
        self.assertAlmostEqual(self.match2.hit.log_evalue, 1.8517964201)
        self.assertAlmostEqual(self.match2.orientation, 1.7475836210)
        self.assertEqual(self.match2.template_vector_list, [res.orientation_vector for res in self.template2.residues])
        self.assertEqual(self.match2.preserved_resid_order, False)
        self.assertEqual(self.match2.complete, False)
        self.assertEqual(self.match2.matched_residues, [('TRP', 'A', '38'),('HIS', 'A', '288'),('ASP', 'A', '289')])

    def test_match_dump(self):
        buffer = io.StringIO()
        self.match1.dump(buffer, header=True)
        with files(test_data).joinpath("results.tsv").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

    def test_match_dumps(self):
        dumps_string = 'query_id	template_id	template_size	template_true_size	template_mcsa_id	template_uniprot_id	template_pdb_id	template_ec	template_cath	template_multimeric	query_multimeric	rmsd	log_evalue	orientation	preserved_order	completeness	matched_residues\n1AMY	285_5_1uh3_1_1_1	5	5	285	Q60053	1uh3	3.2.1.13,3.2.1.10,3.2.1.135	3.20.20.80	False	False	0.32093143180639955	-3.084244780540347	0.18811427657680269	True	True	GLU_A_204,ASP_A_87,ASP_A_179,HIS_A_288,ASP_A_289\n'
        self.assertEqual(self.match1.dumps(header=True), dumps_string)

    def test_match_dump2pdb(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=False, include_query=False) # TODO include checks for either options too
        with files(test_data).joinpath("1AMY_matches.pdb").open() as f:
            self.assertEqual(buffer.getvalue(), f.read())

# TODO test single_query_run and completeness check. This will be refactored anyway I assume

# TODO test matcher_run ?? this has no output to test.