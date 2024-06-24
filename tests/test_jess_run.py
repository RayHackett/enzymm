import errno
import unittest
from importlib.resources import files
from . import test_data

from pathlib import Path
import io
import os

from matcher import template, jess_run
import pyjess

matcher = Path(__package__).parent / 'matcher/'

class TestMatch(unittest.TestCase):

    maxDiff = None

    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb'), 'r') as f:
            template_text1 = f.read()
            cls.template1 = template.AnnotatedTemplate.annotated_loads(template_text1, str(0), warn=False)
            pyjess_template1 = pyjess.Template.loads(template_text1, id=cls.template1.id)
            jess_1 = pyjess.Jess([pyjess_template1])

        with files(test_data).joinpath("1AMY.pdb").open() as f:
            molecule = pyjess.Molecule.load(f)

        query = jess_1.query(molecule, 2, 1.5, 1.5, max_candidates=10000, best_match=True)
        best_hits = list(query)
        
        cls.match1 = jess_run.Match(hit=best_hits[0], template=cls.template1)
        cls.match1.complete = True  # this is false because we did not run the completeness check here. we manually set it to true

        with open(Path(matcher, 'jess_templates_20230210/3_residues/results/csa3d_0415/csa3d_0415.cluster_1_1_2.1be0_A124-A175-A125-A289-A260.template.pdb'), 'r') as f:
            template_text2 = f.read()
            cls.template2 = template.AnnotatedTemplate.annotated_loads(template_text2, str(0), warn=False)
            pyjess_template2 = pyjess.Template.loads(template_text2, id=cls.template2.id)
            jess_2 = pyjess.Jess([pyjess_template2])

        query = jess_2.query(molecule, 2, 1, 1, max_candidates=10000, best_match=True)
        best_hits = list(query)
        
        cls.match2 = jess_run.Match(hit=best_hits[0], template=cls.template2)

    def test_match(self):
        self.assertEqual(self.match1.hit.molecule.id, '1AMY')
        self.assertEqual(self.match1.template.id, self.template1.id)
        self.assertEqual(self.match1.template.pdb_id, self.template1.pdb_id)
        self.assertEqual(self.match1.template.effective_size, self.template1.effective_size)
        self.assertEqual(self.match1.template.dimension, self.template1.dimension)
        self.assertEqual(self.match1.template.mcsa_id, self.template1.mcsa_id)
        self.assertEqual(self.match1.template.uniprot_id, self.template1.uniprot_id)
        self.assertEqual(self.match1.template.ec, self.template1.ec)
        self.assertEqual(self.match1.template.cath, self.template1.cath)
        self.assertEqual(self.match1.template.multimeric, self.template1.multimeric)
        self.assertEqual(self.match1.multimeric, False)
        self.assertEqual(self.match1.query_atom_count, 3339)
        self.assertEqual(self.match1.query_residue_count, 558)
        self.assertAlmostEqual(self.match1.hit.rmsd, 0.32093143)
        self.assertAlmostEqual(self.match1.hit.log_evalue, -3.08424478)
        self.assertAlmostEqual(self.match1.orientation, 0.15327054322)
        self.assertEqual(self.match1.template_vector_list, [res.orientation_vector for res in self.template1.residues])
        self.assertEqual(self.match1.match_vector_list, [template.Vec3(x=0.2290067979141952, y=-0.3853409610281773, z=0.377114677867322), template.Vec3(x=0.4249816660862038, y=-0.21966898402981627, z=-0.3540863184957992), template.Vec3(x=0.45459385444007694, y=-0.34869961601989985, z=0.10687378206512577), template.Vec3(x=-0.8733960645698886, y=0.2563504028143271, z=-0.9840695023070225), template.Vec3(x=-0.510183600042339, y=-0.1958417994791759, z=0.18963368325429997)])
        self.assertEqual(self.match1.preserved_resid_order, True)
        self.assertEqual(self.match1.complete, True)
        self.assertEqual(self.match1.matched_residues, [('GLU', 'A', '204'),('ASP', 'A', '87'),('ASP', 'A', '179'),('HIS', 'A', '288'),('ASP', 'A', '289')])

        self.assertEqual(self.match2.template.id, self.template2.id)
        self.assertEqual(self.match2.template.pdb_id, self.template2.pdb_id)
        self.assertEqual(self.match2.template.effective_size, self.template2.effective_size)
        self.assertEqual(self.match2.template.dimension, self.template2.dimension)
        self.assertEqual(self.match2.template.mcsa_id, self.template2.mcsa_id)
        self.assertEqual(self.match2.template.uniprot_id, self.template2.uniprot_id)
        self.assertEqual(self.match2.template.ec, self.template2.ec)
        self.assertEqual(self.match2.template.cath, self.template2.cath)
        self.assertEqual(self.match2.template.multimeric, self.template2.multimeric)
        self.assertEqual(self.match2.multimeric, False)
        self.assertAlmostEqual(self.match2.hit.rmsd, 1.7353479120)
        self.assertAlmostEqual(self.match2.hit.log_evalue, 1.8517964201)
        self.assertAlmostEqual(self.match2.orientation, 1.6503123465442575)
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
        dumps_string = 'query_id	template_pdb_id	template_pdb_chains	template_cluster_id	template_cluster_member	template_cluster_size	template_effective_size	template_dimension	template_mcsa_id	template_uniprot_id	template_ec	template_cath	template_multimeric	query_multimeric	query_atom_count	query_residue_count	rmsd	log_evalue	orientation	preserved_order	completeness	matched_residues\n1AMY	1uh3	A	1	1	1	5	5	285	Q60053	3.2.1.10,3.2.1.135	2.60.40.10,2.60.40.1180,3.20.20.80	False	False	3339	558	0.32093143180639955	-3.084244780540347	0.1532705432273403	True	True	GLU_A_204,ASP_A_87,ASP_A_179,HIS_A_288,ASP_A_289\n'
        self.assertEqual(self.match1.dumps(header=True), dumps_string)

    def test_match_dump2pdb(self):
        buffer = io.StringIO()
        self.match1.dump2pdb(buffer, transform=False, include_query=False)
        with files(test_data).joinpath("1AMY_matches_query.pdb").open() as f:
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

class Test_single_query_run(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_1.1uh3_A396-A262-A356-A471-A472.template.pdb'), 'r') as f:
            template_text1 = f.read()

        with open(Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0045/csa3d_0045.cluster_1_1_1.2cxg_A227-A229-A257-A327-A328.template.pdb'), 'r') as f:
            template_text2 = f.read()

        with open(Path(matcher, 'jess_templates_20230210/3_residues/results/csa3d_0914/csa3d_0914.cluster_4_2_2.3gvr_A41-A164-A37-A50-A86.template.pdb'), 'r') as f:
            template_text3 = f.read()

        with open(Path(matcher, 'jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_1_2.1uh3_A396-A262-A356-A471-A472.template.pdb'), 'r') as f:
            template_text4 = f.read()

        with open(Path(matcher, 'jess_templates_20230210/3_residues/results/csa3d_0285/csa3d_0285.cluster_1_2_2.1uh3_A396-A262-A356-A471-A472.template.pdb'), 'r') as f:
            template_text5 = f.read()

        template_list = []
        for index, txt in enumerate([template_text1, template_text2, template_text3, template_text4, template_text5]):
            template_list.append(template.AnnotatedTemplate.annotated_loads(txt, str(index), warn=False))
        cls.template_list = template_list

        with files(test_data).joinpath("bad_templates/no_remark.pdb").open() as f:
            template_no_remark = f.read()
        cls.unannotated_templates = [template.AnnotatedTemplate.annotated_loads(template_no_remark, str(index), warn=False)]

        cls.bad_template_list = [template.AnnotatedTemplate.annotated_loads(template_text1, str(0), warn=False), template.AnnotatedTemplate.annotated_loads(template_text2, str(0), warn=False)]
            
        with files(test_data).joinpath("1AMY.pdb").open() as f:
            cls.molecule = pyjess.Molecule.load(f)

    def test_single_query_run_function(self):
        # we only have to check that filtering works, completeness checking works and that templates are correctly mapped back again
        matches = jess_run.single_query_run(molecule=self.molecule, templates=self.template_list, rmsd_threshold=2, distance_cutoff=1.5, max_dynamic_distance=1.5)
        pdb_ids = []
        for match in matches:
            pdb_ids.append(match.template.pdb_id)
        pdb_ids.sort()

        self.assertEqual(len(matches), 5)
        self.assertEqual(pdb_ids, ['1uh3', '1uh3', '1uh3', '2cxg', '3gvr'])

        for match in matches:
            if match.template.id == '0':
                self.assertEqual(match.complete, True)
                self.assertEqual(match.template.pdb_id, '1uh3')
            elif match.template.id == '1':
                self.assertEqual(match.complete, True)
                self.assertEqual(match.template.pdb_id, '2cxg')
            elif match.template.id == '2':
                self.assertEqual(match.complete, False)
                self.assertEqual(match.template.pdb_id, '3gvr')
            elif match.template.id == '3':
                self.assertEqual(match.complete, True)
                self.assertEqual(match.template.pdb_id, '1uh3')
            elif match.template.id == '4':
                self.assertEqual(match.complete, True)
                self.assertEqual(match.template.pdb_id, '1uh3')

        matches2 = jess_run.single_query_run(molecule=self.molecule, templates=self.unannotated_templates, rmsd_threshold=2, distance_cutoff=1.5, max_dynamic_distance=1.5)
        for match in matches2:
            self.assertEqual(match.complete, True)

    def test_nonunique_template_id_error(self):
        with self.assertRaises(KeyError):
            jess_run.single_query_run(molecule=self.molecule, templates=self.bad_template_list, rmsd_threshold=2, distance_cutoff=1.5, max_dynamic_distance=1.5)

class TestMatcher(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        res5_templates = list(template.load_templates(template_dir=Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0285/'), warn=False, verbose=False))
        res4_templates = list(template.load_templates(template_dir=Path(matcher, 'jess_templates_20230210/4_residues/results/csa3d_0285'), warn=False, verbose=False))
        res3_templates = list(template.load_templates(template_dir=Path(matcher, 'jess_templates_20230210/3_residues/results/csa3d_0344'), warn=False, verbose=False))

        selected_templates = res5_templates + res4_templates
        with_smaller_templates = selected_templates + res3_templates
        cls.template_matcher = jess_run.Matcher(templates=selected_templates, cpus=2, warn=True)
        cls.template_matcher2 = jess_run.Matcher(templates=selected_templates, skip_smaller_hits=True)
        with cls.assertWarns(cls, Warning):
            cls.template_matcher3 = jess_run.Matcher(templates=with_smaller_templates, match_small_templates=True, warn=True, cpus=-1)

        with files(test_data).joinpath("1AMY.pdb").open() as f:
            cls.molecule = pyjess.Molecule.load(f)
        with files(test_data).joinpath("AF-P0DUB6-F1-model_v4.pdb").open() as f:
            cls.molecule2 = pyjess.Molecule.load(f)
            cls.molecule3 = cls.molecule2.conserved(80)

    def test_init(self):
        jess_params = {
            3: {'rmsd': 2, 'distance': 0.9, 'max_dynamic_distance': 0.9},
            4: {'rmsd': 2, 'distance': 1.7, 'max_dynamic_distance': 1.7},
            5: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0},
            6: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0},
            7: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0},
            8: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0}}

        self.assertEqual(self.template_matcher.jess_params, jess_params)
        self.assertEqual(self.template_matcher.cpus, 2)
        self.assertEqual(self.template_matcher.template_effective_sizes, [5, 4])
        self.assertEqual(self.template_matcher2.cpus, (os.cpu_count() or 1)-2)
        self.assertEqual(self.template_matcher3.cpus, (os.cpu_count() or 1))

    def test_Matcher_run(self):
        processed_molecules_1 = self.template_matcher.run(molecules=[self.molecule, self.molecule2])
        self.assertEqual(list(processed_molecules_1.keys())[0], self.molecule)
        self.assertEqual(list(processed_molecules_1.keys())[1], self.molecule2)
        self.assertEqual(len(processed_molecules_1[self.molecule]), 2)
        self.assertEqual(len(processed_molecules_1[self.molecule2]), 2)
        self.assertEqual(processed_molecules_1[self.molecule2][0].query_residue_count, 511)

        processed_molecules_2 = self.template_matcher2.run(molecules=[self.molecule, self.molecule3])
        self.assertEqual(len(processed_molecules_2[self.molecule]), 1) # skip smaller hits
        self.assertEqual(processed_molecules_2[self.molecule3][0].query_residue_count, 494) # conservation cutoff applied

        processed_molecules_3 = self.template_matcher3.run(molecules=[self.molecule])
        self.assertEqual(len(processed_molecules_3[self.molecule]), 3) # match-small-templates

# TODO test conservation

class Test_main(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        molecule_path = str(Path(files(test_data).joinpath("1AMY.pdb")).resolve())

        list_path = str(Path(files(test_data).joinpath("input_list.txt")).resolve())
        list_path2 = str(Path(files(test_data).joinpath("input_dir.txt")).resolve())

        selected_template_dir = str(Path(matcher, 'jess_templates_20230210/5_residues/results/csa3d_0285/').resolve())

        cls.arguments_normal = ['-i', molecule_path, '-t', selected_template_dir, '-o', 'temp_result.tsv']
        cls.arguments_list = ['-l', list_path, '-t', selected_template_dir, '-o', 'temp_result.tsv', '-j', '5', '1.1', '1.5', '-c', '2']
        cls.arguments_both = ['-i', molecule_path, '-l', list_path, '-t', selected_template_dir, '-o', 'temp_result.tsv']

        cls.bad_argument_1 = ['-l', molecule_path, '-o', 'temp_result.tsv'] # not a file in list file
        cls.bad_argument_2 = ['-l', list_path2, '-o', 'temp_result.tsv'] # passing a dir in list file

        # TODO either remove the temp_result.tsv file or avoid writing it? I cant see how that oculd work though

    def test_default_main(self):
        self.assertEqual(jess_run.main(self.arguments_normal), 0)
        self.assertEqual(jess_run.main(self.arguments_list), 0)
        self.assertEqual(jess_run.main(self.arguments_both), 0)

        retcode = jess_run.main(self.bad_argument_1, stderr=io.StringIO())
        self.assertEqual(retcode, errno.ENOENT)
        retcode2 = jess_run.main(self.bad_argument_2, stderr=io.StringIO())
        self.assertEqual(retcode2, errno.EISDIR)