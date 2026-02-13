EnzyMM Output
=============

EnzyMM produces a `TSV` table as output.
Optionally `PDB` structures with matched residues can be written too.

🖹 TSV Table
^^^^^^^^^^^

Each row shows data for a match between one of our catalytic templates and the query structure.
The table contains the following columns:

- **query_id**: `str` The id of the query - either user provided or derived from the header section of the structure.
- **pairwise_distance**: `float` The pairwise distance threshold at which this match was found in Ångstrom.
- **match_index**: `int` A running index for each match to a given query.
- **template_pdb_id**: `str` The PDB identifier to the experimental structure from which the template was derived.
- **template_pdb_chain**: `str` The PDB chain or chains from which the template was derived. Refers to the automatically assigned chain identifier in the biological assembly.
- **template_cluster_id**: `int` The id of the conformational cluster to which the template structure belongs.
- **template_cluster_member**: `int` The member index within the cluster. A member might be a partial catalytic site.
- **template_cluster_size**: `int` The total number of members within the cluster.
- **template_effective_size**: `int` The number of specific residues defining the template. Usually equal to the number of side chain interacting residues with specific match codes.
- **template_dimension**: `int` The toal number of residues, including unspecific ones, in the template.
- **template_mcsa_id**: `int` The entry number in the M-CSA to which the template refers. The M-CSA entry can help provide more information about the enzymatic mechanism in question and to trace annotations of the template back to sources in scientific literature.
- **template_uniprot_id**: `str` The UniProt identifier to the protein from which the template was derived. Each template comes from a single protein but some represent homo-multimers.
- **template_ec**: `list` of `str` Enzyme Commission numbers associated with the template which categorize the enzyme function(s) of the template.
- **template_cath**: `list` of `str` CATH identifiers to the Protein Structure Classification database describing domains of the template structure.
- **template_multimeric**: `bool` Wether the template contains multiple chains.
- **query_multimeric**: `bool` Wether residues from multiple chains in the query were matched.
- **query_atom_count**: `int` The number of atoms in the query model.
- **query_residue_count**: `int` The number of residues in the query model.
- **rmsd**: `float` Atom-wise Root-mean-square distance in Ångstrom between atoms matched between the template and the query structure. This metric shows how well the template superposes with the query structure. This metric contributes to filtering 3- and 4-residue matches. Results above 2 Å are never returned.
- **log_evalue**: `float` Statistical measure based on RMSD and template size which should be used with caution. Anything less than -4 should be a very good hit and -3 is OK. (E is the expected number of hits at random).
- **orientation**: `float` The mean of pairwise orientation angles in radians between corresponding residues between the template and the query structure. While related to superposition between template and query it is sensitive to changes to important chemical angles determining electrostatic interactions.
- **preserved_order**: `bool` Whether the relative order of residues in the protein sequence of the template and query is identical.
- **completeness**: `bool` Wether all members of the same cluster also matched the query structure. `True` if there is only one member.
- **predicted_correct**: `bool` Wether the match was predicted to be correct.
- **matched_residues**: `str` The matched residues in the query. The format is ['3-letter-code']_['chain-identifier']_['residue-number'].
- **number_of_mutated_residues**: `int` The number of mutated residues in the template.
- **number_of_side_chain_residues_(template,reference)**: `tuple` (`int`, `int`) The number of residues in the template which interact through their side chain and (all) the total number of residues including unspecific residues interacting through their main chain.
- **number_of_metal_ligands_(template,reference)**: `tuple` (`int`, `int`) The number of residues which contribute to metal binding or coordination in the template and the reference. A template composed of mostly metal binding residues is likely less predictive of catalytic function but might indicate a metal binding site.
- **number_of_ptm_residues_(template, reference)**: `tuple` (`int`, `int`) The number of posttranslationally modified residues in the template and the reference.
- **total_reference_residues**: `int` The number of catalytic residues annotated in the structure the template was derived from. Since a template might only represent a partial catalytic site, this shows how many residues the total site might be composed of.


PDB Structures
^^^^^^^^^^^^^^

Optionally you can write `PDB` files for the matches,
if you supply the `--pdbs` flag pointing to a directory.
What `PDB` files will be written depends on if the `--transform` flag is set.

By default (the `--transform` flag is not set), one `PDB` file will be written
per query structure with the matches in the coordinate system of the query molecule.

Alternatively if the `--transform` flag is set, one `PDB` file will be written per
template pdb structure with the matches in the coordinate system of the respective template.
This can make superposing multiple matches of different queries with the same template easier.

You can further configure the output in each `PDB` file with the following two flags:
- `--include-query`: will include the query molecule(s) in the `PDB` file
- `--include-template`: will include the template(s) in the `PDB` file

Each match is then added as a subsequent model. By default, matches are written in the
query reference frame so that many matches to the same query can we viewed together.

A single match in `PDB` format will look like this:

.. code::

    REMARK True MATCH 1AMY 0
    REMARK TEMPLATE_PDB 1uh3_A
    REMARK TEMPLATE CLUSTER 1_1_1
    REMARK TEMPLATE RESIDUES 1uh3_A396-A262-A356-A471-A472
    REMARK MOLECULE_ID 1AMY
    REMARK MATCH INDEX 0
    REMARK QUERY COORDINATE FRAME
    ATOM   1602  CD  GLU A 204       4.241  63.910  32.378  1.00  7.85           C 
    ATOM   1603  OE1 GLU A 204       3.145  64.416  32.618  1.00  6.76           O 
    ATOM   1604  OE2 GLU A 204       5.295  64.536  32.437  1.00 11.92           O 
    ATOM    673  CG  ASP A  87       0.258  65.949  25.285  1.00  6.24           C 
    ATOM    675  OD2 ASP A  87       1.453  66.078  25.016  1.00  2.00           O 
    ATOM    674  OD1 ASP A  87      -0.145  65.786  26.442  1.00  7.47           O 
    ATOM   1409  CG  ASP A 179       0.050  67.050  30.742  1.00 11.90           C 
    ATOM   1411  OD2 ASP A 179      -0.640  67.800  31.432  1.00  9.48           O 
    ATOM   1410  OD1 ASP A 179       1.266  67.144  30.660  1.00 16.59           O 

.. note::
    Multiple matches are seperated by `MDL` and `ENDMDL` lines!

⚠️ Limitations and Caveats
^^^^^^^^^^^^^^^^^^^^^^^^^^

Be critical and apply common sense when interpreting **EnzyMM** results.
Think of template matching as a fast and deterministic approach to identify any
structural pattern which satisfies the constraints of the template.
Since each template encodes different constraints, interpreting results from many
matches is not straightforward. Especially with smaller (and thus less specific)
templates, false positives are a real issue. We try to mitigate this issue by
applying stringent filtering based on RMSD and residue orientation. Nonetheless,
interpret results with caution. Our testing shows a precision in identifying
catalytic residues of roughly 80-95% depending on template size.

Given that a template approach is intrinsically limited by the coverage of the
template library — or in the case of **EnzyMM** the coverage of the M-CSA —
recall outside of these enzyme families might be limited.

Consider also that **EnzyMM** makes no prediction of EC number. It simply shows the EC
numbers associated with the template. For smaller templates which specify for
example metal binding sites or common motifs such as a protein transfer chain
or an oxyanion hole, EC prediction might be unreliable. EC numbers particularly
at the 4th level become substrate/product specific rather than categorizing
enzymatic mechanisms. Often templates are only descriptive up to the 3rd EC level.

Keep in mind that **EnzyMM** does not include checks for pocket volume or
substrate/solvent accessibility. Matches to random non-catalytic arrangements 
of residues obstructed by other parts of the protein structure are possible.
When analyzing predicted protein structures, particular caveats apply.
Template matching is particularly sensitive to conformation.
Poorly predicted side chain rotamers, unrealistic tertiary structure arrangements
and overall protein conformations can severely limit matches **EnzyMM** can find.
Keep in mind that pLDDT does not necessarily reflect on side-chain conformations.
Often structures will be relaxed through energy minimization which may impact results.
