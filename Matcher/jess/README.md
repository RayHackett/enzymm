## Jess: A 3D template searching program

Jess is Copyright (c) of Jonathan Barker, 2002.  
Now maintained and modified by Ioannis Riziotis (riziotis@ebi.ac.uk).  
For details on the algorithm, please read the [paper](10.1093/bioinformatics/btg226).  

Last updated 14/12/2022

### Installation

Compile with something like:

`cd src`  
`gcc -c *.c`  
`gcc -o jess *.o -lm `  
`sudo mv jess /usr/local/bin`  

### Usage

`jess [template-list] [target-list] [rmsd] [distance] [max-dynamic-distance] [conservation-cutoff] [flags]`

* `template-list`: a list of filenames of TESS templates
* `target-list`: is a list of filenames of PDB files to search
* `rmsd`: the RMSD cutoff at which results are reported
* `distance`: the global distance cutoff used to guide the search
* `max-dynamic-distance`: the maximum allowed template/query atom distance 
			  after adding the global distance cutoff and the 
			  individual atom distance cutoff defined in the
			  temperature field of the ATOM record in the template
			  file.
* `conservation-cutoff`: a cutoff defined in the b-factor field of the query structure.
                         Atoms whose value in the b-factor field is lower than this cutoff
                         will be ignored.

As a rough estimate the distance cutoff should be 1.5-4 times
the RMSD cutoff. But if you make it very large execution 
time will suffer. The smaller the better!

[flags] : optional flags as a string with no spaces:  
* `f` : see PDB filenames in progress on stderr  
* `n` : do not transform coordinates of hit into	the template coordinate frame  
* `i` : include matches composed of residues belonging to
      multiple chains (if template is single-chain), or
	  matches with residues from a single chain
	  (if template has residues from multiple chains)  
* `q` : write filename of query instead of PDB ID from HEADER  
* `e` : parse atoms from all models separated by ENDMDL (use with
	  care). By default, Jess will only parse the first model

Example:

`cd examples`  
`../jess templates testfiles 2 3 3 0 > output`  

The output file is a flat file containing PDB fragments 
preceded by a single record of the form

`REMARK pdb-code rmsd template-file [some debugging info]`

where pdb-code is the PDB code of the hit file, rmsd is 
the rmsd after optimal superposition with the template and
template-file is the file containing the template which 
was hit.

The debugging info currently contains Det=number and
log(E)~number. If Det is not 1.0 then the superposition
is not valid (tell me about it please!). log(E) is a
preliminary statistical measure which should be used
with caution. Anything less than -4 should be a very
good hit and -3 is OK. (E is the expected number of hits
at random).

The PDB fragment which follows the remark is transformed 
into the coordinate frame of the template by optimal
superposition (if [n] flag is not there).

Each hit is followed by ENDMDL and a blank line.

### Filtering the output

Please note that in some cases, Jess performs multiple 
optimal aligments at a specific atom set in a given 
template-target pair. If you want to keep only the best hit,
pipe the output to the 'filter_jessout.py' script.

