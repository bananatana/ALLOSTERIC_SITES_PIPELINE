# GPCR ALLOSTERIC SITE SEARCH PIPELINE #

A framework that uses a combination of different scripts to analyze allosteric connectivity within GPCR proteins and and identifies possible allosteric sites. 
The protocol is divided into several stages:

**1.** System preparation for MD simulations  
**2.** MD simulations  
**3.** Processing trajectories for analysis  
**4.** Calculation of configurational entropy and mutual information  
**5.** Analysis of mutual information  
**6.** Identification of allosteric sites.  


## 1. System preparation for MD simulations ##

The first thing that is needed is a high-quality (and complete) structure of the receptor we want to map for allosteric sites. 
I used refined structure ([Alphafold](https://alphafold.ebi.ac.uk/)) with [GPCRdb](https://gpcrdb.org/). Be sure to check if the structure makes sense (hyper structuring may occur when using Alphafold).
[Modeller](https://salilab.org/modeller/) can also be used for smaller reconstructions. If you need to add cysteine ​​bridges while using Modeller, the script for modeling looks like this:

~~~
from modeller import *
from modeller.automodel import *    # Load the AutoModel class
log.verbose()
env = Environ()
# directories for input atom files
env.io.atom_files_directory = ['.','.']
class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self.residue_range('150:A', '168:A'),      #missing loop
                         self.residue_range('218:A', '223:A'))      #missing loop
# Redefine the special_patches routine to include the additional disulfides
# (this routine is empty by default):
class MyModel(AutoModel):
    def special_patches(self, aln):
        # A disulfide between residues 163 and 69 in chain A:
        self.patch(residue_type='DISU', residues=(self.residues['163:A'],
                                                  self.residues['69:A']))
        # A disulfide between residues 151 164 in chain A:
        self.patch(residue_type='DISU', residues=(self.residues['151:A'],
                                                  self.residues['164:A']))
a = MyModel(env, alnfile = './alignment.ali',
            knowns = 'Structure', sequence = 'Stucture_fill')
a.starting_model= 1
a.ending_model  = 5
a.make()
~~~

If you want to map the active conformation of GPCR, make sure that the initial model should contain G protein (at least alpha subunit). From experience, the active conformation is not stable during MD simulation without G protein.


