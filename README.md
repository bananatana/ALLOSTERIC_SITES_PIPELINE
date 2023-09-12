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

### 1.1. Assembling system for MD simulations ###
When you have all the necessary (and complete) molecules for MD simulation (GPCR, G protein, orthosteric ligand, allosteric ligand, ions, crystal waters, cholesterol molecules), first you need to align everything using the [PPM web server](https://opm.phar.umich.edu/ppm_server) so that the placement in the membrane is correct. 
After that you can continue with the [PYMEMDYN](https://github.com/GPCR-ModSim/pymemdyn) protocol.

## 1. MD simulations ##
After implementing the [PYMEMDYN](https://github.com/GPCR-ModSim/pymemdyn) protocol, it is necessary to perform  at least 1μs long MD simulation using [GROMACS](https://www.gromacs.org/). The input file looks like this:

~~~
; title =  Production run for PARENT
;include             = -I./posre_hoh.itp
;define              =  -DPOSRES ; Restraint heavy atoms of the protein (all except H's)
; Run Parameters
integrator          =  md        ; leap-frog integrator
dt                  =  0.002     ; in ps !  = 2 fs
nsteps              =  100000000 ; total 200 ns
; Bond Parameters
constraints         =  all-bonds ; constraint all bonds using LINCS ; SETTLE for water
; Output Parameters
nstxout             =  50000     ; save coordinates every 100 ps
nstvout             =  25000     ; save velocities every 50 ps
nstenergy           =  25000     ; save energies every 50 ps
nstxtcout           =  2500      ; save xtc trajectory every 100 ps
; Neighbor Searching Parameters
nstlist             =  5         ; update the neighbor list every 10 fs
                                 ; This works with twin-cutoff (if rlist < rcoulomb)
ns_type             =  grid      ; It makes a grid in the box for neighbor list searching
rlist               =  1.2       ; = rcoulomb with PME
rcoulomb            =  1.2       ; real-space cutoff
rvdw                =  1.2       ; short range vdw cuoff
; Electrostatic Parameters
coulombtype         =  PME       ; Particle mesh Ewald for LR interactions
fourierspacing      =  0.15      ; grid dimensions for FFT
;;ewald_geometry      =  3dc       ; only for slab geometries
pme_order           =  4         ; interpolation order for FFT
ewald_rtol          =  1e-5      ; relative accuracy of direct/reciprocal space
optimize_fft        =  yes
; Temperature Coupling
Tcoupl              =  nose-hoover           ; thermostat
tau_t               =  0.5   0.5  0.5        ; time constaint in ps
tc-grps             =  wation protlig membr  ; coupling groups
ref_t               =  310   310   310       ; reference temperature
; Pressure Coupling
Pcoupl              =  Parrinello-Rahman     ; barostat
tau_p               =  2.0
compressibility     =  4.5e-5       4.5e-5
ref_p               =  1.0          1.0
pcoupltype          =  semiisotropic
; Generation of Velocities
;gen_vel             =  yes
;gen_temp            =  310
;gen_seed            =  -1
~~~

After modifying prod.mdp run do grompp:
~~~
gmx grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr --maxwarn 2
~~~

Take into account that due to the size of the system you probably won't be able to run the entire length of the trajectory at once; it is optimal to do it five times for 200 ns.
Now you can run your simulation. Here I provide example of my running script from Berzelius:
~~~
#!/bin/bash -l
#SBATCH --time=3-00:00:00
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --gpus 1
#SBATCH --job-name=allosteric
export GMX_BIN=/proj/compbiochem/users/bin/gromacs-2022.4/bin/gmx
export RUNDIR=$TMPDIR/$SLURM_JOB_ID
echo "Host: $(hostname)"
echo "Tmpdir: $RUNDIR"
echo "Jobdir: $SLURM_SUBMIT_DIR"
# copy files to scratch dir
rsync -ah $SLURM_SUBMIT_DIR/ $RUNDIR/ --exclude="slurm*" --exclude="*.sh"
cd $RUNDIR
${GMX_BIN} mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c confout.gro -g production.log -x traj_prod.xtc -cpo mdrun.cpt
rsync -ah --update  $RUNDIR/ $SLURM_SUBMIT_DIR/
## cleanup
rm -rf $RUNDIR
echo "Done"
~~~

As first increment of simulation finishes, you extend your simulation with desired number of steps:
~~~
gmx convert-tpr -s topol_prod.tpr -extend 200000 -o prod_new.tpr #extend for 200ns
~~~
And then run from the checkpoint file with the same running script:
~~~
gmx mdrun -s prod_new.tpr -o traj.trr -e ener.edr -c confout_new.gro -g production.log -x traj_prod.xtc -cpo mdrun_new.cpt -cpi mdrun.cpt
~~~

After 1μs trajectory is ready for quality control and analysis. 






