# carrier-mobility
Code to estimate carrier mobility in organic semiconductor dimers using gaussian09 output files by Marcus-Hush and Marcus-levitch-Jortner models.

In order to run this code you need to download and compile other code also from github to generate transfer integrals and follow a naming and location pattern for the input files. Once you have compiled the other software just paste its executable to this root.

Software for transfer integral calculations link(CATNIP): https://github.com/JoshuaSBrown/QC_Tools

**Citing**
If you use this code to get values for a paper it would be nice if you cited the original paper https://doi.org/10.1016/j.cplett.2020.138226 related to this code.

**Running the code**
Running the code needs some Gaussian09 .log files and these files needs to addressed in the code.
As a simple example to illustrate the code functionality we use a ethylene dimmer electron mobility calculation, there you have two input folders. One for the CATNIP software inputs (for more information about them go to software main page) the other folder named "cm_inputs" is the one with .log files needed here.

**Preparing input files**
It is provided Gaussian .com input files used to generate the logs used by the code.

**Neutral state**
First we need to optimize monomers geometries for neutral and excited states, use "freq" keyword to check convergence.

"ethylene-neutral.com"

####################################
%NProcShared=4
%Mem=2GB
%chk=ethylene-neutral

# opt B3LYP/6-31G* freq=noraman int=ultrafine

ethylene neutral

0 1
 C	 0.672749	 0.000000	0.000000
 C	-0.672749	 0.000000	0.000000
 H	 1.242623	 0.934806	0.000000
 H	 1.242623	-0.934806	0.000000
 H	-1.242623	 0.934806	0.000000
 H	-1.242623	-0.934806	0.000000
 ####################################

 **Excited state**

Gaussian input for excited state, anion or cation depending on your needs, in this example we have an anion state (electron mobility). Usually is better to optimize the neutral form and use its geometry as input to optimize the excited state, as you can notice, this input has the oldchk being the .chk output from neutral optimization:

"ethylene-anion.com"

####################################
%NProcShared=4
%Mem=2GB
%oldchk=ethylene-neutral
%chk=ethylene-anion

# opt B3LYP/6-31G* freq=noraman int=ultrafine

ethylene anion

-1 2
C	 0.672749	 0.000000	0.000000
C	-0.672749	 0.000000	0.000000
H	 1.242623	 0.934806	0.000000
H	 1.242623	-0.934806	0.000000
H	-1.242623	 0.934806	0.000000
H	-1.242623	-0.934806	0.000000
####################################

**Displacement vector files**
With geometry optimizations converged we need to calculate the displacement vector in both directions, from neutral to anion and anion to neutral:
The input from neutral to anion is the ethylene-neutral-displacement.com:

####################################
%NProcShared=4
%Mem=2GB
%oldchk=ethylene-neutral.chk
%chk=neutral
# B3LYP/6-31+G(d,p) Freq=(SaveNM,noraman) Geom=Checkpoint guess=read

neutral form

0 1


--Link1--
%oldchk=ethylene-anion.chk
%chk=anion
# B3LYP/6-31+G(d,p) Freq=(SaveNM,noraman) Geom=Checkpoint guess=read

anion

-1 2

--Link1--
%oldchk=neutral
%chk=fc
# B3LYP/6-31+G(d,p) Freq=(FC,ReadFCHT) Geom=Check

photoionization

0,1

PRTMAT=12 MAXC1=1000 MAXBANDS=1 SPECHWHM=0.5

anion.chk
####################################

and from anion to neutral is ethylene-anion-displacement.com:
####################################
%NProcShared=4
%Mem=2GB
%oldchk=ethylene-anion.chk
%chk=anion
# B3LYP/6-31+G(d,p) Freq=SaveNM Geom=Check

anion

-1 2


--Link1--
%oldchk=ethylene-neutral.chk
%chk=neutral
# B3LYP/6-31+G(d,p) Freq=SaveNM Geom=Check

anion

0 1

--Link1--
%oldchk=anion
%chk=fc
# B3LYP/6-31+G(d,p) Freq=(FC,ReadFCHT) Geom=Check

photoionization

-1,2

PRTMAT=12 MAXC1=1000 MAXBANDS=1 SPECHWHM=0.5

neutral.chk
####################################

If you pay attention, all displacement inputs use as starting point checkpoint files (.chk) from the optimizations, so you need to keep them.

**Geometry file**
As you have all the logs make a .xyz file with dimmer geometry using the neutral and excited geometries as the two monomers.
####################################
12
ethylene 0
C	 0.672749	 0.000000	0.000000
C	-0.672749	 0.000000	0.000000
H	 1.242623	 0.934806	0.000000
H	 1.242623	-0.934806	0.000000
H	-1.242623	 0.934806	0.000000
H	-1.242623	-0.934806	0.000000
C	 0.672749	 0.000000	5.000000
C	-0.672749	 0.000000	5.000000
H	 1.242623	 0.934806	5.000000
H	 1.242623	-0.934806	5.000000
H	-1.242623	 0.934806	5.000000
H	-1.242623	-0.934806	5.000000
####################################

**Running the code**
In the main code (carrier_mobility.py) you have to write file paths in input files section of the code.
**Files path**
### carrier mobility input files ###
xyzfile =  'cm_inputs/example/ethylene-0-dimer.xyz'              ###  Path of geometry file .xyz
log     =  'cm_inputs/example/ethylene-neutral-displacement.log' ###  Path of neutral to excited displacement Gaussian .log  file
log1    =  'cm_inputs/example/ethylene-anion-displacement.log'   ###  Path of excited to neutral displacement Gaussian .log  file
### CATNIP input files ### For more details check CATNIP page
pun_file   = 'catnip_inputs/example/ethylene-0-dimer.pun'        ###  Path of Dimer pun file for CATNIP
pun_file_1 = 'catnip_inputs/example/ethylene-0-m1.pun'           ###  Path of 1st monomer pun file for CATNIP
pun_file_2 = 'catnip_inputs/example/ethylene-0-m2.pun'           ###  Path of 2nd monomer pun file for CATNIP

**Molecular orbitals**
Choose molecular orbital types for CATNIP transfer integral calculations (ex: LUMO, HOMO...).

orb_ty_1 = 'LUMO'  ###  1st orbital chosen for coupling
orb_ty_2 = 'LUMO'  ###  2nd orbital chosen for coupling
