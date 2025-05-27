echo "
============================================================================================      CYPWare 1.0 script and HTML tool     ==========================================================================================================================
The oxd-CYPWare is the modified script of CYPWare 1.0. The oxd-CYPWare allows teh user to perform the MD simulations for the CYP450 oxidized state.

Please cite the use of CYPWare 1.0 and the HTML Tools as follows:
Dixit, V. A.; USN Murty, Bajaj, P.; Blumberger, J.; and Sam P. de Visser. Mechanisms of Electron Transfer Rate Modulations in Cytochrome P450 BM3. J. Phys. Chem. B. 2022 126 (47), 9737-9747. DOI: 10.1021/acs.jpcb.2c03967
https://pubs.acs.org/doi/10.1021/acs.jpcb.2c03967

WELCOME to oxd-CYPWare web interface to perform the MD simulations for the CYP450 oxidized state.  
This script requires the following information for successful calculations (also note the dependencies on other software which are required to be in the USER-PATH are given below).

This user-friendly HTML tool (oxd-CYPWare html) can be used to put all the required information in a text file which CYPWare 1.0 (this script) will read direclty.

If you have these requirements ready, then proceed with the next steps, else perform the docking calculation first and then come back to this oxd-CYPWare web interface.
1) Name the ligand (without file extension) that was docked into CYP450 active site using Autodock Vina.

2) Name the protein (without file extension) that was used for the Autodock Vina docking simulations.

3) Unix path on your server/workstation where you wish the new files should be created.  You should have write-permissions for that folder. It could be your home directory or any of the sub-directories.

4) Name the folders of the oxidized state of the Heme center. These directory will be created by the script.

5) Unix path where all the parameter files for the oxidized state are kept.  These are essential for successfully running the MD simulations.  
	These are kept in a folder called parameters on the GitHub page (in case of query you can post a question to the CYPWare 1.0/oxd-CYPWare development team).

6) Posenumber for the ligand docked pose for which a complex will be created followed by MD simulations.  
	Currently, oxd-CYPWare/CYPWare 1.0 takes Autodock Vina output (pdbqt file) to create the complex, but it can be modified by the user to accept other file formats (which obabel can read).

7) Net molecular charge on the ligand.  This is required to successfully generate ligand parameters with the antechamber program and is required for MD simulations.

The methodology underlying the CYPWare is developed by the PI: Dr. Vaibhav A. Dixit, Asst. Prof., Dept. of Med. Chem., NIPER Guwahati in the Advanced Center of Computer-Aided Drug Design (A-CADD).
The web interface (GUI) is developed in collaboration with C-DAC (Dr. Vinod, Mr. Saurabh, and the team)
CYPWare software, GUI, websites, tools, layouts, and logos are subject to copyright protection and are the exclusive property of the Department of Medicinal Chemistry, NIPER Guwahati. 
The name and logos of the NIPER Guwahati website must not be associated with publicity or business promotion without NIPER G's prior written approval. 

CYPWare 1.0 is available under the creative commons license and is free for academic research groups working in a degree-granting university/institute.  
Any work/report/thesis/research-article/review-article resulting from the use of CYPWare 1.0 should properly cite the software and publication associated with the same.


===========================================================================  Dependencies  ===========================================================================================
oxd-CYPWare makes use of the following opensource tools, thus these need to be installed first from GitHub, Sourceforge, or as a conda package.
Ensure that dependencies are satisfied before running CYPWare 1.0, else the calculation will not complete as expected.

1) Openbabel 3.1 or higher 

	Openbabel is available as a conda package and can be installed with one of the following commands.

	conda install -c conda-forge openbabel
	conda install -c conda-forge/label/cf202003 openbabel

If you don't have conda, then install it from the main website https://www.anaconda.com/
Instructions for conda installation can be found on its website.

2) AmberTools and Amber18 or higher
	
	AmberTools is a freely available software used for the setup and analysis of MD simulations. It is available from the http://ambermd.org/ website.
	It can also be installed as a conda package with any one of the following command.
	
	conda install -c conda-forge ambertools
	conda install -c conda-forge/label/cf202003 ambertools

	Amber18, 20 or the latest 24 version is a widely used MD engine and includes sander, pmemd, and their serial, parallel (MPI), and GPU (cuda) versions.
	It is available at a reasonable price from Prof. David Case's group at UCSF http://ambermd.org/GetAmber.php#amber

	AMBERHOME directory should be in your path for oxd-CYPWare to run correctly

3) AmberMdPrep is a wrapper script by Daniel R. Roe used for equilibrating a biomolecular system for MD simulations

	It can be downloaded from https://github.com/drroe/AmberMdPrep
	Untar the AmberMdPrep folder and make sure that the AmberMdPrep.sh script is in your path
	
	AmberMdPrep depends on the latest GitHub version of cpptraj which is available here https://github.com/Amber-MD/cpptraj
	The GitHub version of cpptraj should be installed in a separate directory and sourced before running AmberMdPrep.  CYPWare 1.0 will source cpptraj automatically if the path is set correctly.

4) Statistical tool st
	A simple statistical tool available on GitHub is used to extract the averages and standard deviations in vertical energy gaps.
	This is available at https://github.com/nferraz/st
	st can be installed using the following commands

	git clone https://github.com/nferraz/st
	cd perl5
	perl Makefile.PL
	make
	make test
	make install

5) Extracting Marcus ET parameters
	This can be done by calling additional scripts on the Linux terminal.
	bash get-Marcus-ET-parm-LB.sh file (for ligand-bound state)
	bash get-Marcus-ET-parm-LF.sh file (for ligand-free state)

For commercial usage of CYPWare 1.0, please contact the PI at vaibhavadixit@gmail.com or vaibhav@niperguwahati.in
============================================================================================================================================================================================================================================= "

source /home/$USER/.bashrc
source $AMBERHOME/amber.sh
module load cuda/10.1
module load amber/24
module load openbabel/3.1.1
cwd=$(pwd)
export $(xargs <$1)
echo Reading $1 file for job parameters
echo "ligand and protein used for MD are $ligand $protein " > $protein-$ligand-MD-logfile.txt
echo "$dirpath used in the calculation" >> $protein-$ligand-MD-logfile.txt
echo "$oxddir used in the calculation" >> $protein-$ligand-MD-logfile.txt

mkdir $dirpath/$oxddir

echo "created directories oxd state " >> $protein-$ligand-MD-logfile.txt
echo "parameters oxd states read-from $oxdparm " >> $protein-$ligand-MD-logfile.txt

echo "Extracting ligand and protein Changing directory to Ferric or the oxidized state i.e. $dirpath/$oxddir"

ligresname=$(awk '/^HETATM/ {print $4}' $ligand.pdb | sort | uniq)
if [ -z "$ligresname" ]; then $ligresname=UNK; else echo $ligresname; fi
ligresname=UNK

egrep -v " ALA | ARG | ASN | ASP | CYS | GLN | GLU | GLY | HID | HIE | HIP | HIS | ILE | LEU | LYS | MET | PHE | PRO | SER | THR | TRP | TYR | VAL | HEM | HM1 | CM1 | FE " $protein.pdb > ligand.pdb

cp $ligand.pdb $dirpath/$oxddir
cp $protein.pdb $dirpath/$oxddir
cp 4ZF6_md_mcpbpy.frcmod $dirpath/$oxddir/
cd $dirpath/$oxddir
echo "parameterizing ligand and generating the protein complex"

obabel $ligand.pdb -O $ligand.pdb -m 

for i in $(ls $ligand*.pdb | grep -v $posenumber ); do echo removing file $i; rm $i ; done

echo "pose number $posenumber from the docking output $ligand.pdb file-was used forMD simulations " >> ../$protein-$ligand-MD-logfile.txt

echo "pose number $posenumber from the docking output $ligand.pdb file-was used forMD simulations " 


obabel $ligand"$posenumber".pdb -O $ligand"$posenumber"-h.pdb -h

grep HETATM $ligand"$posenumber"-h.pdb > $ligand-MD.pdb
echo "Give the expeced net-charge on the molecule as an integer e.g. -1, 0 or 1"
#babel test-$ligand.pdb -O test.mol2 --partialcharge gasteiger
#molchrg=$(egrep -A1000 ATOM test.mol2 | grep -v ROOT | grep LIG | sed -e "1d" | awk '{print $9}' | st | awk '{print $4}' | tail -1 | xargs printf "%1.0f")

$AMBERHOME/bin/antechamber -i $ligand-MD.pdb -fi pdb -o $ligand.mol2 -fo mol2 -c bcc -s 2 -nc $molchrg
$AMBERHOME/bin/parmchk2 -i $ligand.mol2 -f mol2 -o $ligand.frcmod
$AMBERHOME/bin/antechamber -i $ligand.mol2 -fi mol2 -o $ligand.ac -fo ac 
$AMBERHOME/bin/antechamber -i $ligand.ac -fi ac -o $ligand-test.pdb -fo pdb 

echo "Antechamber used to prepare ligand parameters and atomtypes  $ligand " >> ../$protein-$ligand-MD-logfile.txt
echo "Antechamber used to prepare ligand parameters and atomtypes  $ligand "

sed 's/[[:space:]]*[^[:space:]]*$//' $ligand-test.pdb | sed 's/[[:space:]]*[^[:space:]]*$//'  > $ligand-MD1.pdb
#rm test.pdb
echo " "
echo "Combining $protien $ligand into a complex" | tee -a $protein-$ligand-MD-logfile.txt

cat $protein.pdb $ligand-MD1.pdb > $protein-$ligand.pdb
echo "unloading abmer20" | tee -a $protein-$ligand-MD-logfile.txt
module unload new_amber/20
echo "loading amber/24" | tee -a $protein-$ligand-MD-logfile.txt
module load amber
echo "Preparing the complex - Amber usage" | tee -a $protein-$ligand-MD-logfile.txt
$AMBERHOME/bin/pdb4amber -i $protein-$ligand.pdb -o $protein-$ligand-comp.pdb
echo "copying files $ligand $protein into $oxddir" | tee -a $protein-$ligand-MD-logfile.txt
cp -f $ligand.mol2 $dirpath/$oxddir
cp -f $ligand.frcmod $dirpath/$oxddir
cp -f $protein-$ligand-comp.pdb $dirpath/$oxddir

echo "loading amber" | tee -a $protein-$ligand-MD-logfile.txt
module load amber

echo "$protein-$ligand complex formed and saved - MD simulations " | tee -a ../$protein-$ligand-MD-logfile.txt
echo "returning to the partent directory" | tee -a $protein-$ligand-MD-logfile.txt
cd ..
echo "Writing tleap input-file the ferric or oxidized state" | tee -a $protein-$ligand-MD-logfile.txt
echo "
source leaprc.protein.ff14SB
source leaprc.gaff2
addAtomTypes { " > $dirpath/$oxddir/$protein-$ligand-tleap.in
echo '        { "M1"  "Fe" "sp3" }' >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo '        { "Y1"  "S" "sp3" }' >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo '        { "Y2"  "N" "sp3" }' >> /$dirpath/$oxddir/$protein-$ligand-tleap.in
echo '        { "Y3"  "N" "sp3" }' >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo '        { "Y4"  "N" "sp3" }' >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo '        { "Y5"  "N" "sp3" }' >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo '}' >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo "CM1 = loadmol2 $oxdparm/CM1.mol2
HM1 = loadmol2 $oxdparm/HM1.mol2
FE1 = loadmol2 $oxdparm/FE1.mol2

$ligresname = loadmol2 $dirpath/$oxddir/$ligand.mol2
loadamberparams $oxdparm/HEM.frcmod
loadamberparams $dirpath/$oxddir/$ligand.frcmod
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams $oxdparm/4ZF6_mcpbpy.frcmod

source leaprc.water.tip3p
mol = loadpdb $dirpath/$oxddir/$protein-$ligand-comp.pdb
bond mol.425.SG mol.482.FE
bond mol.481.NA mol.482.FE
bond mol.481.NB mol.482.FE
bond mol.481.NC mol.482.FE
bond mol.481.ND mol.482.FE
bond mol.424.C mol.425.N
bond mol.425.C mol.426.N
savepdb mol $dirpath/$oxddir/$protein-$ligand-dry.pdb
saveamberparm mol $dirpath/$oxddir/$protein-$ligand-dry.prmtop $dirpath/$oxddir/$protein-$ligand-dry.inpcrd
charge mol
solvatebox mol TIP3PBOX 10.0
addions mol K+ 0
addions mol Cl- 0
charge mol
savepdb mol $dirpath/$oxddir/$protein-$ligand-solv.pdb
saveamberparm mol $dirpath/$oxddir/$protein-$ligand-solv.prmtop $dirpath/$oxddir/$protein-$ligand-solv.inpcrd
quit" >> $dirpath/$oxddir/$protein-$ligand-tleap.in
echo "tleap input file has been created for the ferric or oxidized state. Please check if the system charged is 0, if not make appropriate changes in the type and number of ions to be added in the addions file and rerun this main script"




	echo "Writing MD input - a 20 ns MD production run"
	echo "Production Stage - Explicit Solvent 20 ns
&cntrl
  ntt=3,           ! Temperature scaling (=3, Langevin dynamics)
  gamma_ln=2.0,    ! Collision frequency of the Langevin dynamics in ps-1
  ntc=2,           ! SHAKE constraints (=2, hydrogen bond lengths constrained)
  ntf=2,           ! Force evaluation (=2, hydrogen bond interactions omitted)
  ntb=1,           ! Boundaries (=1, constant volume)
  cut=10.0,        ! Cutoff
  dt=0.002,        ! The time step in picoseconds
  nstlim=10000000, ! Number of MD steps to be performed
  ig=-1,           ! Random seed (=-1, get a number from current date and time)
  ntwr=5000,      ! Restart file written every ntwr steps
  ntwx=5000,      ! Trajectory file written every ntwx steps
  ntpr=5000,      ! The mdout and mdinfo files written every ntpr steps
  ioutfm=1,        ! Trajectory file format (=1, Binary NetCDF)
  iwrap=1,         ! Translate water molecules into the original simulation box
  igb=0,           ! GB model (=0, explicit solvent)
  irest=1,         ! Flag to restart the simulation
  ntx=5,           ! Initial condition (=5, coord. and veloc. read from the inpcrd file)
/
" > $dirpath/$oxddir/prod20ns.in

#Ferrousdir=$(echo Ferrous-$protein-$ligand | sed 's/_out//g')

	echo "
parm $protein-$ligand-solv.prmtop
trajin $protein-$ligand-solv-prod.nc
autoimage
center
reference final.1.ncrst

nativecontacts name $protein-$ligand-LIG-FE :482&!@H= :483&!@H= byresidue out $protein-$ligand-LIG-FEcontacts.dat mindist maxdist distance 10.0 first map mapout $protein-$ligand-LIG-FEresmap.gnu contactpdb $protein-$ligand-LIG-FEcontactmap.pdb series seriesout $protein-$ligand-LIG-FEnativecontacts.dat writecontacts $protein-$ligand-LIG-FEwritecontacts.dat 
nativecontacts name $protein-$ligand-LIG-HEM :481&!@H= :483&!@H= byresidue out $protein-$ligand-LIG-HEMcontacts.dat mindist maxdist distance 10.0 first map mapout $protein-$ligand-LIG-HEMresmap.gnu contactpdb $protein-$ligand-LIG-HEMcontactmap.pdb series seriesout $protein-$ligand-LIG-HEMnativecontacts.dat writecontacts $protein-$ligand-LIG-HEMwritecontacts.dat " > $dirpath/$oxddir/contacts.in

echo "
parm $protein-$ligand-solv.prmtop
trajin $protein-$ligand-solv-prod-rp.nc
autoimage
center
reference final.1.ncrst

nativecontacts name $protein-$ligand-LIG-FE :482&!@H= :483&!@H= byresidue out $protein-$ligand-LIG-FEcontacts.dat mindist maxdist distance 10.0 first map mapout $protein-$ligand-LIG-FEresmap.gnu contactpdb $protein-$ligand-LIG-FEcontactmap.pdb series seriesout $protein-$ligand-LIG-FEnativecontacts.dat writecontacts $protein-$ligand-LIG-FEwritecontacts.dat 
nativecontacts name $protein-$ligand-LIG-HEM :481&!@H= :483&!@H= byresidue out $protein-$ligand-LIG-HEMcontacts.dat mindist maxdist distance 10.0 first map mapout $protein-$ligand-LIG-HEMresmap.gnu contactpdb $protein-$ligand-LIG-HEMcontactmap.pdb series seriesout $protein-$ligand-LIG-HEMnativecontacts.dat writecontacts $protein-$ligand-LIG-HEMwritecontacts.dat " > $dirpath/$oxddir/contacts-rp.in

	echo "
parm $protein-$ligand-solv.prmtop
trajin $protein-$ligand-solv-prod-rp1.nc
autoimage
center
reference final.1.ncrst

nativecontacts name $protein-$ligand-LIG-FE :482&!@H= :483&!@H= byresidue out $protein-$ligand-LIG-FEcontacts.dat mindist maxdist distance 10.0 first map mapout $protein-$ligand-LIG-FEresmap.gnu contactpdb $protein-$ligand-LIG-FEcontactmap.pdb series seriesout $protein-$ligand-LIG-FEnativecontacts.dat writecontacts $protein-$ligand-LIG-FEwritecontacts.dat 
nativecontacts name $protein-$ligand-LIG-HEM :481&!@H= :483&!@H= byresidue out $protein-$ligand-LIG-HEMcontacts.dat mindist maxdist distance 10.0 first map mapout $protein-$ligand-LIG-HEMresmap.gnu contactpdb $protein-$ligand-LIG-HEMcontactmap.pdb series seriesout $protein-$ligand-LIG-HEMnativecontacts.dat writecontacts $protein-$ligand-LIG-HEMwritecontacts.dat " > $dirpath/$oxddir/contacts-rp1.in


leapinput=$protein-$ligand-tleap.in
leapoutput=$protein-$ligand-tleap.out
#echo "give base file names to be used for MD prmtop, inpcrd, out, rst7, mdinfo and trajectory (nc) files"
basefile=$protein-$ligand-solv
echo $leapinput $leapoutput $basefile | tee -a $protein-$ligand-MD-logfile.txt
echo "creating MD prmtop inpcrd files for the Ferric or the oxidized state in $dirpath/$oxddir " | tee -a $protein-$ligand-MD-logfile.txt
echo "creating gpu-job-head" | tee -a $protein-$ligand-MD-logfile.txt
echo "#!/bin/bash
#SBATCH -A meity
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --nodelist=node1
#SBATCH --time=95:50:20
#SBATCH --job-name=$ligand-$protein
#SBATCH --error=job.%J.err_node_40
#SBATCH --output=job.%J.out_node_40
#SBATCH --partition=GPU_NODES     ##partition name
#SBATCH --gres=gpu:1
myname=$USER

source /home/software/amber24/amber.sh
module load openmpi/4.1.4
module load amber

" > gpu-job-head.txt


echo "creating cpu-job-head " | tee -a $protein-$ligand-MD-logfile.txt
echo "#!/bin/bash
#SBATCH -A meity     ##account name
#SBATCH --nodes=1  ## number of nodes
#SBATCH -p CPU_NODES   #partition name
#SBATCH --ntasks-per-node=8  ## number os  processor core(s) per node
#SBATCH --time=995:50:20      ##time optional
#SBATCH --job-name=$ligand-$protein   ##job name
#SBATCH --error=job.%J.err    ## any error during job submission or execution
#SBATCH --output=job.%J.out   ##any output after job execution
myname=$USER


" > cpu-job-head.txt
echo "creating tleap input and sbatch - $oxddir state" | tee -a $protein-$ligand-MD-logfile.txt
cat cpu-job-head.txt > $dirpath/$oxddir/$protein-$ligand-tleap.batch
echo "module load amber" >> $dirpath/$oxddir/$protein-$ligand-tleap.batch
echo "source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/$protein-$ligand-tleap.batch
echo "$AMBERHOME/bin/tleap -s -f $dirpath/$oxddir/$leapinput > $dirpath/$oxddir/$leapoutput " >> $dirpath/$oxddir/$protein-$ligand-tleap.batch


oxdtleapjobid=$(sbatch $dirpath/$oxddir/$protein-$ligand-tleap.batch | awk '{print $4}')
echo "sleeping - 1m"
sleep 1m
oxdtleapjobstatus=$(sacct | grep $oxdtleapjobid | awk '{print $6}' | head -1)



echo "creating input and submitting equil jobs" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt

while [[ $redtleapjobstatus == RUNNING ]];
do
	sleep 1m
done
echo "finished equiliration of $oxddir structure" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt

source /home/vaibhav/cpptraj-master/cpptraj.sh
echo "changing directory to $oxddir" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cd $dirpath/$oxddir
cp ../gpu-job-head.txt .
echo "creating MDequil sbatch files and submitting the job $oxddir state" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cat gpu-job-head.txt > $dirpath/$oxddir/$protein-$ligand-MDequil.batch
echo "module load amber" >> $dirpath/$oxddir/$protein-$ligand-MDequil.batch
echo "source $AMBERHOME/amber.sh
export CUDA_VISIBLE_DEVICES=0" >> $dirpath/$oxddir/$protein-$ligand-tleap.batch
echo "source /home/vaibhav/cpptraj-master/cpptraj.sh" >> $dirpath/$oxddir/$protein-$ligand-MDequil.batch
echo "/home/anila/anila/md/AmberMdPrep.sh -O -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/$basefile.inpcrd --temp 300 --ares HM1 --ares FE1 --ares CM1 --ares LIG " >> $dirpath/$oxddir/$protein-$ligand-MDequil.batch
oxdequiljobid=$(sbatch $dirpath/$oxddir/$protein-$ligand-MDequil.batch | awk '{print $4}')
oxdequiljobstatus=$(sacct | grep $oxdequiljobid | awk '{print $6}' | head -1)


echo "sleeping till the equiljobs are finish"

while [[ $(sacct | grep $oxdequiljobid | awk '{print $6}' | head -1 ) == RUNNING ]];
do
	sleep 1m
done
echo "finished equiliration of $oxddir structure" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt

echo "changing directory to $oxddir and submitting 20ns MD-jobs - triplicate (prod, prod-rp and prod-rp1) after equilibriation is success" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cd $dirpath/$oxddir
cat gpu-job-head.txt > $dirpath/$oxddir/oxdMD-rp-rp1.batch
echo "module load amber" >> $dirpath/$oxddir/oxdMD-rp-rp1.batch
echo "source $AMBERHOME/amber.sh
export CUDA_VISIBLE_DEVICES=0" >> $dirpath/$oxddir/oxdMD-rp-rp1.batch
echo "pmemd.cuda -O -i $dirpath/$oxddir/prod20ns.in -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/final.1.ncrst -ref $dirpath/$oxddir/final.1.ncrst -r $dirpath/$oxddir/$basefile-prod.rst7 -inf $dirpath/$oxddir/$basefile-prod.mdinfo -o $dirpath/$oxddir/$basefile-prod.mdout -x $dirpath/$oxddir/$basefile-prod.nc 
pmemd.cuda -O -i $dirpath/$oxddir/prod20ns.in -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/final.1.ncrst -ref $dirpath/$oxddir/final.1.ncrst -r $dirpath/$oxddir/$basefile-prod-rp.rst7 -inf $dirpath/$oxddir/$basefile-prod-rp.mdinfo -o $dirpath/$oxddir/$basefile-prod-rp.mdout -x $dirpath/$oxddir/$basefile-prod-rp.nc 
pmemd.cuda -O -i $dirpath/$oxddir/prod20ns.in -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/final.1.ncrst -ref $dirpath/$oxddir/final.1.ncrst -r $dirpath/$oxddir/$basefile-prod-rp1.rst7 -inf $dirpath/$oxddir/$basefile-prod-rp1.mdinfo -o $dirpath/$oxddir/$basefile-prod-rp1.mdout -x $dirpath/$oxddir/$basefile-prod-rp1.nc 
cpptraj.cuda -i $dirpath/$oxddir/contacts.in
cpptraj.cuda -i $dirpath/$oxddir/contacts-rp.in
cpptraj.cuda -i $dirpath/$oxddir/contacts-rp1.in 

" >> $dirpath/$oxddir/oxdMD-rp-rp1.batch 
oxdMDjobid=$(sbatch --dependency=afterok:$oxdequiljobid $dirpath/$oxddir/oxdMD-rp-rp1.batch | awk '{print $4}')
oxdMDjobstatus=$(sacct | egrep $oxdMDjobid | awk '{print $6}'  | head -1)




while [[ $(sacct | grep $oxdMDjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 5m
done



echo "changing directory to $oxddir and submitting cpptraj, MD trajectories (prod, prod-rp and prod-rp1) after production runs is success" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cd $dirpath/$oxddir


oxdspMDjobstatus=$(sacct | egrep $oxdspMDjobid | awk '{print $6}'  | head -1 )



oxdspMDprotjobstatus=$(sacct | egrep $oxdspMDprotjobid | awk '{print $6}'  | head -1 )







