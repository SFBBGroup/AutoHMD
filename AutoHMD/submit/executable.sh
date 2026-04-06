#!/bin/bash

# Check if the environment modules system is available and load it
if [ -r /usr/share/modules/init/bash ]; then
        source /usr/share/modules/init/bash
fi

# Load the AMBER molecular dynamics software package (version A24)
module load amber/A24

# Define input files for the simulation
prmtop="HMR_solvated_system.prmtop"  # Topology file (contains system structure)
rst="solvated_system.rst7"           # Initial restart/coordinate file

# Create the main replicate directory
mkdir -p rep
folder="$(pwd)/rep"
workingdir="${folder}"

# Navigate to the working directory
cd "${folder}"

# Create subdirectories for organizing simulation outputs
mkdir -p out restart trajectory analysis

# Define paths for output organization
outdir="${workingdir}/out"        # Standard output files (.out)
rstdir="${workingdir}/restart"    # Restart files (.rst)
ncdir="${workingdir}/trajectory"  # Trajectory files (.nc)
analdir="${workingdir}/analysis"   # Analysis files

# Return to the original directory
cd -

echo 'Creating input files for MD simulations...'

# Create 1_water_minimization.in
cat > 1_water_minimization.in << 'EOF'
&cntrl
 imin = 1,
 maxcyc = 10000,
 ncyc = 1000,
 ntr=1, restraintmask='!:WAT',
 restraint_wt=1.0,
&end
EOF

# Create 1_part2_system_minimization.in
cat > 1_part2_system_minimization.in << 'EOF'
&cntrl
 imin = 1,
 maxcyc = 20000,
 ntpr = 100,
 ntr = 0,
 ncyc = 10000,
/
EOF

# Create 2_heat.in
cat > 2_heat.in << 'EOF'
Heating at V constant
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=1,
 ntxo=2,
 irest=0,
 tempi=0.0,
 temp0=310.0,
 dt=0.001,
 nstlim=500000,

 ntb=1,
 ntp=0,

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,

 nmropt=1,
 ntr=1,restraintmask='!:WAT&!@H=',
 restraint_wt=100.0,
/
 &wt type='REST', istep1=0, istep2=0, value1=1.0, value2=1.0 /
 &wt type='TEMP0', istep1=0, istep2=500000, value1=0.0, value2=310.0 /
 &wt type='END',
/
EOF

# Create 3_Water_density.in
cat > 3_Water_density.in << 'EOF'
Water density equilibration - NPT
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=5,
 ntxo=2,
 irest=1,
 tempi=310.0,
 temp0=310.0,
 dt=0.001,
 nstlim=1000000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0
 taup=0.5

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,

 ntr=1,
 restraintmask='!:WAT&!@H=',
 restraint_wt=100.0,
/
EOF

# Create 4_Water_densitY_soft.in
cat > 4_Water_densitY_soft.in << 'EOF'
Heating at V constant
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=5,
 ntxo=2,
 irest=1,
 tempi=310.0,
 temp0=310.0,
 dt=0.001,
 nstlim=1000000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0
 taup=0.5

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,

 ntr=1,
 restraintmask='!:WAT&!@H=',
 restraint_wt=10.0,
/
EOF

# Create 5_2nd_Minimization.in
cat > 5_2nd_Minimization.in << 'EOF'
&cntrl
 imin = 1,
 maxcyc = 10000,
 ncyc = 0,
 ntr=1,
 restraintmask='@N,CA,C,O',
 restraint_wt=10.0,
&end
EOF

# Create 6_First_equil_10kcalmolA2.in
cat > 6_First_equil_10kcalmolA2.in << 'EOF'
First Equilibration - NPT
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=1,
 ntxo=2,
 irest=0,
 tempi=310.0,
 temp0=310.0,
 dt=0.001,
 nstlim=1000000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0
 taup=0.5

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,

 ntr=1,
 restraintmask='@N,CA,C,O',
 restraint_wt=10.0,
/
EOF

# Create 7_Second_equilibration_1kcalmolA2.in
cat > 7_Second_equilibration_1kcalmolA2.in << 'EOF'
Second Equilibration - NPT
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=5,
 ntxo=2,
 irest=1,
 tempi=310.0,
 temp0=310.0,
 dt=0.001,
 nstlim=1000000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0
 taup=0.5

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,

 ntr=1,
 restraintmask='@N,CA,C,O',
 restraint_wt=1.0,
/
EOF

# Create 8_Third_equilibration_01kcalmolA2.in
cat > 8_Third_equilibration_01kcalmolA2.in << 'EOF'
Third Equilibration - NPT
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=5,
 ntxo=2,
 irest=1,
 tempi=310.0,
 temp0=310.0,
 dt=0.001,
 nstlim=1000000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0
 taup=0.5

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,

 ntr=1,
 restraintmask='@N,CA,C,O',
 restraint_wt=0.1,
/
EOF

# Create 9_Unrestr_equilibration.in
cat > 9_Unrestr_equilibration.in << 'EOF'
Unrestrained Equilibration - NPT
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001

 imin=0,
 ig=-1,
 ntx=5,
 ntxo=2,
 irest=1,
 tempi=310.0,
 temp0=310.0,
 dt=0.001,
 nstlim=1000000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0
 taup=0.5

 ntpr=500,
 ntwx=500,
 ntwr=2500,
 nscm=1000,
/
EOF

# Create mdHMR310.in
cat > mdHMR310.in << 'EOF'
Production - NPT - 5 ns 4fs
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001,

 imin=0,
 ig=-1,
 ntx=1,
 ntxo=2,
 irest=0,
 tempi=310.0,
 temp0=310.0,
 dt=0.004,
 nstlim=7500000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0,
 taup=0.5,

 ntpr=3125,
 ntwx=3125,
 ntwr=5000,
 nscm=3125,
/
EOF

# Create mdHMR330.in
cat > mdHMR330.in << 'EOF'
Production - NPT - 5 ns 4fs
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001,

 imin=0,
 ig=-1,
 ntx=1,
 ntxo=2,
 irest=0,
 tempi=330.0,
 temp0=330.0,
 dt=0.004,
 nstlim=3125000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0,
 taup=0.5,

 ntpr=3125,
 ntwx=3125,
 ntwr=5000,
 nscm=3125,
/
EOF

# Create mdHMR360.in
cat > mdHMR360.in << 'EOF'
Production - NPT - 5 ns 4fs
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001,

 imin=0,
 ig=-1,
 ntx=1,
 ntxo=2,
 irest=0,
 tempi=360.0,
 temp0=360.0,
 dt=0.004,
 nstlim=3125000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0,
 taup=0.5,

 ntpr=3125,
 ntwx=3125,
 ntwr=5000,
 nscm=3125,
/
EOF

# Create mdHMR390.in
cat > mdHMR390.in << 'EOF'
Production - NPT - 5 ns 4fs
 &cntrl
 cut=8.0,
 ntt=3,
 gamma_ln=1.0,
 ntc=2,
 ntf=2,
 tol=0.00001,

 imin=0,
 ig=-1,
 ntx=1,
 ntxo=2,
 irest=0,
 tempi=390.0,
 temp0=390.0,
 dt=0.004,
 nstlim=3750000,

 ntb=2,
 ntp=1,
 barostat=1,
 pres0=1.0,
 taup=0.5,

 ntpr=3125,
 ntwx=3125,
 ntwr=5000,
 nscm=3125,
/
EOF

echo 'All input files created successfully!'
echo 'Starting MD simulations...'


# Step 1: Water minimization
# Minimizes only water molecules while keeping solute fixed
pmemd.cuda -O -i 1_water_minimization.in              -p $prmtop -o ${outdir}/water_min_1.out  -c $rst             -r ${rstdir}/water_min_1.rst                      -ref $rst

# Step 2: Full system minimization
# Minimizes the entire system (solute + solvent)
pmemd.cuda -O -i 1_part2_system_minimization.in                               -p $prmtop -o ${outdir}/min.out          -c ${rstdir}/water_min_1.rst  -r ${rstdir}/min.rst                    -ref ${rstdir}/water_min_1.rst

# Step 3: Heating phase
# Gradually raises system temperature to target value
pmemd.cuda -O -i 2_heat.in                            -p $prmtop -o ${outdir}/heat_2.out       -c ${rstdir}/min.rst          -r ${rstdir}/heat_2.rst                           -ref ${rstdir}/min.rst

# Step 4: Water density equilibration
# Adjusts water density at target temperature
pmemd.cuda -O -i 3_Water_density.in                   -p $prmtop -o ${outdir}/water_dens_3.out -c ${rstdir}/heat_2.rst       -r ${rstdir}/water_dens_3.rst -x ${ncdir}/water_dens_3.nc  -ref ${rstdir}/heat_2.rst

# Step 5: Soft water density equilibration
# Continues density equilibration with softer restraints
pmemd.cuda -O -i 4_Water_densitY_soft.in              -p $prmtop -o ${outdir}/water_soft_4.out -c ${rstdir}/water_dens_3.rst -r ${rstdir}/water_soft_4.rst -x ${ncdir}/water_soft_4.nc  -ref ${rstdir}/water_dens_3.rst

# Step 6: Second minimization
# Additional minimization after density equilibration
pmemd.cuda -O -i 5_2nd_Minimization.in                -p $prmtop -o ${outdir}/min_5.out        -c ${rstdir}/water_soft_4.rst -r ${rstdir}/min_5.rst        -x ${ncdir}/min_5.nc         -ref ${rstdir}/water_soft_4.rst

# Step 7: First equilibration (10 kcal/mol/Å²)
# Equilibrates system with strong positional restraints
pmemd.cuda -O -i 6_First_equil_10kcalmolA2.in         -p $prmtop -o ${outdir}/first_equ_6.out  -c ${rstdir}/min_5.rst        -r ${rstdir}/first_equ_6.rst  -x ${ncdir}/first_equ_6.nc   -ref ${rstdir}/min_5.rst

# Step 8: Second equilibration (1 kcal/mol/Å²)
# Equilibrates with moderate restraints
pmemd.cuda -O -i 7_Second_equilibration_1kcalmolA2.in -p $prmtop -o ${outdir}/second_equ_7.out -c ${rstdir}/first_equ_6.rst  -r ${rstdir}/second_equ_7.rst -x ${ncdir}/second_equ_7.nc  -ref ${rstdir}/first_equ_6.rst

# Step 9: Third equilibration (0.1 kcal/mol/Å²)
# Equilibrates with weak restraints
pmemd.cuda -O -i 8_Third_equilibration_01kcalmolA2.in -p $prmtop -o ${outdir}/third_equ_8.out  -c ${rstdir}/second_equ_7.rst -r ${rstdir}/third_equ_8.rst  -x ${ncdir}/third_equ_8.nc   -ref ${rstdir}/second_equ_7.rst

# Step 10: Unrestrained equilibration
# Final equilibration without restraints
pmemd.cuda -O -i 9_Unrestr_equilibration.in           -p $prmtop -o ${outdir}/unres_equ_9.out  -c ${rstdir}/third_equ_8.rst  -r ${rstdir}/unres_equ_9.rst  -x ${ncdir}/unres_equ_9.nc   -ref ${rstdir}/third_equ_8.rst

# Production MD simulations at different temperatures
# HMR stands for "Hydrogen Mass Repartitioning"

# Production MD at 310K (37°C - physiological temperature)
pmemd.cuda -O -i mdHMR310.in                          -p $prmtop -o ${outdir}/mdProd310.out    -c ${rstdir}/unres_equ_9.rst  -r ${rstdir}/mdProd310.rst    -x ${ncdir}/mdProd310.nc

# Production MD at 330K (57°C)
pmemd.cuda -O -i mdHMR330.in                          -p $prmtop -o ${outdir}/mdProd330.out    -c ${rstdir}/mdProd310.rst    -r ${rstdir}/mdProd330.rst    -x ${ncdir}/mdProd330.nc

# Production MD at 360K (87°C)
pmemd.cuda -O -i mdHMR360.in                          -p $prmtop -o ${outdir}/mdProd360.out    -c ${rstdir}/mdProd330.rst    -r ${rstdir}/mdProd360.rst    -x ${ncdir}/mdProd360.nc

# Production MD at 390K (117°C)
pmemd.cuda -O -i mdHMR390.in                          -p $prmtop -o ${outdir}/mdProd390.out    -c ${rstdir}/mdProd360.rst    -r ${rstdir}/mdProd390.rst    -x ${ncdir}/mdProd390.nc

# Load Python environment for trajectory analysis
module load python/3.12.3/gcc13.2/miniconda_24.4.0

# Compact trajectories and calculate RMSD (Mean Square Displacement)
# This removes water molecules and calculates RMSD (Root Mean Square Deviation)
python3 compactTraj_CalcRmsd.py solvated_system.pdb --prmtop-in HMR_solvated_system.prmtop  --traj-in ${ncdir}/mdProd*.nc --prmtop-out ${analdir}/rep0_noWat_hMD.prmtop --traj-out ${analdir}/rep0_noWat_hMD.nc --reference nowat_system.rst7 --save ${analdir}/cpptraj_full.in --output ${analdir}/iRMSD_r0.dat --run

# Delete trajectory files to save disk space
echo 'Deleting trajectory files'
rm ${ncdir}/*nc
echo 'Trajectory files deleted'

# Generate RMSD plot
echo 'Generating plot'
module load python/3.12.3/gcc13.2/miniconda_24.4.0
python3 RMSD_plot.py
echo 'iRMSD plot generated'
