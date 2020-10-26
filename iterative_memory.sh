#bin/bash

#
#  Program to iterativly reconstruct the memory of a system
#
#  Necessary input:
# 	- Lammps executable with "fix gle" command
# 	- Results for the force-autocorrelation function
# 

lammps_bin=/home/gerhard/Documents/bin/lammps-master/src/lmp_mpi

# iterative parameter
run=1_reconstruct_005
start=1
iterations=125

# regularization
correction_length=16; 		#how many steps the correction should be applied
incr_step=4; 			#minimal increment between each iteration
incr_N=1;

# simulation
Nmem=500; 			#length of memory
dtorig=0.005; 			#timestep
Torig=1.0; 			#temperature
morig=80;			#mass

# seed
seed=12345

# variables
let incr=-1

# functions
function create_input { #creates a LAMMPS-input script
  let freq=4*Nmem
  echo "
    ##############################
    # set input values
    ##############################

    # box
    variable		box_size equal 50

    # system paramter
    variable    	T equal ${Torig}
    variable		Nc equal 100
    variable 		M equal ${morig}

    # tuneable parameters
    variable		skin equal 1.5

    # interactions
    variable		rc equal 1.112

    # simulation
    variable	dt equal ${dtorig}
    variable	eq_steps equal 5000
    variable	sim_steps equal 100000

    ##############################
    # system setup
    ##############################

    # lj 3d solvent
    units		lj
    dimension	3
    atom_style	atomic

    # create lattice (to place particles), simulation box
    lattice         sc 1.0
    region          simbox block 0 "'${box_size}'" 0 "'${box_size}'" 0 "'${box_size}'"
    boundary 	p p p
    create_box      1 simbox

    # create particles
    # colloid
    create_atoms    1 random "'${Nc}'" ${seed} simbox
    group		colloid type 1

    # initialize particles
    mass            1 "'${M}'"
    timestep	"'${dt}'"

    # interaction
    pair_style      none

    # define nighbourlist
    neighbor	0.3 nsq
    neigh_modify	once yes
    atom_modify	sort 1000 1.0

    ##############################
    # define integration
    ##############################

    # assign random velcities velocities
    velocity    all create "'$T'" ${seed}

    # define integration
    fix	        2 all gle "'$T'" $1 ${Nmem} ${seed}

    ##############################
    # eq run
    ##############################

    thermo      1000
    run	        "'${eq_steps}'"

    ##############################
    # define on-the-fly calculations
    ##############################

    # self and pair time correlation function of velocities and forces
    variable	fx atom fx
    variable	fy atom fy
    variable	fz atom fz
    variable	vx atom vx
    variable	vy atom vy
    variable	vz atom vz
    fix 	3 colloid ave/correlate/peratom 1 ${Nmem} ${freq} v_vx v_vy v_vz v_fx v_fy v_fz file $2 ave running overwrite switch pergroup restart

    ##############################
    # production run
    ##############################

    thermo          1000
    run	        "'${sim_steps}'"
  "
}

function create_next { # creates an awk-script to calculate memory for the next step
  echo "
   function abs(value){return (value<0?-value:value);}
   BEGIN{
    Nmem = ${Nmem};
    dtorig = ${dtorig};
    Torig = ${Torig};
    m = ${morig}
    
    correction_length = ${correction_length};
    incr = ${incr}; 
  
    counter = -1
   }
   {
    if (FILENAME == "'"old_memory.cor"'") data_mem[FNR-1] = "'$3/dtorig/dtorig'";
    if (FILENAME == "'"old_result.cor"'") data_result_gle[FNR-5] = ("'$4'"+"'$6'"+"'$8'")/3;
    if (FILENAME == "'"md_result.cor"'") data_result_md[FNR-1] = "'$2'";
   }
   END{

    for (t=0; t<Nmem; t++) {
      data_diff[t] = - data_result_gle[t] + data_result_md[t];
    }
    # differentiation
    for (t=0; t<Nmem; t++) {
      left = data_diff[t-1];
      if (left == 0) left = data_diff[t]
      mid = data_diff[t]
      right = data_diff[t+1]
      if (right == 0) right = data_diff[t]
    
      data_cor[t]=-m*m*(left-2*mid+right)/dtorig/dtorig
    }
    
    for (t=0; t<Nmem; t++) {
      if ( counter == -1 && t >= incr) {
	counter = 0;
      }
      if (counter >= 0 && counter < correction_length) {
	data_cor[t] *= (correction_length - counter)/correction_length;
	counter++;
      }
      if (counter>=correction_length) data_cor[t] = 0;
    }
    
    for (t=0; t<Nmem; t++) {
      print t*dtorig, data_mem[t]+data_cor[t];
    }
    
   }
  "
}

# if start=1
if [ $start -eq 1 ]; then 
  # if folder exists
  if [ -d run$run ]; then
    read -p "run$run exists! Are you sure to overwrite it? " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      rm -r run$run
    else
      exit 1;
    fi
  fi
  mkdir run$run
fi

if [ $start -gt 1 ] && [ ! -d run$run ]; then 
  exit 1;
fi

cd run$run

if [ $start -eq 1 ]; then 
  # create first folder
  mkdir gle1
  cd gle1
  create_input gle1_memory.cor gle1_result.cor > in.colloid.gle
  let seed=seed
  cp ../../md_result.cor .
  cp ../../md_result_ff.cor .
  
  # start the iteration
  awk '{print $1,$2/'"${Torig}}" md_result_ff.cor > gle1_memory.cor
  echo "Start step 1: Incr=${incr}"
  ${lammps_bin} -in in.colloid.gle
  cd ..
else
  cd gle$start
  create_input gle${start}_memory.cor gle${start}_result.cor > in.colloid.gle
  ${lammps_bin} -in in.colloid.gle
  cd ..
fi

let start=start+1

for ((i=$start; i<=$iterations; i=i+1))
do
  let oldstep=$i-1
  if [ -d gle$i ]; then
    rm -r gle$i
  fi
  mkdir gle$i
  cd gle$i
  create_input gle${i}_memory.cor gle${i}_result.cor > in.colloid.gle
  let seed=seed
  create_next gle${oldstep}_result.cor > next_step.awk
  cp ../../md_result.cor .
  cp ../gle${oldstep}/gle${oldstep}_result.cor ./old_result.cor
  cp ../gle${oldstep}/ansatz.dat ./old_memory.cor
  awk -f next_step.awk md_result.cor old_result.cor old_memory.cor > gle${i}_memory.cor
  if (( $i % $incr_N == 0 )) 
  then
    let incr=$incr+$incr_step
  fi
  echo "Start step $i: Incr=${incr}"
  ${lammps_bin} -in in.colloid.gle
  cd ..
done

# create scripts to plot results
# v(t)v
printf "reset\nplot " > plot_vv.g
for ((i=1; i<=$iterations; i=i+1))
do
  printf "'gle$i/gle${i}_result.cor' u 2:(\$4+\$6+\$8)," >> plot_vv.g
done

# f(t)f
printf "reset\nplot " > plot_ff.g
for ((i=1; i<=$iterations; i=i+1))
do
  printf "'gle$i/gle${i}_result.cor' u 2:(\$10+\$12+\$14)," >> plot_ff.g
done
