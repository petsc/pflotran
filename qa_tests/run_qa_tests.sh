#!/bin/bash

nproc=1

# Notes:
#
# To print the %-error at each grid cell in the QA test simulation, set
# the argument "print_error=true," otherwise set it to "print_error=false".
# 
# To print a plot which visually compares the pflotran solution to the
# analytical solution, set the argument "plot_flag=true," otherwise set it
# to "plot_flag=false". If plots are generated, the QA test script will pause
# and will not re-start until you close the plot.
# 
# BC = boundary condition
# BC_1st_kind = dirichlet boundary condition type
# BC_2nd_kind = neumann boundary condition type

echo ' '
echo '========================================================================='
echo '====== Running STEADY THERMAL tests ====================================='
echo '========================================================================='
echo '====== TH Mode =========================================================='
echo '========================================================================='
cd thermal/th_mode/steady
pf=$(cd ../../../../src/pflotran; pwd)
echo $PWD
echo ' '

echo '============================='
echo '  1D_conduction_BC_1st_kind'
echo '============================='
cd 1D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_kind > screen.txt
python run_1D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  1D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 1D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_2nd_kind > screen.txt
python run_1D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '=============================='
echo '  2D_conduction_BC_1st_kind'
echo '=============================='
cd 2D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 2D_conduction_BC_1st_kind > screen.txt
python run_2D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  2D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 2D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 2D_conduction_BC_1st_2nd_kind > screen.txt
python run_2D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '=============================='
echo '  3D_conduction_BC_1st_kind'
echo '=============================='
cd 3D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 3D_conduction_BC_1st_kind > screen.txt
python run_3D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo ' '
echo '========================================================================='
echo '====== Running STEADY THERMAL tests ====================================='
echo '========================================================================='
echo '====== GENERAL Mode ====================================================='
echo '========================================================================='
cd ../../general_mode/steady
pf=$(cd ../../../../src/pflotran; pwd)
echo $PWD
echo ' '

echo '============================='
echo '  1D_conduction_BC_1st_kind'
echo '============================='
cd 1D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_kind > screen.txt
python run_1D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  1D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 1D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_2nd_kind > screen.txt
python run_1D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '=============================='
echo '  2D_conduction_BC_1st_kind'
echo '=============================='
cd 2D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 2D_conduction_BC_1st_kind > screen.txt
python run_2D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  2D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 2D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 2D_conduction_BC_1st_2nd_kind > screen.txt
python run_2D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '=============================='
echo '  3D_conduction_BC_1st_kind'
echo '=============================='
cd 3D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 3D_conduction_BC_1st_kind > screen.txt
python run_3D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo ' '
echo '========================================================================='
echo '====== Running TRANSIENT THERMAL tests =================================='
echo '========================================================================='
echo '====== TH Mode =========================================================='
echo '========================================================================='
cd ../../th_mode/transient
pf=$(cd ../../../../src/pflotran; pwd)
echo $PWD
echo ' '

echo '============================='
echo '  1D_conduction_BC_1st_kind'
echo '============================='
cd 1D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_kind > screen.txt
python run_1D_conduction_BC_1st_kind.py print_error=false plot_flag=false
#rm *.vtk *.out screen.txt
cd ..

echo '============================='
echo '  1D_conduction_BC_2nd_kind'
echo '============================='
cd 1D_conduction_BC_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_2nd_kind > screen.txt
python run_1D_conduction_BC_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  1D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 1D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_2nd_kind > screen.txt
python run_1D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  2D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 2D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 2D_conduction_BC_1st_2nd_kind > screen.txt
python run_2D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo ' '
echo '========================================================================='
echo '====== Running TRANSIENT THERMAL tests =================================='
echo '========================================================================='
echo '====== GENERAL Mode ====================================================='
echo '========================================================================='
cd ../../general_mode/transient
pf=$(cd ../../../../src/pflotran; pwd)
echo $PWD
echo ' '

echo '============================='
echo '  1D_conduction_BC_1st_kind'
echo '============================='
cd 1D_conduction_BC_1st_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_kind > screen.txt
python run_1D_conduction_BC_1st_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '============================='
echo '  1D_conduction_BC_2nd_kind'
echo '============================='
cd 1D_conduction_BC_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_2nd_kind > screen.txt
python run_1D_conduction_BC_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  1D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 1D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 1D_conduction_BC_1st_2nd_kind > screen.txt
python run_1D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo '================================='
echo '  2D_conduction_BC_1st_2nd_kind'
echo '================================='
cd 2D_conduction_BC_1st_2nd_kind
echo 'Running PFLOTRAN simulation . . .'
mpirun -np $nproc $pf/pflotran -input_prefix 2D_conduction_BC_1st_2nd_kind > screen.txt
python run_2D_conduction_BC_1st_2nd_kind.py print_error=false plot_flag=false
rm *.vtk *.out screen.txt
cd ..

echo ' '
echo '========================================================================='
echo '====== Running FLOW tests ==============================================='
echo '========================================================================='
cd ../../../flow
pf=$(cd ../../src/pflotran; pwd)
echo $PWD
echo ' '
echo 'No tests to run. . .'

echo ' '
echo ' DONE.'