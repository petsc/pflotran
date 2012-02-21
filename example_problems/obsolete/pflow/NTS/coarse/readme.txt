The input files in this directory are for running a simulation of flow 
around a melt glass resulting from an underground test at the Nevada 
Test Site.  Two "pflow.in" files are required.  The first is necessary 
to find the steady-state solution WITHOUT the melt glass in place.  To 
find this, do (assuming that the working directory contains the 'pflow' 
executable or a link thereto):

mpirun -np N ./pflow -use_thc -use_analytical -ksp_type bcgs -pflowin pflow1.in

The -use_thc option tells PFLOW to use the "THC" system: non-isothermal, 
fully saturated flow with a nonreactive tracer.  The -use_analytical 
option specifies that an analytically computed Jacobian should be used.
Because it takes a long time to reach steady state, it is probably useful 
to add something like '-chkptfreq 50' to write a checkpoint file every 50 
flow steps.

When a steady state is reached, pflow will terminate and will write the 
state to a file name 'pflow_init0.dat'.  Copy or rename the final 
pflow_init0.dat file to pflow_init.dat, and then add the following at the 
bottom of pflow_init.dat, before the final line containing only a '/' 
character:

:glass
  11  31  16  35  41  45   1.d7  75.   1.   1.

The '75.' above specifies the initial temperature (degrees Celcius) of 
the melt glass after re-wetting.  Since this temperature is poorly 
constrained, you may want to experiment with varying this temperature.

With pflow_init.dat edited, you can run pflow with the melt glass in place:

mpirun -np N ./pflow -use_thc -use_analytical -ksp_type bcgs -pflowin pflow2.in
