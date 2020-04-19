Usage:

* Set up shell environment with preconfigured compiler, MPI, etc. (other systems will be added to the script later) 

. environment.sh idun

* Compile everything (batch job does this over again for now, but why not)

make

* Run on 8 EPIC nodes x 18 ranks x 2 OMP-threads

qsub dambreak_idun.pbs

* Animate output

make dambreak.mp4
