# Analysis of high-speed drop impact onto deep liquid pool (DropImpact_JFM2023) 
This repository contains codes related to paper [Analysis of high-speed drop impact onto deep liquid pool](https://doi.org/10.1017/jfm.2023.701) that has been published in the _Journal of Fluid Mechanics_ in 2023. 

# Basilisk
In order to run the simulation, please check and install the open-source finite-volume adaptive [Basilisk](http://basilisk.fr/) solver.  
Please be aware that the code may fail to compile or run depending on the software version.

For compiling the code:  
`CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 drop.c -o drop -L$BASILISK/gl -lglutils -lfb_tiny -lm`  

Run the code in parallel:  
`mpirun -np 4 ./drop`
