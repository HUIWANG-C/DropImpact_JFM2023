# Analysis of high-speed drop impact onto deep liquid pool (DropImpact_JFM2023) 
This repository contains codes related to paper [**Analysis of high-speed drop impact onto deep liquid pool**](https://doi.org/10.1017/jfm.2023.701) that has been published in the _Journal of Fluid Mechanics_ in 2023.   
## 📚 Citation
If you use this code or data, please cite the following:  

```bibtex
@article{wang2023analysis,  
  title={Analysis of high-speed drop impact onto deep liquid pool},  
  author={Wang, Hui and Liu, Shuo and Bayeul-Lain{\'e}, Annie-Claude and Murphy, David and Katz, Joseph and Coutier-Delgosha, Olivier},  
  journal={Journal of Fluid Mechanics},  
  volume={972},  
  pages={A31},  
  year={2023},  
  publisher={Cambridge University Press}  
}
```

## Basilisk
To run the simulation, please check and install the open-source finite-volume adaptive [Basilisk](http://basilisk.fr/) solver.  
Please note that the code may not compile or run correctly depending on the software version.

### 🛠️ For compiling the code in parallel:  
`CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 drop.c -o drop -L$BASILISK/gl -lglutils -lfb_tiny -lm`  

### ▶️ Run the code in parallel:  
`mpirun -np 4 ./drop`
