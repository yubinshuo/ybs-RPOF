[README.md](https://github.com/yubinshuo/ybs-RPOF/files/7098097/README.md)
#  ybs-RPOF
## Mixed ReaxFF-Parameters Optimization Framework for reaction force field through SA, GA, PAO, SOPPI etc.
In this repository:
1. We provide python scripts for the optimization of reaction force field parameters.
2. We provide a example for the optimization of reaction force field parameters.

## Requirements 
* Python 3.6.9
* Scipy 1.4.1
* Pandas 1.0.2 
* Numpy 1.18.1
* matplotlib 3.3.4
* lammps.so 2018.Aug.22


## example Data
1. The data for example are available at: ```./example/```


## Run the script
#### Example:
```bash
python ybs.py -np 10 -mnp 4 -pso_new
```

### output files  
The output files are available at: ```./outdata/```
Final Optimization model of ReaxFF named as ffield.reax.start under ```./outdata/``` folder

Note: The input file in the working folder contains four files in the format shown in the example.
