# Source Code for 3D Mesh Deformation

## Reference

This repository provides source code for: 

[our paper citation]

The code is mainly modified based on the Codimensional Incremental Potentail Contact simulator (C-IPC). The reference to the original simulator is:

Minchen Li, Danny M. Kaufman, Chenfanfu Jiang, [Codimensional Incremental Potential Contact](https://ipc-sim.github.io/C-IPC/), ACM Transactions on Graphics (SIGGRAPH  2021)

All the modified codes are highlighted in the code between:

```
// Smock Added Start

*replaced or added code in Cpp*

// Smock Added End
```

Or:

```
## Smock Added Start

*replaced or added code in Python*

## Smock Added End
```



## Installation (same as original C-IPC)

### UBUNTU
Install python
```
sudo apt install python3-distutils python3-dev python3-pip python3-pybind11 zlib1g-dev libboost-all-dev libeigen3-dev freeglut3-dev libgmp3-dev
pip3 install pybind11
```
Build library
```
python build.py
```

### MacOS
Build library
```
python build_Mac.py
```



## Run the scripts

### input
1. Mesh .obj with stitching information (represented by "l vIdx1 vIdx2"; "l vIdx2 vIdx3" …)

   Stitching Information Format:
    ```
    l vIdx1 vIdx2 # l: stitching line vIdx: vertex index for the stitching line
    ```

2. Positional constraints information file (xx.constraintInfo)

   Constraints Format:

   ```
   v vIdx vX vY vZ # position for the vIdx-th vertices
   s sIdx len 		# length for the sIdx-th stitching lines
   ```
   
   
### script

To run the examples in the paper:

```
cd Projects/FEMShell
python xx.py
```



## Reproduce the examples

To avoid redundant files, necessary input files and scripts to reproduce our results are prepared in the C-IPC.



All the scripts are provided in *Projects/FEMShell/*.



## Note

The code only tested in Ubuntu, and the constraint information are only used with cloth using isotropic membrane model without strain limiting. The force analysis under other settings remain unchanged as the original C-IPC.