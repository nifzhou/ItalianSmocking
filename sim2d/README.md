# Source code for 2D Pattern Simulation

## Requirement

Python version >=3.10

pip version >=20.3

necessary packages: numpy, scipy, matplotlib



## Run the scripts

The main structure for the 2D simulator is shown in "main_sim.py". The hyper-parameters are shown in the main file with detailed descriptions.

To run some samples, the main file requires <u>two input file paths</u>. The first file path refers to the mesh with complete stitching lines. The second file path refers to the same mesh with only front stitching lines. A sample script to run it is:

```
python .\main_sim.py ..\data\curve\row3_d1.obj ..\data\curve\row3_d1_half.obj
```



All the mesh files shown in the paper are prepared in the *data* folder. 



## Output

To support the downstream C-IPC, the 2D deformer will output the necessary files in the original folder of the first input mesh:

1. original mesh for C-IPC (suffix: **"_2c"**): rescale the input mesh and change the OBJ format into C-IPC standard
2. original mesh with configuration (suffix: "**_def**"): original mesh with positive(negative) height for midpoints [optional]
3. original mesh with configuration for C-IPC (suffix: "**_def2c**"): rescaled and revised mesh with positive(negative) height for midpoints to feed into the C-IPC
4. constraints information to feed into C-IPC (suffix: "**.constraintInfo**")



The output files can be pasted to the folder that aligns with the input path of the downstream C-IPC. To avoid redundant files, the necessary files to reproduce the C-IPC results in the paper are provided in *deform3d* folder (details see the README in it)



## Note

The code only tested in Ubuntu and windows 11.







