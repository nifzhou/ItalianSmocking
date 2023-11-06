# Italian Smocking



## Reference

This repository provides source codes for paper:

Ningfeng Zhou, [Jing Ren](https://ren-jing.com/), [Olga Sorkine-Hornung](https://igl.ethz.ch/people/sorkine). Computational Smocking through Fabric-Thread Interaction.

We formalize the Italian smocking pattern and develop a simple method to represent and deform the planar mesh into the smocked ones given the initial stitching lines.

More details can be found at:



## Implementation

The source code contains two stages:

1. 2D Simulation: estimate the expected position of critical points after smocking and prepare inputs for the downstream simulator. (detailed description in folder *sim2d*)
2. 3D Deformation: guide the simulator to deform the planar mesh into a plausible smocked one. (detailed description in folder *deform3d*)



## Our results

|     Pattern      |                             Rect                             |                            vwave                             |                            square                            |                            zcurve                            |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|    our result    | ![image-20231105234742728](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234742728.png) | ![image-20231105234754337](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234754337.png) | ![image-20231105234801212](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234801212.png) | ![image-20231105234807580](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234807580.png) |
| real fabrication | ![image-20231105235307655](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235307655.png) | ![image-20231105235312129](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235312129.png) | ![image-20231105235323022](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235323022.png) | ![image-20231105235412321](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235412321.png) |

|     Pattern      |                            zigzag                            |                          half sharp                          |                        sharp (layer2)                        |                            curve                             |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|    our result    | ![image-20231105235723349](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235723349.png) | ![image-20231105234852426](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234852426.png) | ![image-20231105235717205](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235717205.png) | ![image-20231105235811735](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105235811735.png) |
| real fabrication | ![image-20231105234834626](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234834626.png) | ![image-20231105234856118](C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234856118.png) | <img src="C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234903414.png" alt="image-20231105234903414" style="zoom:80%;" /> | <img src="C:\Users\60242\AppData\Roaming\Typora\typora-user-images\image-20231105234949547.png" alt="image-20231105234949547" style="zoom:150%;" /> |



## Notes

TBD

