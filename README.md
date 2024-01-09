# Italian Smocking



## Reference

This repository provides source codes for paper:

[Ningfeng Zhou](https://github.com/nifzhou), [Jing Ren](https://ren-jing.com/), [Olga Sorkine-Hornung](https://igl.ethz.ch/people/sorkine). Computational Smocking through Fabric-Thread Interaction.

We formalize the Italian smocking pattern and develop a simple method to represent and deform the planar mesh into the smocked ones given the initial stitching lines.

More details can be found at: [project page], [paper], [[suppl. video]](https://youtu.be/L6AdmSCmbFc).



## Implementation

The source code contains two stages:

1. 2D Simulation: estimate the expected position of critical points after smocking and prepare inputs for the downstream simulator. (detailed description in folder *sim2d*)
2. 3D Deformation: guide the simulator to deform the planar mesh into a plausible smocked one. (detailed description in folder *deform3d*)



## Our results

|     Pattern      |                             Rect                             |                            vwave                             |                            square                            |                            zcurve                            |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|   input stitch   |                  ![rect](./assets/rect.png)                  |                 ![vwave](./assets/vwave.png)                 |                ![square](./assets/square.png)                |                ![zcurve](./assets/zcurve.png)                |
|    our result    | ![image-20231105234742728](./assets/image-20231105234742728.png) | ![image-20231105234754337](./assets/image-20231105234754337.png) | ![image-20231105234801212](./assets/image-20231105234801212.png) | ![image-20231105234807580](./assets/image-20231105234807580.png) |
| real fabrication | ![image-20231105235307655](./assets/image-20231105235307655.png) | ![image-20231105235312129](./assets/image-20231105235312129.png) | ![image-20231105235323022](./assets/image-20231105235323022.png) | ![image-20231105235412321](./assets/image-20231105235412321.png) |

|        Pattern         |                            zigzag                            |                          half sharp                          |                        sharp (layer2)                        |                            curve                             |
| :--------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| input stitch (partial) |           ![zigzagLike0](./assets/zigzagLike0.png)           |           ![zigzagLike1](./assets/zigzagLike1.png)           |           ![zigzagLike2](./assets/zigzagLike2.png)           |                 ![curve](./assets/curve.png)                 |
|       our result       | ![image-20231105235723349](./assets/image-20231105235723349.png) | ![image-20231105234852426](./assets/image-20231105234852426.png) | ![image-20231105235717205](./assets/image-20231105235717205.png) | ![image-20231105235811735](./assets/image-20231105235811735.png) |
|    real fabrication    | ![image-20231105234834626](./assets/image-20231105234834626.png) | ![image-20231105234856118](./assets/image-20231105234856118.png) | <img src="./assets/image-20231105234903414.png" alt="image-20231105234903414" style="zoom:80%;" /> | <img src="./assets/image-20231105234949547.png" alt="image-20231105234949547" style="zoom:150%;" /> |



## Notes

### Acknowledgements

We would like to thank the anonymous reviewers for their insightful feedback. We extend our gratitude to M. Rifad (YouTube channel ``DIY Stitching``), F. Shanas (YouTube channel ``handiworks``), and S. Fyms (YouTube channel ``FymsEmbroidery``) for generously granting us permission to use the images of their remarkable fabrication results.
This work was supported in part by the ERC Consolidator Grant No.101003104 (MYCLOTH).

