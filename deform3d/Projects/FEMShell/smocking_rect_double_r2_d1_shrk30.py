import sys
sys.path.insert(0, "../../Python")
import Drivers
from JGSL import *
import time

if __name__ == "__main__":
    sim = Drivers.FEMDiscreteShellBase("double", 3)
    clothI = 0
    membEMult = 1
    bendEMult = 100
    sim.frame_num = 100

    time_init = time.time()

    # load garment rest shape
    weight_pos = 0.01
    sim.add_garment_3D("input/final/rect/double_row2_d1_2c.obj", Vector3d(0, 0, 0), Vector3d(1, 1, 1), \
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    sim.initialize_garment()

    ## if comment add_constraint_file(): org C-IPC
    sim.add_constraint_file("input/final/rect/double_row2n5_d1_refineH_shrk30.constraintInfo", weight_pos)

    ## uncomment this to eliminate non-zero rest length of stitching lines
    # sim.stitch_rstLen = StdVectorXd()      # vector 1-1 correspondence (per stitching line)

    ## assign 0 to k_stitch to give up stitching constraints (act as there is no stitching lines)
    sim.k_stitch = 0.1      # default: 0.1

    ## uncomment this to eliminate positional constraints
    # sim.weight_tgt_pos = 0.0               # weight for constraint energy
    # sim.target_pos = StdVectorVector3d()      
    # sim.target_pos_idx = StdVectorXi()

    sim.gravity = Vector3d(0, 0, 0)
    sim.dt = 10
    sim.frame_dt = 0.04
    sim.withCollision = True
    sim.mu = 0

    # iso, no strain limit
    sim.initialize(sim.cloth_density_iso[clothI], membEMult,
        sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0)
    sim.bendingStiffMult = bendEMult / membEMult
    sim.kappa_s = Vector2d(0, 0)


    # load draped garment as initial config
    # sim.load_frame("input/final/rect/double_row2n5_d1_def2c.obj")

    sim.initialize_OIPC(1e-3, 0)

    time_sim = time.time()

    sim.run()

    time_end = time.time()

    print("total time: %.4fs, init time: %.4fs, sim time: %.4fs" %(time_end-time_init, time_sim-time_init, time_end-time_sim))
