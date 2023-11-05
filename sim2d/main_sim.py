
import sys
import numpy as np
from smockPattern import EngSmockPattern
from patchSpringNetwork import MassSpringNetwork
import time
from utils import write_obj_file, read_obj_file, generate_eng_smock_from_all_and_front_obj_file, write_constraints, write_cipc_input_mesh


## python version 3.10(.11)
## pip version >= 20.3
## numpy matplotlib

## sys.argv:
## [1] obj file to be opened (full smocked)
## [2] obj file to be opened (front half smocked)



#########################################################
######################## Main ###########################
#########################################################
if __name__ == "__main__":

    '''hyper parameters start'''
    shrink_percentage = 0.2     # [0,1] control the percentage of the pulled-out stitching line (roughly)
    spring_stiffness = 1.0
    edge_stiffness = 5.0
    spring_rest_length = 0.01   # >=epsilon in spring network
    coeff_damp = 0.9
    dt = 0.1
    sim_num_iter = 200000        # upper bound for number of simulation iteration
    sim_anime = False   # 2D animation
    '''hyper parameters end'''

    '''hyper parameters for arap 3d deformer ()'''
    
    ## False(default): output Constraints for downstream C-IPC
    ## True: Apply ARAP to preview deformed mesh (no constraints output)
    useArapPreview = False      
    arap_weight = 1.0
    arap_num_iter = 100
    

    time_init = time.time()

    if (len(sys.argv) == 2):
        V, F, StitchLines = read_obj_file(sys.argv[1]) 
        smock = EngSmockPattern(V, F, StitchLines)
    elif(len(sys.argv) == 3):
        smock = generate_eng_smock_from_all_and_front_obj_file(sys.argv[1], sys.argv[2], arap_weight, arap_num_iter)


    ## generate regular network for force analysis
    sim_spring = MassSpringNetwork(smock, shrink_percentage, spring_stiffness, edge_stiffness, spring_rest_length, coeff_damp, dt, sim_num_iter)

    rest_edges, rest_stitches = sim_spring.pre_calculate_net()
    smock.stitchLenInit = rest_stitches

    time_sim = time.time()


    underlay_pos, idx_graph2mesh, idx_mesh2graph, num_iter = sim_spring.run(rest_edges, rest_stitches, withAnime=sim_anime)
        
    time_deform = time.time()

    if useArapPreview:
        deformed_V = smock.deform(underlay_pos, idx_graph2mesh, idx_mesh2graph)
        print("------ARAP deformation done!------")
    else:
        constraints = smock.extract_constraints(underlay_pos, idx_graph2mesh, idx_mesh2graph)
        write_constraints(smock, sim_spring, constraints, sys.argv[1][:-4]+'_refineH_shrk'+str(int(shrink_percentage*100))+'.constraintInfo')
        print("------Constraints Output done!------")

        ## prepare obj input for C-IPC (scaled and revise format)
        write_cipc_input_mesh(smock, sys.argv[1])


    time_end = time.time()
    print("total 2D sim time: %.4fs" %(time_end-time_init))


