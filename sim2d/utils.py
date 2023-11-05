from smockPattern import EngSmockPattern
from patchSpringNetwork import MassSpringNetwork
import numpy as np

#########################################################
##################### IO Functions#######################
#########################################################

def write_obj_file(V, F, StitchLines=None, file_path="./deformed.obj"):
    with open(file_path, 'w') as f:
        f.write("# OBJ file\n")

        # write vertices
        for id_v in range(len(V)):
            f.write("v %.4f %.4f %.4f\n" % (V[id_v, 0], V[id_v, 1], V[id_v, 2]))

        # write faces
        for id_f in range(len(F)):
            f.write("f")
            p = F[id_f]
            for id_v in range(len(p)):
                f.write(" %d" % (p[id_v] + 1))
            f.write("\n")
        
        if StitchLines is not None:
            for id_l in range(len(StitchLines)):
                f.write("l")
                stitch = StitchLines[id_l]
                for id_v in range(len(stitch)):
                    f.write(" %d" % (stitch[id_v] + 1))
                f.write("\n")


def read_obj_file(file_path):
    '''
    direct read obj file (w/o front & back stitch line info)
    '''

    V, F, StitchLines = [], [], []

    with open(file_path, 'r') as f:
        lines = f.readlines()

        for line in lines:
            elems = line.split()

            if elems[0] == "v":
                pos = []
                for i in range(1, len(elems)):
                    pos.append(float(elems[i]))    # v x y z
                V.append(pos)

            elif elems[0] == "f":
                face = []
                for j in range(1, len(elems)):
                    vf = elems[j].split('/')[0]    # f vtxIdx/texIdx/normalIdx
                    face.append(int(vf) - 1)
                F.append(np.array(face))

            elif elems[0] == "l":
                stitch = []
                for k in range(1, len(elems)):
                    stitch.append(int(elems[k]) - 1)
                StitchLines.append(stitch)         # l v1 v2 v3 ...
            

    V = np.array(V)
    return V, F, StitchLines


def write_constraints(smockPattern:EngSmockPattern, sim_spring:MassSpringNetwork, constraints:dict, file_path):
    '''
    '' write stitch rest length (follow stitch order) & positional constraints of underlay&mid-stitch nodes
    '' format:
    '' s idx_stitch rest_len
    '' v idx_v x y z
    '''
    stitch_vec = sim_spring.pos[sim_spring.stitchStartIdx] - sim_spring.pos[sim_spring.stitchEndIdx]
    stitch = np.linalg.norm(stitch_vec, axis=1, keepdims=False)

    # write original obj to cipc standard (can rescale ==> set "scale")
    # C-IPC works better with the largest dimension around 1.0
    # rescale
    bbx_net = np.max(smockPattern.org_V, axis=0) - np.min(smockPattern.org_V, axis=0)
    scale = 1.0/max(bbx_net)    # align with C-IPC mesh resize
    print("constraint scale: %f \n" %(scale))

    with open(file_path, 'w') as f:
        f.write("# constraint info \n")

        for idx_s in range(len(stitch)):
            f.write("s %d %.6f\n" % (idx_s, scale * stitch[idx_s]))

        for idx_v, pos_v in constraints.items():
            f.write("v %d %.4f %.4f %.4f\n" % (idx_v, scale * pos_v[0], scale * pos_v[1], scale * pos_v[2]))
    
    print("write constraint info to '%s'" %(file_path))


def write_cipc_input_mesh(smock:EngSmockPattern, input_file_path:str):
    
    # write original obj to cipc standard (can rescale ==> set "scale")
    # C-IPC works better with the largest dimension around 1.0
    # rescale
    bbx_net = np.max(smock.org_V, axis=0) - np.min(smock.org_V, axis=0)
    scale = 1.0/max(bbx_net)

    cipc_file_path = input_file_path[:-4]+'_2c.obj'

    file_org = open(input_file_path, 'r')
    lines_org = file_org.readlines()

    with open(cipc_file_path, 'w') as f:

        for line in lines_org:
            elems = line.split()
            if elems:
                if elems[0] == "v":
                    f.write('v %f %f %f\n' %(float(elems[1]) * scale, float(elems[2]) * scale, float(elems[3]) * scale))
                elif elems[0] == "l":
                    f.write('stitch %d %d %d 0\n' %(int(elems[1]), int(elems[2]), int(elems[2])))
                else:
                    f.write(line)
    file_org.close()
    
    print("saved file: '%s'" %(cipc_file_path))
    print("------original CIPC pattern generation done!------")

    sim_spring = MassSpringNetwork(smock)
    rest_edges, rest_stitches = sim_spring.pre_calculate_net()
    smock.stitchLenInit = rest_stitches
    V, F, StitchLines = smock.generate_init_config()

    # write obj at original standard
    new_file_path = input_file_path[:-4]+'_def.obj'
    write_obj_file(V * scale, F, StitchLines, new_file_path)
    print("saved file: '%s'" %(new_file_path))
    print("scaling factor: %f" %(scale))
    print("------shift of pattern generation done!------")

    # write obj at CIPC convention (support only two vert at a line segment)
    file_input = open(new_file_path, 'r')
    lines_input = file_input.readlines()

    with open(input_file_path[:-4]+'_def2c.obj', 'w') as f:

        for line in lines_input:
            elems = line.split()
            if elems:
                if elems[0] == "l":
                    f.write('stitch %d %d %d 0\n' %(int(elems[1]), int(elems[2]), int(elems[2])))
                else:
                    f.write(line)


    file_input.close()
    print("saved file: '%s'" %(input_file_path[:-4]+'_def2c.obj'))
    print("------initial config of CIPC pattern generation done!------")




#########################################################
################### Pattern Generation ##################
#########################################################

def generate_eng_smock_from_all_and_front_obj_file(file_path_all, file_path_front, arap_weight=1.0, arap_num_iter=1000):
    '''
    merge read files suitable for english smocking (define front and back lines)

    the "l" in obj files contains only 2 vertices each (l v1 v2)
    '''
    V_all, F_all, StitchLines_all = read_obj_file(file_path_all)
    V_front, F_front, StitchLines_front = read_obj_file(file_path_front)

    assert len(V_all) == len(V_front) and len(F_all) == len(F_front), "two obj file did not have the same number of vertices and faces."

    front_line_idx = []

    for l_id in range(len(StitchLines_all)):
        line = StitchLines_all[l_id]
        line_reverse = [line[1], line[0]]
        if line in StitchLines_front:
            front_line_idx.append(l_id)
        elif line_reverse in StitchLines_front:
            front_line_idx.append(l_id)
    
    smock = EngSmockPattern(V_all, F_all, StitchLines_all, front_line_idx, arap_weight, arap_num_iter)

    print("------smocking pattern generation done!------")
    print("--Total Line Num: %d, Front Line Num: %d--" %(len(StitchLines_all), len(front_line_idx)))
    print("number of vertices in mesh: %d" %(len(V_all)))
    return smock
