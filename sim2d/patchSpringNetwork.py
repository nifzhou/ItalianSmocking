import numpy as np
from smockPattern import EngSmockPattern
import copy
from arap import ARAP

import matplotlib.pyplot as plt

thred_len_stop = 0.0      # threshold of stop criterion (0.1 for smocking pattern, 0.01 for simple stitch)

def gen_network(smockPattern:EngSmockPattern):
    '''
    for attempt 5: use sparse spring network for analysis
    generate regular sparse network for inextensibility constraints
    sparse network consists:
        1. stitching line: small rest length (adaptive?)
        2. fabric connection: stiff for stretch; no resistance for compression

    extra nodes are only used for the constraints of the stitching line, no additional constraint position is added to the ARAP afterwards (uncertain y)
    
    initally try simple shrink both node of the stretched spring. (no matter it is stitching node or free node)

    output:
        pos: position of the network nodes
        idx_graph2mesh (list): idx_graph2mesh[graph_idx] = vtx_idx
        stitchStartIdx, stitchEndIdx: start & end node idx of the stitching lines
        edgeStartIdx, edgeEndIdx: start & end node idx of the regular network graph (with possible stitching line pair)
        len_all_stitch: sum of all the initial length of the stitching lines

    '''

    # extract stitching nodes from smocking pattern
    idx_vts = []
    for lines in smockPattern.StitchLines:
        for vts in lines:
            idx_vts.append(vts)
    idx_stitch = list(set(idx_vts))

    print("number of stitching nodes: %d" %(len(idx_stitch)))

    pos_stitch = smockPattern.V[idx_stitch]

    # extract the stitch springs from the smocking pattern
    stitchStartIdx_mesh = []
    stitchEndIdx_mesh = []
    for lines in smockPattern.StitchLines:
        for vts in range(len(lines)-1):    # adjacent vert pairs in stitching line form a spring
            stitchStartIdx_mesh.append(lines[vts])
            stitchEndIdx_mesh.append(lines[vts+1])

    '''
    Generate spring net
    Stitching lines -> stitch spring
    original fabric -> edge spring (sparsely distributed if fabric exists there)
    '''
    # initial: find the shortest initial stitching line length as the unit length of the network (also rest length for - & | edges)
    # but: generated stitched pattern did not preserve the same length as the unsubdivided smocking pattern
    # thus: change equal-length net to "stitch-node-aligned" net (generate net with same x & z value of the stitch nodes, interpolate if gap is over threshold)
    diff_L_vec = smockPattern.V[stitchStartIdx_mesh] - smockPattern.V[stitchEndIdx_mesh]
    diff_L = np.linalg.norm(diff_L_vec, axis=1, keepdims=True)
    len_all_stitch = np.sum(diff_L)

    dL = np.min(np.max(np.abs(diff_L_vec), axis=1))        # grid size is the smallest dx or dz     

    # extract all x & z value of stitch nodes (sorted)
    x_net = np.sort(np.unique(pos_stitch[:,0]))
    z_net = np.sort(np.unique(pos_stitch[:,2]))

    ######################################
    #########Add nodes inside pattern#####
    ######################################
    # interpolate if dist above threshold
    diff_x = x_net[1:] - x_net[:-1]
    add_x_mask = diff_x > (dL * 1.5)   # dL * scale: threshold of interpolation
    # add in inverse order to prevent index confusion
    for add_idx in reversed(range(len(add_x_mask))):
        if add_x_mask[add_idx] == True:
            num = round(diff_x[add_idx] / dL)
            step = diff_x[add_idx] / num
            add_x = np.arange(x_net[add_idx], x_net[add_idx+1], step)
            x_net = np.insert(x_net, add_idx+1, add_x[1:])

    diff_z = z_net[1:] - z_net[:-1]
    add_z_mask = diff_z > (dL * 1.5)
    for add_idx in reversed(range(len(add_z_mask))):
        if add_z_mask[add_idx] == True:
            num = round(diff_z[add_idx] / dL)
            step = diff_z[add_idx] / num
            add_z = np.arange(z_net[add_idx], z_net[add_idx+1], step)
            z_net = np.insert(z_net, add_idx+1, add_z[1:])


    ######################################
    ########Add nodes around boundary#####
    ######################################
    
    # find the bounding box of the network
    bbx_net = [np.min(smockPattern.V, axis=0), np.max(smockPattern.V, axis=0)]

    #interpolate the boundary nodes to net
    num_x_low = int((x_net[0] - bbx_net[0][0]) / dL) + 1
    num_x_high = int((bbx_net[1][0] - x_net[-1]) / dL) + 1
    num_y_low = int((z_net[0] - bbx_net[0][2]) / dL) + 1
    num_y_high = int((bbx_net[1][2] - z_net[-1]) / dL) + 1

    # np.arange auto omit the last one if can be evenly divided
    add_x_low = np.arange(bbx_net[0][0], x_net[0], (x_net[0] - bbx_net[0][0]) / num_x_low)  # s2b
    add_x_high = np.arange(bbx_net[1][0], x_net[-1], (x_net[-1]- bbx_net[1][0]) / num_x_high)   #b2s
    add_z_low = np.arange(bbx_net[0][2], z_net[0], (z_net[0] - bbx_net[0][2]) / num_y_low)  # s2b
    add_z_high = np.arange(bbx_net[1][2], z_net[-1], (z_net[-1] - bbx_net[1][2]) / num_y_high)  #b2s

    # delete additional axis if too close to the temporary net
    if x_net[0] - add_x_low[-1] <= 1e-5:
        add_x_low = add_x_low[:-1]
    if add_x_high[-1] - x_net[-1] <= 1e-5:
        add_x_high = add_x_high[:-1]
    if z_net[0] - add_z_low[-1] <= 1e-5:
        add_z_low = add_z_low[:-1]
    if add_z_high[-1] - z_net[-1] <= 1e-5:
        add_z_high = add_z_high[:-1]

    x_net = np.insert(x_net, 0, add_x_low)
    x_net = np.insert(x_net, len(x_net), add_x_high[::-1])
    z_net = np.insert(z_net, 0, add_z_low)
    z_net = np.insert(z_net, len(z_net), add_z_high[::-1])




    # find vtx idx for the whole network, skip if dist(nearest point in V) is above threshold (may have gap in between)
    pos_init = np.array(np.meshgrid(x_net,z_net)).T.reshape(-1,2)       # order: (xmin,ymin) -> (xmin,ymax) -> (xmax,ymin) -> (xmax,ymax)
    pos_init = np.insert(pos_init, 1, np.zeros(pos_init.shape[0]), axis=1)
    
    dist_net, idx_net = smockPattern.tree.query(pos_init, k=1)

    ## TODO:filter if dist above threshold (gap in cloth fabric)
    idx_net = idx_net.tolist()
    pos = smockPattern.V[idx_net]
    idx_graph2mesh = idx_net
    #idx_graph2mesh = sorted(idx_net, key=idx_net.index)
    idx_mesh2graph = {mesh_idx:graph_idx for graph_idx, mesh_idx in enumerate(idx_graph2mesh)}

    stitchStartIdx = []
    stitchEndIdx = []
    for idx_stitch in range(len(stitchStartIdx_mesh)):
        stitchStartIdx.append(idx_mesh2graph[stitchStartIdx_mesh[idx_stitch]])
        stitchEndIdx.append(idx_mesh2graph[stitchEndIdx_mesh[idx_stitch]])

    # tile the network (- | \ / four kinds of edges)
    # follow the order of the pos (xmin,ymin) -> (xmin,ymax) -> (xmax,ymin) -> (xmax,ymax)
    edgeStartIdx = []
    edgeEndIdx = []
    nx = len(x_net)
    nz = len(z_net)
    for idx_x in range(nx):

        for idx_z in range(nz):

            idx_p = idx_x * nz + idx_z

            # add |^ (skip last nz)
            if idx_z != nz-1:
                edgeStartIdx.append(idx_p)
                edgeEndIdx.append(idx_p+1)

            if idx_x == nx-1:
                continue
            # following edges are skipped for the last nx (xmax, ymin) -> (xmax, ymax)

            # add -> (skip last nx)
            edgeStartIdx.append(idx_p)
            edgeEndIdx.append(idx_p+nz)

            # add /> (skip last nz & nx)
            if idx_z != nz-1:
                edgeStartIdx.append(idx_p)
                edgeEndIdx.append(idx_p+nz+1)

            # add \> (skip first nz & last nx)
            if idx_z != 0:
                edgeStartIdx.append(idx_p)
                edgeEndIdx.append(idx_p+nz-1)

    print("number of total 2D nodes: %d" %(len(pos)))

    return pos, idx_graph2mesh, idx_mesh2graph, stitchStartIdx, stitchEndIdx, edgeStartIdx, edgeEndIdx, len_all_stitch, nx, nz


class MassSpringNetwork:
    """
    solve stitching line graph as mass spring system

    Parameters:

    k - stiffness
    restL: rest length of spring
    coeff_damp: velocity damping coefficient
    dt: time step
    t: total time
    num_iter: total number of iteration

    idx_graph2mesh (list): idx_graph2mesh[graph_idx] = vtx_idx
    idx_mesh2graph (dict): idx_mesh2graph[vtx_idx] = graph_idx


    """
    def __init__(self, smockPattern:EngSmockPattern, shrink_percentage=0.5, spring_stiffness=1.0, edge_stiffness=10.0, spring_rest_length=0.1,\
                 coeff_damp=0.9, dt=0.1, num_iter=1000):

        self.k = spring_stiffness
        self.k_e = edge_stiffness
        self.restL = spring_rest_length
        self.coeff_damp = coeff_damp
        self.dt = dt
        self.t = 0.0
        self.num_iter = num_iter


        # generate regular network & connect them to the initial mesh
        self.pos, self.idx_graph2mesh, self.idx_mesh2graph, self.stitchStartIdx, \
            self.stitchEndIdx, self.edgeStartIdx, self.edgeEndIdx, self.len_all_stitch, self.net_x, self.net_z = gen_network(smockPattern)     #self.len_all_stitch, 

        self.init_pos = copy.copy(self.pos)

        self.stop_len_stitch = shrink_percentage * self.len_all_stitch
        print("len_all_stitch: %.4f" %(self.len_all_stitch))


        self.vel = np.zeros_like(self.pos)

        self.restL_adaptive = None

        ## differentiate front and back stitching lines for visualization
        ## can be commented if not using animation
        self.frontStitchStartIdx = []
        self.frontStitchEndIdx = []
        self.backStitchStartIdx = []
        self.backStitchEndIdx = []
        allLineIdx = [idx for idx in range(len(smockPattern.StitchLines))]
        backLineIdx = list(set(allLineIdx) - set(smockPattern.frontLineIdx))

        for idx_f in smockPattern.frontLineIdx:
            self.frontStitchStartIdx.append(self.stitchStartIdx[idx_f])
            self.frontStitchEndIdx.append(self.stitchEndIdx[idx_f])
        for idx_b in backLineIdx:
            self.backStitchStartIdx.append(self.stitchStartIdx[idx_b])
            self.backStitchEndIdx.append(self.stitchEndIdx[idx_b])



        

    
    def cal_acc(self, springStartIdx=None, springEndIdx=None):
        '''
        calculate the acceleration of the spring system at each step
        TODO: modify in Gauss-Seidel fashion
        '''
        if springStartIdx is None:
            springStartIdx = self.stitchStartIdx
        if springEndIdx is None:
            springEndIdx = self.stitchEndIdx
        diff_L_vec = self.pos[springEndIdx] - self.pos[springStartIdx]
        diff_L = np.linalg.norm(diff_L_vec, axis=1, keepdims=True)

        # dL = diff_L - self.restL
        if self.restL_adaptive is None:
            dL = diff_L - self.restL
        else:
            dL = diff_L - self.restL_adaptive

        acc_end = - self.k * dL * np.divide(diff_L_vec, diff_L, out=np.zeros_like(diff_L_vec), where=diff_L!=0)     # acc for point at end idx

        acc = np.zeros_like(self.pos)   # acceleration for each point aligned with self.pos

        for sprI in range(acc_end.shape[0]):
            acc[springStartIdx[sprI], :] -= acc_end[sprI]
            acc[springEndIdx[sprI], :] += acc_end[sprI]

        return acc


    def pre_calculate_net(self):
        '''
        store the initial length of the network as the rest length (for edges only)
        '''
        rest_edge_vec = self.pos[self.edgeStartIdx] - self.pos[self.edgeEndIdx]
        rest_edge = np.linalg.norm(rest_edge_vec, axis=1, keepdims=True)

        rest_stitch_vec = self.pos[self.stitchStartIdx] - self.pos[self.stitchEndIdx]
        rest_stitch = np.linalg.norm(rest_stitch_vec, axis=1, keepdims=False)


        return rest_edge, rest_stitch

    def pre_calculate_stitch(self, rest_edge):
        '''
        refine 2: adaptive rest length
        store the rest length for stitching lines
        '''
        epsilon = 0.01   # thickness (extreme rest length) 0.01 for smocking pattern, 0.001 for simple stitch
        edge_vec = self.pos[self.stitchStartIdx] - self.pos[self.stitchEndIdx]
        edge = np.linalg.norm(edge_vec, axis=1, keepdims=False)
        cos_edge = (np.abs(edge_vec[:,2])).reshape(-1, 1)       # * self.restL   * rest_edge
        
        # if cos_edge < epsilon: replace with epsilon (lower bound)
        self.restL_adaptive = np.zeros_like(cos_edge)
        np.clip(cos_edge, epsilon, None, out=self.restL_adaptive)

        return np.sum(edge)
        


    def solve_inextensibility(self, rest_edges, acc):
        '''
        for network edges that are inextensible
        option 1. add additional acc for stretch (selected now) possibly unstable
        option 2. directly modify position for stretch

        also set a lower bound to (kind of) avoid self intersection?
        '''

        epsilon = 0.01   # thickness (extreme rest length)   # TODO: really needed? pleat may overlap after projection to 2D
        #epsilon = 0.001      # for simple stitch
        diff_L_vec = self.pos[self.edgeEndIdx] - self.pos[self.edgeStartIdx]
        diff_L = np.linalg.norm(diff_L_vec, axis=1, keepdims=True)
        dL = diff_L - rest_edges
        
        stretch_mask = (dL > 0.0)  # skip compression force
        # acc_mask = stretch_mask.astype(float)

        ## refine 1: add extra compression force
        ## TODO: need to be eliminated for stitching lines??
        compress_mask = (diff_L < epsilon)    # enforce force for too close nodes
        acc_mask = np.logical_or(stretch_mask, compress_mask).astype(float)
        # refine 1.1: for stability: if compress, modify the dL to (l - epsilon) to decrease the possible compression force
        dL = (diff_L - epsilon) * compress_mask.astype(float) + dL * stretch_mask.astype(float)

        acc_end = - acc_mask * self.k_e * dL * np.divide(diff_L_vec, diff_L, out=np.zeros_like(diff_L_vec), where=diff_L!=0)     # acc for point at end idx

        for sprI in range(acc_end.shape[0]):
            acc[self.edgeStartIdx[sprI], :] -= acc_end[sprI]
            acc[self.edgeEndIdx[sprI], :] += acc_end[sprI]

        return acc


    def sim_step(self, acc):
        '''
        step calculation used for possible per-step-deformation
        '''

        self.pos += (self.vel + (acc * self.dt / 2.0)) * self.dt
        self.vel += acc * self.dt 
        self.vel *= self.coeff_damp
        self.t += self.dt
    

    def run(self, rest_edges, rest_stitches, withAnime=False):
        '''
        simulate the ultimate static result of the system (may require a stop threshold as hyper param?)
        '''

        if withAnime:
            fig = plt.figure()
            # y-axis equals to zero (omit)
            # axis fixed
            x_min = np.min(self.pos[:,0])-0.1
            x_max = np.max(self.pos[:,0])+0.1
            ## for correct mesh align with input pattern
            # y_min = np.min(self.pos[:,2])-1.0
            # y_max = np.max(self.pos[:,2])+1.0

            ## for blender save with inverted Z-axis
            y_min = -np.max(self.pos[:,2])-0.1
            y_max = -np.min(self.pos[:,2])+0.1

        #acc = np.zeros_like(self.pos)
        total_iter = self.num_iter
        len_stitch_prev = self.len_all_stitch
        img_cnt = 0

        for n_iter in range(self.num_iter):

            ## add adaptive stitch rest length iteratively
            len_stitch_now = self.pre_calculate_stitch(rest_stitches)

            # edge_vec = self.pos[self.stitchStartIdx] - self.pos[self.stitchEndIdx]
            # edge = np.linalg.norm(edge_vec, axis=1, keepdims=False)
            # len_stitch_now = np.sum(edge)

            if len_stitch_now - self.stop_len_stitch < thred_len_stop:      # threshold of stop criterion
                total_iter = n_iter+1
                print("======early stopping at "+str(total_iter)+" iteration======")
                break

            acc = self.cal_acc()

            # add network structure constraints
            acc = self.solve_inextensibility(rest_edges, acc)

            self.sim_step(acc)

            # if withAnime and (n_iter == self.num_iter-1 or n_iter == 0):
            # if withAnime and ((n_iter+1) % 20 == 0 or n_iter == 0):
            # if withAnime and (len_stitch_prev - len_stitch_now >= 0.5 or n_iter == 0):
            if withAnime:
                ## visualization 
                plt.cla()
                plt.plot(self.pos[[self.edgeStartIdx, self.edgeEndIdx], 0], -self.pos[[self.edgeStartIdx, self.edgeEndIdx], 2], color="silver", alpha=0.5, linewidth=0.8, solid_capstyle='butt')
                plt.plot(self.pos[[self.frontStitchStartIdx, self.frontStitchEndIdx], 0], -self.pos[[self.frontStitchStartIdx, self.frontStitchEndIdx], 2], color="#38BBA1", linewidth=2.5, solid_capstyle='butt')
                plt.plot(self.pos[[self.backStitchStartIdx, self.backStitchEndIdx], 0], -self.pos[[self.backStitchStartIdx, self.backStitchEndIdx], 2], color="#AB8AE0", linewidth=2.5, solid_capstyle='butt')

                # plt.scatter(self.pos[:,0], self.pos[:,2])  # vertices 
                plt.xlim((x_min, x_max))
                plt.ylim((y_min, y_max))
                fig.set_size_inches(8, 7)
                # fig.suptitle('iter: '+str(n_iter)+'\nlen(%): '+str(len_stitch_now*100.0/self.len_all_stitch))
                plt.axis('off')
                plt.rcParams["figure.facecolor"] = "#f7f2ee"
                
                # plt.savefig('..\data\springGraph\sprNet2D'+str(img_cnt).zfill(5)+'.png',dpi=500, bbox_inches='tight', facecolor="#f7f2ee")
                # img_cnt += 1

                plt.pause(0.001)

                # len_stitch_prev = len_stitch_now

        if withAnime:

            plt.cla()
            plt.plot(self.pos[[self.edgeStartIdx, self.edgeEndIdx], 0], -self.pos[[self.edgeStartIdx, self.edgeEndIdx], 2], color="silver", alpha=0.5, linewidth=0.8, solid_capstyle='butt')
            plt.plot(self.pos[[self.frontStitchStartIdx, self.frontStitchEndIdx], 0], -self.pos[[self.frontStitchStartIdx, self.frontStitchEndIdx], 2], color="#38BBA1", linewidth=2.5, solid_capstyle='butt')
            plt.plot(self.pos[[self.backStitchStartIdx, self.backStitchEndIdx], 0], -self.pos[[self.backStitchStartIdx, self.backStitchEndIdx], 2], color="#AB8AE0", linewidth=2.5, solid_capstyle='butt')
            # plt.scatter(self.pos[:,0], self.pos[:,2])   
            plt.xlim((x_min, x_max))
            plt.ylim((y_min, y_max))
            fig.set_size_inches(8, 7)
            plt.axis('off')
            # fig.suptitle('iter: '+str(n_iter)+'\nlen(%): '+str(len_stitch_now*100.0/self.len_all_stitch))

            # plt.savefig('..\data\springGraph\sprNet2D'+str(img_cnt).zfill(5)+'.png',dpi=500, bbox_inches='tight', facecolor="#f7f2ee")
            plt.show()
        

        # use only stitched node positions as output constraints
        stitch_idx = list(set(self.stitchStartIdx + self.stitchEndIdx)) 
        stitch_pos = self.pos[stitch_idx]
        stitch_idx_graph2mesh = []
        for idx in stitch_idx:
            stitch_idx_graph2mesh.append(self.idx_graph2mesh[idx])
        stitch_idx_mesh2graph = {mesh_idx:graph_idx for graph_idx, mesh_idx in enumerate(stitch_idx_graph2mesh)}
        print("------spring system simulation done!------")
        print(len_stitch_now)
        return stitch_pos, stitch_idx_graph2mesh, stitch_idx_mesh2graph, total_iter


    def get_net_pos(self):
        return self.pos, self.idx_graph2mesh
    
    def get_net_size(self):
        return self.net_x, self.net_z
    