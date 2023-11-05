import numpy as np
from arap import ARAP
from scipy import spatial
import itertools
from collections import Counter
import copy

import igl
import scipy as sp
import os

height = 1.0        # height value for mid stitch nodes (around 1.0 for regular pattern, 0.05 for simple stitches)
height_percent = 0.5    # height percentage according to the initial stitching line (<0.5)

class EngSmockPattern:
    '''
    smocking pattern for english smocking
    '''
    
    def __init__(self, V, F, StitchLines, frontLineIdx=None, arap_weight=1.0, num_iter=1):
        self.V = copy.copy(V)
        self.org_V = copy.copy(V)

        # arap method use triangular mesh
        self.F = np.array(F) 

        self.StitchLines = StitchLines
        self.arap_weight = arap_weight
        self.num_iter = num_iter

        if frontLineIdx is None:
            self.frontLineIdx = [i for i in range(0, len(StitchLines), 2)]  # default: even idx of lines are front stitch lines)
        else:
            self.frontLineIdx = frontLineIdx    # line ids for front sewing

        self.LineEnds = []

        # prepare initial length of stitching lines for height constraint
        self.stitchLenInit = np.array([height]*len(StitchLines))

        # prepare KDTree to search mid point & possible net point
        self.tree = spatial.KDTree(self.org_V)

        # search for all the mid point vertices (pre-step for pleat constraints)
        # here the mid_points position belongs to the original mesh (not introduce simulation position)
        mid_points = []
        for stitch in self.StitchLines:
            mid_points.append(np.mean(self.org_V[stitch], axis=0, keepdims=False))
        mid_points = np.array(mid_points)

        _, self.idx_p_mid = self.tree.query(mid_points, k=1)     # mid point's vert idx in mesh (follow stitching line's rank)


    @property
    def n_pts(self):
        return self.V.shape[0]

    @property
    def n_lines(self):
        return len(self.StitchLines)
    
    @property
    def n_lines_front(self):
        return len(self.frontLineIdx)

   
    def extract_constraints(self, constraint_pos:np.ndarray, idx_graph2mesh:list, idx_mesh2graph:dict, idx_p_mid=None):
        '''
        generate constraints for underlay and pleat nodes
        Underlay constraint: position same as simulator output
        Pleat constraint: select (nearest) mid point of stitching line && shift in y-direction (both up & down)
        '''
        if idx_p_mid is None:
            idx_p_mid = self.idx_p_mid

        constraints = {}

        # underlay constraint
        for i in range(len(idx_graph2mesh)):
            constraints[idx_graph2mesh[i]] = constraint_pos[i]

        ## pleat constraint
        allLineIdx = [idx for idx in range(len(self.StitchLines))]
        backLineIdx = list(set(allLineIdx) - set(self.frontLineIdx))

        
        # for front line, shift y downward (<0)
        for idx_l_f in self.frontLineIdx:
            idx_mesh_f = idx_p_mid[idx_l_f]     # mid point's vert idx for idx_l_f-th stitching line
            # position needs calculation based on simulated position
            pos_f = np.zeros(3)
            stitch_f = self.StitchLines[idx_l_f]
            for stitch_f_node in stitch_f:    # stitch node in mesh idx
                pos_f += constraint_pos[idx_mesh2graph[stitch_f_node]]
            pos_f /= float(len(stitch_f))   # avg the position

            # assume two nodes for each stitchig line
            stitch_vec_current = constraint_pos[idx_mesh2graph[stitch_f[0]]] - constraint_pos[idx_mesh2graph[stitch_f[1]]]
            stitch_len_current = np.linalg.norm(stitch_vec_current)

            pos_f[1] = -(self.stitchLenInit[idx_l_f] - stitch_len_current)* height_percent  # front line means deformation is backward
            constraints[idx_mesh_f] = pos_f

        # for back line, shift y upward (>0)
        for idx_l_b in backLineIdx:
            idx_mesh_b = idx_p_mid[idx_l_b]     # mid point's vert idx for idx_l_b-th stitching line

            # position needs calculation based on simulated position
            pos_b = np.zeros(3)
            
            stitch_b = self.StitchLines[idx_l_b]
            for stitch_b_node in stitch_b:    # stitch node in mesh idx
                pos_b += constraint_pos[idx_mesh2graph[stitch_b_node]]
            pos_b /= float(len(stitch_b))   # avg the position

            # assume two nodes for each stitchig line
            stitch_vec_current = constraint_pos[idx_mesh2graph[stitch_b[0]]] - constraint_pos[idx_mesh2graph[stitch_b[1]]]
            stitch_len_current = np.linalg.norm(stitch_vec_current)

            pos_b[1] = (self.stitchLenInit[idx_l_b] - stitch_len_current)* height_percent  # back line means deformation is upward
            constraints[idx_mesh_b] = pos_b
        
        return constraints



    def deform(self, constraint_pos:np.ndarray, idx_graph2mesh:list, idx_mesh2graph:dict):
        '''
        1. deform the mesh using the constraint_pos (position of simulated spring ends) and according verts idx
        2. apply ARAP
        '''

        self.constraints = self.extract_constraints(constraint_pos, idx_graph2mesh, idx_mesh2graph)

        arap = ARAP(self.V.T, self.F, self.constraints.keys(), self.arap_weight)
        print("------ARAP initialization done!------")
        return arap(self.constraints, self.num_iter, self.V).T


    def generate_init_config(self):
        '''
        add y-direction shift for mid points of each stitching line
        '''
        
        allLineIdx = [idx for idx in range(len(self.StitchLines))]
        backLineIdx = list(set(allLineIdx) - set(self.frontLineIdx))

        # prepare KDTree to search mid point
        tree = spatial.KDTree(self.V)

        # search for all the mid point vertices (pre-step for pleat constraints)
        # here the mid_points position belongs to the original mesh (not introduce simulation position)
        mid_points = []
        for stitch in self.StitchLines:
            mid_points.append(np.mean(self.V[stitch], axis=0, keepdims=False))
        mid_points = np.array(mid_points)

        _, idx_p_mid = tree.query(mid_points, k=1)     # mid point's vert idx in mesh (follow stitching line's rank)

        V_init = np.array(self.V)    # deep copy of V
        for idx_l_f in self.frontLineIdx:
            idx_mesh_f = idx_p_mid[idx_l_f]     # mid point's vert idx for idx_l_f-th stitching line
            V_init[idx_mesh_f, 1] = -self.stitchLenInit[idx_l_f]* height_percent
        
        for idx_l_b in backLineIdx:
            idx_mesh_b = idx_p_mid[idx_l_b]     # mid point's vert idx for idx_l_f-th stitching line
            V_init[idx_mesh_b, 1] = self.stitchLenInit[idx_l_b]* height_percent
            
        return V_init, self.F, self.StitchLines


