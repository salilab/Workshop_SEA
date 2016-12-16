import sys, os, os.path
import math
from math import log, sqrt
import glob
import numpy as np
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros

import IMP.rmf
import RMF

def compute_precision(path, component):
    conform = []
    name = []
    num = 0
    dof = 0
    for file in glob.glob("%s/*.rmf3" % path):
        name.append(file.split("/")[-1])
        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        partcoord = []

        for c in h.get_children():
            for s in c.get_children():
                if str(component) in s.get_name():
                    for ps in IMP.atom.get_leaves(s):
                        p = IMP.core.XYZ(ps)
                        partcoord.append(p.get_x())
                        partcoord.append(p.get_y())
                        partcoord.append(p.get_z())
        conform.append(partcoord)
        partcoord = []
        num = num + 1
    std = np.array(conform).std(0)*np.array(conform).std(0)
    return np.sum(std)

def get_dof(file0, component):
    for file in glob.glob(file0):
        partcoord = []
        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]

        for c in h.get_children():
            for s in c.get_children():
                if str(component) in s.get_name():
                    for ps in IMP.atom.get_leaves(s):
                        p = IMP.core.XYZ(ps)
                        partcoord.append(p.get_x())
                        partcoord.append(p.get_y())
                        partcoord.append(p.get_z())
    return len(partcoord)

def get_alphak(alpha_kminus1,k,ndof):
    if k==2:
        alpha_k=1-0.75/ndof
    elif k>2:
        alpha_k=alpha_kminus1+(1-alpha_kminus1)/6.0
    return alpha_k

total_num_clusters=int(sys.argv[1]) # e.g. kmeans 1 through 8 will be 8
path0=sys.argv[2]
component=sys.argv[3]

s_k_list=[]
a_k_list =[]

ndof = get_dof(path0+"1/cluster.0/0-*.rmf3", component)
for k in range(1,total_num_clusters+1):
    s_k=0.0 # sum of distortions
    for i in range(k):
        curr_precision = compute_precision(path0+str(k)+"/cluster."+str(i), component)
        s_k += curr_precision

    if k==1:
        a_k_list.append(1-0.75/ndof)
        s_k_list.append(s_k)
        print k, 1-0.75/ndof, 1
        continue
    else:
        a_k = get_alphak(a_k_list[-1],k,ndof)
        a_k_list.append(a_k)
        print k, a_k, s_k/(a_k*s_k_list[-1])
        s_k_list.append(s_k)
        continue
