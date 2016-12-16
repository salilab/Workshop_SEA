#!/usr/bin/env python
#####################################################
# Last Update: December 4, 2016
# by Ilan E. Chemmama
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################

import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf
import os, sys

import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.proteomics

def create_rotational_symmetry2(rps, maincopy, copies, rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0),
                                nSymmetry=None, skip_gaussian_in_clones=False):
    '''
    The copies must not contain rigid bodies.
    The symmetry restraints are applied at each leaf.
    '''
    from math import pi
    rps.representation_is_modified = True
    ncopies = len(copies) + 1
    main_hiers = IMP.atom.get_leaves(rps.hier_dict[maincopy])

    for k in range(len(copies)):
        if (nSymmetry is None):
            rotation_angle = 2.0 * pi / float(ncopies) * float(k + 1)
        else:
            if ( k % 2 == 0 ):
                rotation_angle = 2.0 * pi / float(nSymmetry) * float((k + 2) / 2)
            else:
                rotation_angle = -2.0 * pi / float(nSymmetry) * float((k + 1) / 2)
        rotation3D = IMP.algebra.get_rotation_about_axis(rotational_axis, rotation_angle)
        sm = IMP.core.TransformationSymmetry(rotation3D)
        clone_hiers = IMP.atom.get_leaves(rps.hier_dict[copies[k]])

        lc = IMP.container.ListSingletonContainer(rps.m)
        for n, p in enumerate(main_hiers):
            if (skip_gaussian_in_clones):
                if (IMP.core.Gaussian.get_is_setup(p)) \
                            and not (IMP.atom.Fragment.get_is_setup(p) or \
                                             IMP.atom.Residue.get_is_setup(p)):
                    continue
            pc = clone_hiers[n]
            IMP.core.Reference.setup_particle(pc.get_particle(), p.get_particle())
            lc.add(pc.get_particle().get_index())

            c = IMP.container.SingletonsConstraint(sm, None, lc)
            rps.m.add_score_state(c)
            print("Completed setting " + str(maincopy) + \
                          " as a reference for " + str(copies[k]) + \
                          " by rotating it in " + str(rotation_angle / 2.0 / pi * 360) + \
                          " degree around the " + str(rotational_axis) + " axis.")
        rps.m.update()

##Input files and parameters
datadirectory = "../data/"
topology_file = datadirectory+"topology_all.txt"

##Trajectory options
num_frames = 30
num_mc_steps = 10

rb_max_trans = 2.00
rb_max_rot = 0.04
bead_max_trans = 3.00


#####################################################
# Create hierarchies and rigid bodies and flexible parts
# for bead representations
#####################################################
m = IMP.Model()

topology = IMP.pmi.topology.TopologyReader(topology_file)
domains = topology.component_list

rigid_bodies=[["Npr2_2"],
              ["Npr2_4"],
              ["Npr2_6"],
              ["Npr3_1"],
              ["Npr3_3", "Npr3_5"],
              ["Npr3_7"],
              ["Npr3_9"],
              ["Npr3_11"],
              ["Sec13_2", "Sec13_4"],
              ["SEA1_2"],
              ["SEA1_4", "SEA1_6", "SEA1_8"],
              ["SEA1_10"],
              ["SEA2_2", "SEA2_4", "SEA2_6", "SEA2_8"],
              ["SEA2_10"],
              ["SEA3_2", "SEA3_4", "SEA3_6", "SEA3_8"],
              ["SEA3_10"],
              ["SEA3_12"],
              ["Seh1.1_2", "Seh1.1_4"],
              ["SEA4.1_2", "SEA4.1_4", "SEA4.1_6", "SEA4.1_8", "SEA4.1_10"],
              ["SEA4.1_12", "SEA4.1_14"],
              ["SEA4.1_16", "SEA4.1_18"]
              ]

bm = IMP.pmi.macros.BuildModel(m,
                               component_topologies=domains,
                               list_of_rigid_bodies=rigid_bodies)
representation = bm.get_representation()

for nc,component in enumerate(domains):
    name = component.name
    sel = IMP.atom.Selection(representation.prot,molecule=name)
    ps = sel.get_selected_particles()
    clr = IMP.display.get_rgb_color(float(nc)/len(domains))
    for p in ps:
        if not IMP.display.Colored.get_is_setup(p):
            IMP.display.Colored.setup_particle(p,clr)
        else:
            IMP.display.Colored(p).set_color(clr)

#Symmetry Restraints - z-axis
create_rotational_symmetry2(representation, "SEA4.1", ["SEA4.2","SEA4.3"])
create_rotational_symmetry2(representation, "Seh1.1", ["Seh1.2","Seh1.3"])
for rt in ['SEA4.2', 'SEA4.3', 'Seh1.2', 'Seh1.3']:
    hs = IMP.pmi.tools.select_by_tuple(representation, rt)
    representation.remove_floppy_bodies(hs)

#--------------------------
# Define Degrees of Freedom
#--------------------------
# Add default mover parameters to simulation
representation.set_rigid_bodies_max_rot(rb_max_rot)
representation.set_floppy_bodies_max_trans(bead_max_trans)
representation.set_rigid_bodies_max_trans(rb_max_trans)

outputobjects = [] # reporter objects (for stat files)
sampleobjects = [] # sampling objects

# Add the movers to the sample and output object lists
outputobjects.append(representation)
sampleobjects.append(representation)

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------
# Excluded Volume Restraint
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(representation,
                                                             resolution=10)
ev.add_to_model()
outputobjects.append(ev)

# Crosslinks
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
kw.set_id_score_key(None)
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb.create_set_from_file(datadirectory+'xlinks_SEA.csv')

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(representation=representation,
                                                                            CrossLinkDataBase=xldb,
                                                                            length=21,
                                                                            label="SeaComplex",
                                                                            resolution=1.0,
                                                                            slope=0.02)

xl1.rs.set_weight(25.0)
xl1.add_to_model()

sampleobjects.append(xl1)
outputobjects.append(xl1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)

#Composite Restraints
weight = 1.0
crd={}

crd["Npr2_dNpr2_497_615_P"]=[(1,496,"Npr2"),"Npr3"]
crd["SEA4_dSEA4:931-1038_P"]=[(1,930,"SEA4.1"),(1,930,"SEA4.2"),(1,930,"SEA4.3"),"Seh1.1","Seh1.2","Seh1.3"]

for key in crd:
    cr=IMP.pmi.restraints.proteomics.ConnectivityRestraint(representation,
                                                           crd[key],
                                                           resolution=100.0)

    cr.add_to_model()
    cr.set_label(key)
    outputobjects.append(cr)

#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol

# Randomize the initial configuration before sampling
representation.shuffle_configuration(50)
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    representation,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1,],    # allows XLs to be drawn in the RMF files
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=1,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="output")

# Start Sampling
mc1.execute_macro()
