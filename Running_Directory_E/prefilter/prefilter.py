from __future__ import print_function
import IMP
import IMP.pmi
import sys,os
import glob
import argparse

parser = argparse.ArgumentParser(description='generate clusters of the RMF solutions')
parser.add_argument('--mpi', action="store", dest="is_mpi", help="mpi enabled")
parser.add_argument('--preload', action="store", dest="load_distance_matrix_file", help="skip the matrix calcuklation and read the precalculated matrix")
parser.add_argument('--nmods', action="store", dest="nbestscoringmodels", help="number of models to be clustered")
parser.add_argument('--nclusters', action="store", dest="nclusters", help="number of clusters to be used by kmeans algorithm")
parser.add_argument('--prefilter', action="store", dest="prefiltervalue", help="prefilter the models by score")
parser.add_argument('--dir', action="store", dest="is_dir", help="input directory")
inputs = parser.parse_args()

if (inputs.is_mpi=="True") or (inputs.is_mpi=="true") or (inputs.is_mpi=="Yes") or (inputs.is_mpi=="yes") :
    inputs.is_mpi = True
else:
    inputs.is_mpi = False
if (inputs.load_distance_matrix_file=="True") or (inputs.load_distance_matrix_file=="true") or (inputs.load_distance_matrix_file=="Yes") or (inputs.load_distance_matrix_file=="yes") :
    inputs.load_distance_matrix_file = True
else:
    inputs.load_distance_matrix_file = False
if inputs.nbestscoringmodels==None:
    inputs.nbestscoringmodels = 50
if inputs.nclusters==None:
    inputs.nclusters = 1
if inputs.prefiltervalue==None:
    inputs.prefiltervalue = 565.0
if inputs.is_dir==None:
    exit
print(inputs)

is_mpi = inputs.is_mpi
load_distance_matrix_file = inputs.load_distance_matrix_file
nbestscoringmodels = int(inputs.nbestscoringmodels)
nclusters = int(inputs.nclusters)
prefiltervalue = float(inputs.prefiltervalue)
is_dir = inputs.is_dir

import macros_sea
model=IMP.Model()

print(is_dir)
mergedirectories= []
mergedirectories.append(is_dir)
mc=macros_sea.AnalysisReplicaExchange0(model,
                                       stat_file_name_suffix="stat",     # don't change
                                       merge_directories=mergedirectories,
                                       global_output_directory='output')

feature_list = ["SimplifiedModel_Total_Score_None",
                "ExcludedVolumeSphere_None",
                "ConnectivityRestraint_None",
                "CrossLinkingMassSpectrometryRestraint_Data_Score_",

                "CrossLinkingMassSpectrometryRestraint_Linear_Score_SeaComplex",
                "CrossLinkingMassSpectrometryRestraint_PriorPsi_Score_SeaComplex",
                "CrossLinkingMassSpectrometryRestraint_PriorSig_Score_SeaComplex 0.0",
                "CrossLinkingMassSpectrometryRestraint_Distance_"
]

mc.clustering("SimplifiedModel_Total_Score_None",
              "rmf_file",
              "rmf_frame_index",
              rmsd_calculation_components=None,
              alignment_components=None,
              number_of_clusters=nclusters,
              feature_keys = feature_list,
              prefiltervalue=prefiltervalue,
              number_of_best_scoring_models=nbestscoringmodels,
              distance_matrix_file=None,
              outputdir="kmeans_"+str(nbestscoringmodels)+"_"+is_dir.split("/")[-1]+"/",
              skip_clustering=True,
              get_every=1
              )
