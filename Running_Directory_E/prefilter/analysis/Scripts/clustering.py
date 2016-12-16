import IMP
import IMP.pmi
import macros_e29
import sys

#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='generate clusters of the RMF solutions')
parser.add_argument('-mpi', action="store", dest="is_mpi", help="mpi enabled")
parser.add_argument('-preload', action="store", dest="load_distance_matrix_file", help="skip the matrix calcuklation and read the precalculated matrix")
parser.add_argument('-nmods', action="store", dest="nbestscoringmodels", help="number of models to be clustered")
parser.add_argument('-nclusters', action="store", dest="nclusters", help="number of clusters to be used by kmeans algorithm")
parser.add_argument('-prefilter', action="store", dest="prefiltervalue", help="prefilter the models by score")
inputs = parser.parse_args()

# mpi enabled
if (inputs.is_mpi=="True") or (inputs.is_mpi=="true") or (inputs.is_mpi=="Yes") or (inputs.is_mpi=="yes") :
    inputs.is_mpi = True
else:
    inputs.is_mpi = False

# skip the matrix calcuklation and read the precalculated matrix
if (inputs.load_distance_matrix_file=="True") or (inputs.load_distance_matrix_file=="true") or (inputs.load_distance_matrix_file=="Yes") or (inputs.load_distance_matrix_file=="yes") :
    inputs.load_distance_matrix_file = True
else:
    inputs.load_distance_matrix_file = False

# number of models to be clustered
if inputs.nbestscoringmodels==None:
    inputs.nbestscoringmodels = 500

# number of clusters to be used by kmeans algorithm
if inputs.nclusters==None:
    inputs.nclusters = 1

# prefilter the models by score
if inputs.prefiltervalue==None:
    inputs.prefiltervalue = 565.0

print inputs

is_mpi = inputs.is_mpi                                          # mpi enabled
load_distance_matrix_file = inputs.load_distance_matrix_file    # skip the matrix calcuklation and read the precalculated matrix
nbestscoringmodels = int(inputs.nbestscoringmodels)             # number of models to be clustered
nclusters = int(inputs.nclusters)                               # number of clusters to be used by kmeans algorithm
prefiltervalue = float(inputs.prefiltervalue)                   # prefilter the models by score


#####################################################
# initialize the macro
#####################################################
model=IMP.Model()

mergedirectories= []
mergedirectories.append("./all_models.999")


mc=macros_e29.AnalysisReplicaExchange0(model,
                                        stat_file_name_suffix="stat",
                                        merge_directories=mergedirectories,
                                        global_output_directory="./")

feature_list=["Total_Score",
              "ConnectivityRestraint_None",
              "CrossLinkingMassSpectrometryRestraint_Data_Score",
              "ExcludedVolumeSphere_all",
              "ExcludedVolumeSphere_ecm29"
              ]

reduced_density_dict={"SEA1":["SEA1"],
                      "SEA2":["SEA2"],
                      "SEA3":["SEA3"],
                      "SEA4.1":["SEA4.1"],
                      "SEA4.2":["SEA4.2"],
                      "SEA4.3":["SEA4.3"],

                      "Npr2":["Npr2"],
                      "Npr3":["Npr3"],
                      "Sec13":["Sec13"],
                      "Seh1.1":["Seh1.1"],
                      "Seh1.2":["Seh1.2"],
                      "Seh1.3":["Seh1.3"],

                      "SEA4.A":["SEA4.1", "SEA4.2", "SEA4.3"],
                      "Seh1.A":["Seh1.1", "Seh1.2", "Seh1.3"],

                      "Whole":["SEA1", "SEA2", "SEA3",
                               "Npr2", "Npr3", "Sec13",
                               "SEA4.1", "SEA4.2", "SEA4.3",
                               "Seh1.1", "Seh1.2", "Seh1.3"]
                      }
# list of component names needed to calculate the RMSD for the clustering
components_names_a={"SEA4.1_A":"SEA4.1",
                    "SEA4.2_A":"SEA4.2",
                    "SEA4.3_A":"SEA4.3",
                    "Seh1.1_A":"Seh1.1",
                    "Seh1.2_A":"Seh1.2",
                    "Seh1.3_A":"Seh1.3",
                    }
components_names_r={"SEA1":"SEA1",
                    "SEA2":"SEA2",
                    "SEA3":"SEA3",
                    "Sec13":"Sec13",
                    "Npr3":"Npr3",
                    "SEA4.1_A":"SEA4.1",
                    "SEA4.2_A":"SEA4.2",
                    "SEA4.3_A":"SEA4.3",
                    "Seh1.1_A":"Seh1.1",
                    "Seh1.2_A":"Seh1.2",
                    "Seh1.3_A":"Seh1.3",
                    }

mc.clustering_rmfs_no_stat("SimplifiedModel_Total_Score_None",
                           "rmf_file",
                           "rmf_frame_index",
                           rmfsdir="./all_models.999",
                           alignment_components=components_names_r,
                           number_of_best_scoring_models=nbestscoringmodels,
                           rmsd_calculation_components=components_names_r,
                           outputdir="kmeans_"+str(nbestscoringmodels)+"_"+str(nclusters)+"/",
                           feature_keys=feature_list,
                           load_distance_matrix_file=load_distance_matrix_file,
                           skip_clustering=False,
                           display_plot=True,
                           exit_after_display=False,
                           get_every=1,
                           number_of_clusters=nclusters,
                           voxel_size=5.0,
                           density_custom_ranges=reduced_density_dict
                           )
