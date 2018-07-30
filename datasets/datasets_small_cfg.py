

import os

figure_directory = 'figures/testing-ntl9'

if not os.path.exists(figure_directory):
    os.makedirs(figure_directory)

figure_filenames = {
                    "cluster_rate"       : "cluster_exploration.png",
                    "cluster_enhancement": "cluster_enhancement.png",
                   }

figure_paths = {k:os.path.join(figure_directory,v)
                for k,v in figure_filenames.items()}

# files navigation
topdir = './'


#anton_parameters = {
#                    "timestep"  : 0.2,   # 200 picoseconds
#                    "topfile"   : 'NTL9-0-protein_fixed_noH.pdb',
#                    "directory" : '',#NTL9-0-protein_fixed',
#                    "filename"  : 'NTL9-3-protein_all.dcd',
#                    "trajfiles" : list(),
#                    "n_trajs"   : 1,
#                    "stride"    : None,
#                    "epochsize" : None,
#                    "totaldata" : None,
#                    "n_clusters": list(),
#                    "datashape" : list(),
#                   }

xm_parameters = {
                    "timestep"  : 0.1, # 50 picoseconds
                    "topfile"   : 'ntl9.protein.pdb',
                    "directory" : 'demo-xm-shortsegs/trajs/less',
                    "filename"  : '*/protein.dcd',
                    "trajfiles" : list(),
                    "n_trajs"   : 5,
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }

lt_parameters = {
                    "timestep"  : 0.025, # 50 picoseconds
                    "topfile"   : 'ntl9.protein.pdb',
                    "directory" : 'demo-lt-shortsegs/trajs/less',
                    "filename"  : '*/protein.dcd',
                    "trajfiles" : list(),
                    "n_trajs"   : 10,
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }

analysis_parameters = {
                        "master_topofile" : "ntl9.ca.pdb",
                        "master_source"   : "Long Trajectories",
                        "atomselection"   : "name CA",
                        "master_dataset"  : {
                                             #"Single Trajectory",
                                             "Explore Microstates",
                                             "Long Trajectories",
                                            },
                        "dim_reduction"   : { "tica" : {
                                                         "lag"    : 500,
                                                         "stride" : 20,
                                                         "dim"    : 4,
                                                       }
                                            },
                        "clustering"      : { "cluster_kmeans" : {
                                                      "k"         : None,
                                                      "max_iter"  : 50,
                                                      "stride"    : 5,
                                                                 }
                                            },
                      }

datasets = {
            #"Single Trajectory"  : anton_parameters,
            "Explore Microstates": xm_parameters,
            "Long Trajectories"  : lt_parameters,
            #"Random Frames"      : rs_parameteres,
            "analysis"           : analysis_parameters,
           }


namekey = '--stride-'

