

import os

figure_directory = 'figures/chignolin_test'
if not os.path.exists(figure_directory):
    os.makedirs(figure_directory)
figure_filenames = {
                    "cluster_rate"       : "cluster_exploration.png",
                    "cluster_enhancement": "cluster_enhancement.png",
                   }
figure_paths = {
                k:os.path.join(figure_directory,v)
                for k,v in figure_filenames.items()
               }
# files navigation
topdir = './'
ls_extensions    = {
                    "timestep"  : 0.004,   # 4 picoseconds
                    "topfile"   : 'chignolin.protein.pdb',
                    "directory" : 'data/',
                    "filename"  : '../../jrossyra/admdrp/admdrp/projects/chignolin_longstart_AS/trajs/*/protein.dcd',
                    "n_trajs"   : 25,
                    "first"     : 1,
                    #"last"      : ,
                    "trajfiles" : list(),
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }
ex_extensions = {
                    "timestep"  : 0.004,   # 4 picoseconds
                    "topfile"   : 'chignolin.protein.pdb',
                    "directory" : 'data/',
                    "filename"  : '../../jrossyra/admdrp/admdrp/projects/chignolin_x/trajs/*/protein.dcd',
                    "n_trajs"   : 50,
                    "first"     : 501,
                    #"last"      : ,
                    "trajfiles" : list(),
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }
ex_parameters = {
                    "timestep"  : 0.004,   # 4 picoseconds
                    "topfile"   : 'chignolin.protein.pdb',
                    "directory" : 'data/',
                    "filename"  : '../../jrossyra/admdrp/admdrp/projects/chignolin_x/trajs/*/protein.dcd',
                    "n_trajs"   : 50,
                    "first"     : 0,
                    "last"      : 500,
                    "trajfiles" : list(),
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }
ls_parameters = {
                    "timestep"  : 0.004, # 4 picoseconds
                    "topfile"   : 'chignolin.protein.pdb',
                    "directory" : 'data/',
                    "filename"  : '../../jrossyra/admdrp/admdrp/projects/chignolin_longstart/trajs/*/protein.dcd',
                    "n_trajs"   : 1,
                    "first"     : 0,
                    "last"      : 1,
                    "trajfiles" : list(),
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }
st_parameters = {
                    "timestep"  : 0.200, # 4 picoseconds
                    "topfile"   : 'CLN025-0-protein_fixed_noH.pdb',
                    "directory" : 'data/',
                    "filename"  : 'CLN025-0-protein_all.dcd',
                    "n_trajs"   : 1,
                    "trajfiles" : list(),
                    "stride"    : None,
                    "epochsize" : None,
                    "totaldata" : None,
                    "n_clusters": list(),
                    "datashape" : list(),
                   }
analysis = {
              "topfile"         : "data/chignolin.ca.pdb",
              #"master_source"   : "Single Trajectory",
              "master_source"   : "Single Trajectory",
              "master_analysis" : {
                  "atomselection"   : "name CA",
                  "topfile"         : "data/chignolin.ca.pdb",
                  "dim_reduction"   : { "tica" : {
                                                   "lag"    : 75, # 15 ns in Anton Frames
                                                   "stride" : 5,  # 1  ns in Anton Frames
                                                   "dim"    : 6,
                                      }          },
                  "clustering"      : { "cluster_kmeans" : {
                                                #"k"         : None,
                                                "k"         : 1000, # Opt was ~900 w/ 22ns lag
                                                "max_iter"  : 50,
                                                "stride"    : 5,
                                      }                    },
                  "msm"      : { "estimate_markov_model" : {
                                                #"reversible": True,
                                                "lag"       : 75,
                                                           }
                                      },
                  },
              "atomselection"   : "name CA",
              "analyze"         : {
                                   #"Long Trajectory",
                                   "Long Extensions",
                                   #"Explore Macrostates",
                                   "Single Trajectory",
                                   "Explore Extensions",
                                  },
              "dim_reduction"   : { "tica" : {
                                               "lag"    : 75,
                                               "stride" : 2,
                                               "dim"    : 6,
                                             }
                                  },
              "clustering"      : { "cluster_kmeans" : {
                                            #"k"         : None,
                                            "k"         : 500,
                                            "max_iter"  : 50,
                                            "stride"    : 5,
                                                       }
                                  },
              "msm"      : { "estimate_markov_model" : {
                                            "reversible": True,
                                            "lag"       : 75,
                                                       }
                                  },
            }

'''
Name of datasets and ref to each dataset dict. Must contain analysis field
with ref to analysis dict, which describes the processing to do.
'''
datasets = {
            "Single Trajectory"  : st_parameters,
            "Long Trajectory"    : ls_parameters,
            "Long Extensions"    : ls_extensions,
            #"Explore Macrostates": ex_parameters,
            "Explore Extensions" : ex_extensions,
            #"Random Frames"      : rs_parameteres,

            "analysis"           : analysis,
           }
'''
Extension to name for feature trajs used in analysis
'''
namekey = '---stride-'
