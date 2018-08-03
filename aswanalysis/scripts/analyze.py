#!/usr/bin/env python

from pprint import pprint
import os
import copy
import numpy as np
import pyemma
from pyemma.coordinates import assign_to_centers
#import pyemma.msm as msm 
pyemma.version
pyemma.config.show_progress_bars=False
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from aswanalysis import aswtools
from datasets import datasets_chignolin_test as dcfg

#  # THESE GUYS COME ALONG FROM
#  # DATASETS_CFG import above
# TODO promotion lines on import to top registry
#   via list `_import_modules`
#   `[bind(imported) for imported in process_tools._imported_modules]`
mdtraj     = process_tools.mdtraj
subprocess = process_tools.subprocess
shlex      = process_tools.shlex
datasets   = dcfg.datasets
namekey    = dcfg.namekey
topdir     = dcfg.topdir


master_dataset   = copy.deepcopy({ nm:pars for nm,pars in datasets.items() if nm in datasets['analysis']['master_source'] })
master_analysis  = datasets['analysis']['master_analysis']
process_tools.prepare_data(master_dataset, master_analysis, topdir)
master_tica_inp, master_input_order \
                 = process_tools.prepare_tica_inputs(master_dataset, master_analysis['topfile'], True)
# nm to A conversion for inverse distances
master_tica_inp  = master_tica_inp[0] / 10
tica_func    = process_tools.build_call(master_analysis['dim_reduction'], process_tools.coor, master_tica_inp)
master_tica = tica_func()
cluster      = process_tools.build_call(master_analysis['clustering'],    process_tools.coor, master_tica_results.get_output())
master_clustering   = cluster()
master_msm_inp      = master_clustering.dtrajs
emm          = process_tools.build_call(master_analysis['msm'],           process_tools.msm,  master_msm_inp)
master_mm           = emm()
process_tools.bind_dtraj_results(master_dataset, master_input_order, master_clustering.dtrajs)
#master_tica, master_clustering, master_mm \
#                 = process_tools.analyze(master_dataset, master_analysis)


analysis_dataset      = copy.deepcopy({ nm:pars for nm,pars in datasets.items() if nm in datasets['analysis']['analyze'] })
analysis_analysis     = datasets['analysis']
process_tools.prepare_data(analysis_dataset, analysis_analysis, topdir)
analysis_tica_inp, analysis_input_order \
                      = process_tools.prepare_tica_inputs(analysis_dataset, analysis_analysis['topfile'], True)
# THIS DOESN"T WORK
#_i                    = analysis_input_order.index("Single Trajectory")
#analysis_tica_inp[_i] = analysis_tica_inp[_i] / 10
analysis_tica_trajs   = master_tica.transform(analysis_tica_inp)
#analysis_tica_trajs   = master_tica.transform(analysis_tica_inp.get_output())
analysis_dtrajs       = assign_to_centers(analysis_tica_trajs, centers=master_clustering.clustercenters, metric='minRMSD', stride=1, chunk_size=5000)
process_tools.bind_dtraj_results(analysis_dataset, analysis_input_order, analysis_dtrajs)

# pretty IMPORTANT part here
process_tools.calc_eq_explored(analysis_dataset, master_mm)



additional_dataset     = copy.deepcopy({"Single Trajectory": datasets["Single Trajectory"], "Long Extensions": datasets["ls_extensions"]})
process_tools.prepare_data(additional_dataset, analysis_analysis, topdir)
additional_tica_inp, additional_input_order \
                       = process_tools.prepare_tica_inputs(additional_dataset, analysis_analysis['topfile'])
additional_tica_trajs  = master_tica.transform(tica_inp.get_output())
additional_dtrajs      = assign_to_centers(additional_tica_trajs, centers=master_kmeans.clustercenters, metric='minRMSD', stride=1, chunk_size=5000)
process_tools.bind_dtraj_results(additional_dataset, input_order, additional_dtrajs)
process_tools.calc_eq_explored(additional_dataset)



# TODO SAVING
plot_this = analysis_dataset
#plot_this = additional_dataset

for nm,pars in plot_this.items():
    plt.plot( [i*0.2*25/1000 for i in range(len(pars['eq_explored']))], pars['eq_explored'], label=nm )
    #plt.plot( range(len(pars['eq_explored'])), pars['eq_explored'], label=nm )

plt.title( "Equilibrium Distribution Explored with Different MD Workflows")
#xticks = range(0, 3500, 1000*20/40)
#xlabs  = [str(20*xt/1000) for xt in xticks]
#plt.xticks(xticks, xlabs)
plt.xlabel("Total Simulation Time ["+r"$\mu$"+"s]")
plt.ylabel("Total Probability of Explored States")
plt.legend(loc=2,prop={'size':12})
plt.savefig('equilibriumexploration.png', dpi=400)
#plt.savefig(clusterexplorationfile, dpi=400)
plt.close()



for nm,pars in plot_this.items():
    plt.plot( range(len(pars['n_clusters'])), pars['n_clusters'], label=nm )

plt.title( "Number of Clusters Explored with Different MD Workflows")
#xticks = range(0, 3500, 1000*20/40)
#xlabs  = [str(20*xt/1000) for xt in xticks]
#plt.xticks(xticks, xlabs)
plt.xlabel("Total Simulation Time ["+r"$\mu$"+"s]")
plt.ylabel("Number of Clusters Explored")
plt.legend(loc=2,prop={'size':12})
plt.savefig('clusterexploration.png', dpi=400)
#plt.savefig(clusterexplorationfile, dpi=400)
plt.close()




# in process_tools #     for nm, pars in datasets.items():
# in process_tools #         print nm
# in process_tools #         datashape = pars['datashape']
# in process_tools #         j    = -1
# in process_tools #         done = False
# in process_tools #         while not done:
# in process_tools #             j += 1
# in process_tools #             datashape.append(list())
# in process_tools #             try:
# in process_tools #                 for i in range(pars['n_trajs']):
# in process_tools #                     datashape[-1].append(mdtraj.load(
# in process_tools #                             pars['trajfiles'][ i + j*pars['n_trajs'] ],
# in process_tools #                             top=pars['topfile'],
# in process_tools #                             ).n_frames
# in process_tools #                             )
# in process_tools #             except IndexError:
# in process_tools #                 done = True

restarts = [350*i for i in range(1,8)]
for nm,pars in datasets.items():
    if nm == 'explore microstates':
        plt.plot(restarts,[datasets[nm]['clusters'][r] for r in restarts],'o', color='grey')
    if nm != 'single trajectory':
        if nm == 'explore microstates':
            lbl='adaptive sampling'
        else:
            lbl=nm
        plt.plot(
             range(len( datasets[nm]['clusters'] )),
             datasets[nm]['clusters'],
             label=lbl,
            )

plt.title( "Conformation Exploration with Different MD Workflows")

xticks = range(0, 3500, 1000*20/40)
xlabs  = [str(20*xt/1000) for xt in xticks]

plt.xticks(xticks, xlabs)

plt.xlabel("Total Simulation Time ["+r"$\mu$"+"s]")
plt.ylabel("Number of Clusters Explored")

plt.legend(loc=2,prop={'size':12})

plt.savefig(clusterexplorationfile, dpi=400)
plt.close()









