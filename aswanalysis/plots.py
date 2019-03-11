'''
Some functions to create the trickier plots of Adaptive Sampling
Workflow data, along with some ready-to-go formatting like
color sets.
'''

import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pyemma.plots as pymplots

import numpy as np


colors = dict()
colors['tennessee'] = tncolors = dict()
tncolors['smokey'] = [rgb/255. for rgb in [88, 89, 91]]
tncolors['orange'] = [rgb/255. for rgb in [255, 130, 0]]
tncolors['white'] = [rgb/255. for rgb in [0, 0, 0]]
tncolors['eureka'] = [rgb/255. for rgb in [235, 234, 100]]

dssp_colormap = {k:v for k,v in zip(
    ['H','E','C'],
    map(lambda c: (c,tuple(tncolors[c])),
    filter(lambda c: c != 'white', tncolors))
)}


def landscape_creator(XY,):#title='', xlabel='TIC0', ylabel='TIC1'):
    #return lambda kwargs: _landscape_plt(XY, title, xlabel, ylabel, **kwargs)
    return lambda kwargs: _landscape_plot(XY, **kwargs)


def _landscape_plot(XY, title='', xlabel='TIC0', ylabel='TIC1'):
    fig, ax = pymplots.plot_free_energy(*XY, cmap='hot')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return (fig, ax)


def plot_dssp(dssp, filepath, colormap=dssp_colormap):
    # FIXME aspect ratio control is required
    # TODO set ylabels based on frame timestep
    # 
    # dssp data structure
    # 
    #  np array with row x col: n_frames x n_residues
    #    with dtype 1 character "S1"
    # 
    ss_names = {"H":"helix", "E":"sheet", "C":"loop"}
    assert all(ss_nm in list(colormap) for ss_nm in list(ss_names))
    C = np.empty(list(dssp.shape)+[3], dtype=float)
    c = list()
    for k,v in colormap.items():
        _c = v[1]
        C[ dssp == k ] = _c
        c.append(
          (plt.Line2D([0], [0], color=_c, lw=8),
          ss_names[k])
        )
    fig, ax = plt.subplots(figsize=(30,10))
    ax.imshow(C)
    ax.set_aspect('equal')
    #plt.xlabel('Time [ns]')
    ax.set_xlabel('Frame Index')
    ax.set_ylabel('Residue Index')
    ax.set_title('Secondary Structure over Trajectory')
    ax.legend(*zip(*c))
  ##  plt.imshow(C)
  ##  plt.xlabel('Frame Index')
  ##  plt.ylabel('Residue Index')
  ##  plt.suptitle('Secondary Structure over Trajectory')
  ##  plt.legend(*zip(*c))
    plt.savefig(filepath, dpi=600)
    plt.close()
  

def plot_cumvar(filepath, tica):
    plt.plot(range(1, 1+tica.cumvar.shape[0]),
        tica.cumvar, linestyle='', marker='o')
    plt.ylim([0,1.05])
    plt.title("Cumulative Variance Explained by TICs")
    plt.xlabel("TIC Index")
    plt.ylabel("Cumulative Variance")

    if filepath.endswith('png'):
        plt.savefig(filepath, dpi=600)
    else:
        plt.savefig(filepath)

    plt.close()


def plot_clust_explored(filepath, dataset, replicates=dict(), epoch_coeff=1):
    assert all([all([r in dataset for r in rep]) for rep in replicates.values()])
    _repmap = dict()
    for nm,rep in replicates.items():
        for r in rep:
            _repmap[r] = nm
    for nm in filter(lambda nm: nm not in _repmap, dataset):
        _repmap[nm] = nm
        replicates[nm] = {nm}
    for nm,rep in replicates.items():
        _ds = dataset[next(iter(rep))] # assuming same shape for all replicates!
        X = np.array([i*_ds['epochsize']*_ds['n_trajs']*_ds['timestep']*epoch_coeff
             for i in range(len(_ds['clust_explored']))])
        V = np.array([[(100.*nc)/_ds['n_clusters'] for nc in dataset[r]['clust_explored']] for r in rep])
        Y = np.mean(V, axis=0)
        plt.plot(X, Y, label=nm)
        if len(V.shape) > 1:
            Yerr = np.std(V, axis=0)
            plt.fill_between(X, Y-Yerr, Y+Yerr, alpha=0.3)
    #plt.title( "Equilibrium Distribution Explored with Different MD Workflows")
    # Hope that use last is OK
    xticks = [i*_ds['epochsize']*_ds['n_trajs']*_ds['timestep']*epoch_coeff for i in range(len(_ds['clust_explored'])) if not bool(i%50)]
    xlabs  = [str(int(xt/1000)) for xt in xticks]
    plt.ylim = [0, 105] # The n_clusters should be same for all!
    plt.axhline(100, c='k', ls=':')
    plt.xticks(xticks, xlabs)
    plt.xlabel("Total Simulation Time ["+r"$\mu$"+"s]")
    plt.ylabel("Percentage of States Visited")
    plt.legend(loc='best',prop={'size':12})
    plt.savefig(filepath, dpi=400)
    #plt.savefig(clusterexplorationfile, dpi=400)
    plt.close()

def plot_eq_explored(filepath, dataset, replicates=dict(), epoch_coeff=1):
    assert all([all([r in dataset for r in rep]) for rep in replicates.values()])
    _repmap = dict()
    for nm,rep in replicates.items():
        for r in rep:
            _repmap[r] = nm
    for nm in filter(lambda nm: nm not in _repmap, dataset):
        _repmap[nm] = nm
        replicates[nm] = {nm}
    for nm,rep in replicates.items():
        _ds = dataset[next(iter(rep))] # assuming same shape for all replicates!
        X = np.array([i*_ds['epochsize']*_ds['n_trajs']*_ds['timestep']*epoch_coeff
             for i in range(len(_ds['eq_explored']))])
        V = np.array([dataset[r]['eq_explored'] for r in rep])
        Y = np.mean(V, axis=0)
        plt.plot(X, Y, label=nm)
        if len(V.shape) > 1:
            Yerr = np.std(V, axis=0)
            plt.fill_between(X, Y-Yerr, Y+Yerr, alpha=0.3)
    #plt.title( "Equilibrium Distribution Explored with Different MD Workflows")
    # Hope that use last is OK
    xticks = [i*_ds['epochsize']*_ds['n_trajs']*_ds['timestep']*epoch_coeff for i in range(len(_ds['eq_explored'])) if not bool(i%50)]
    xlabs  = [str(int(xt/1000)) for xt in xticks]
    plt.ylim = [0,1.1]
    plt.axhline(1, c='k', ls=':')
    plt.xticks(xticks, xlabs)
    plt.xlabel("Total Simulation Time ["+r"$\mu$"+"s]")
    plt.ylabel("Total Probability of Visited States")
    plt.legend(loc='best',prop={'size':12})
    plt.savefig(filepath, dpi=400)
    #plt.savefig(clusterexplorationfile, dpi=400)
    plt.close()


def plot_eq_explored_walltime():
    for nm,pars in dataset.items():
        #plt.plot(range(len(pars['eq_explored'])), pars['eq_explored'], label='Epoch count' )
        plt.plot([i*pars['epochsize']*pars['timestep']*epoch_coeff
                  for i in range(len(pars['eq_explored']))],
                  pars['eq_explored'], label=nm )
    #plt.title( "Equilibrium Distribution Explored with Different MD Workflows")
    # Use last is OK
    xticks = [i*10*170 for i in range(6)] # 0, 10, 20, etc. days
    xlabs  = [str(i*10) for i in range(6)]
    plt    . xticks(xticks, xlabs)
    plt    . xlabel("Total Wall Time [days]")
    plt    . ylabel("Total Probability of Visited States")
    plt    . legend(loc='best',prop={'size':12})
    fp     = filepath.split('.')
    fp[-2] += '-walltime'
    plt    . savefig('.'.join(fp), dpi=400)
    plt    . close()


def plot_feature_correlations(filename, tica, feat_labels=None):
    # FIXME the colormap isn't working very nicely
    fig, ax = plt.subplots(1,1, figsize=(15,60))
    n_tics  = tica.feature_TIC_correlation.shape[1]
    sortkey = np.argsort(tica.feature_TIC_correlation[:,0])
    for i,corr in enumerate(tica.feature_TIC_correlation[sortkey].T):
        plt.plot(
            abs(corr),
            range(len(corr)),
            linestyle = '',
            marker    = 'o',
            #markersize= 0.5,
            color     = plt.cm.cool(75*(n_tics-i)),
            label     = 'TIC{0}'.format(i)
        )
    if not feat_labels:
        feat_labels = [str(i) for i in sortkey]
    else:
        feat_labels = [feat_labels[i] for i in sortkey]
    [ax.axhline(0.5+offset, c=str(0.85+0.08*(offset%2))) for offset in range(-1,len(feat_labels))]
    ax.legend(loc='upper right')
    ax.set_title("TIC Correlations with Features")
    ax.set_ylabel("Feature")
    ax.set_xlabel("Absolute Correlation")
    ax.set_yticks(range(len(feat_labels)))
    ax.set_yticklabels(feat_labels)
    plt.savefig(filename, dpi=500)
    plt.close()



