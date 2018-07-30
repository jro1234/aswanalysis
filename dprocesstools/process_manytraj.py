#!/usr/bin/env python


from __future__ import print_function
import itertools
from glob import glob # glob supports wildcard selections, very useful ;-)
import numpy as np
import pyemma
import pyemma.coordinates as coor
pyemma.version
pyemma.config.show_progress_bars=False
import analtools

# TODO CHECK if this is ok and delete
#            extra imports
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import time


printout = lambda s: print('{0} ||::|| {1}'.format(time.strftime("%Y-%m-%d %H:%M:%S"), s))

#workshopfolder = 'wsdata_full'
workshopfolder = 'wsdata_full2'


printout("INITIALIZING SETUP")

# NUBMER OF CLUSTERS
# ie microstates
kk   = 1000
n_TICs = 12
n_metastable_groups = 12

rowstring = lambda x: '  '.join(['{{{0}}}'.format(i) for i in range(x)])


doallplots = False
#lagidx = 8
lagidx = 0
indir = './'
topfile =  indir + 'NTL9-0-protein_fixed.pdb'
traj_list = sorted(glob(indir + 'NTL9-0-protein-*'))
#print(topfile)
feat = coor.featurizer(topfile)
feat.add_distances(feat.select_Ca())
#feat.add_contacts(feat.select('(rescode K or rescode R or rescode H) and mass > 13'),
#                  feat.select('(rescode D or rescode E) and mass > 13')
#                 )

#feat.add_backbone_torsions(cossin=True, periodic=True)
#feat.add_residue_mindist()
#feat.add_contacts(feat.select("not water and element H"), feat.select("not water and (element O or element N)"))
index_labels = feat.describe()
feat_labels = list()

for i,l in enumerate(index_labels):
    if l.startswith("DIST"):
        ll = ' - '.join([' '.join(l.split()[:-1])
                         for l in l.split('-')])
    elif l.startswith("RES_DIST"):
        ll = ':'.join(l.split(' (closest-heavy)'))
    else:
        ll = l
    ll += ' [{0}]'.format(i)
    feat_labels.append(ll)


fig, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(featuretrajs[0][:,155])
ax[1].plot(scorestraj[0][:,0])
plt.savefig(workshopfolder+'/stacked_feat_tic.png', dpi=600)
plt.close()


inp = coor.source([traj_list], feat, chunksize=1000)
print('number of trajectories = ', inp.number_of_trajectories())
print('trajectory length = ', inp.trajectory_length(0))
#print('trajectory time step = ', 500.0 / (inp.trajectory_length(0) - 1), 'ns')
print('number of dimension = ', inp.dimension())



featuretrajs = inp.get_output()
header = '--'.join(feat.describe())
rowsize = featuretrajs[0].shape[1]
format_template = ' '.join(['{{{0}:.3f}}'.format(i) for i in range(rowsize)])
#print(header)
#print(format_template)
#inp.write_to_csv('trajs.dat', header=header, delimiter='  ', )

printout("WRITING FEATURE TRAJ")
for i,ft in enumerate(featuretrajs):
    ftdata = ft[::1112]
    analtools.dump_data(workshopfolder+'/traj_{0}.dat'.format(i), ftdata, format_template, header)

# This used when coor.source(traj_list) has traj_list not embedded in outer list
#  ie each traj is treated individually instead of as a part of larger traj
#with open(workshopfolder+'/last-first.overlap.dat', 'w') as f:
#    overlaps = list()
#    for i, ft in enumerate(featuretrajs):
#        overlapts.append(
#            np.dot( ft[-1]                / np.linalg.norm(ft[-1]),
#                    featuretrajs[i+1][200]  / np.linalg.norm(featuretrajs[i+1][200])
#                  )
#    f.write(overlaps)

# TODO
# flatten lags to list of values
# lagtime = lag [step] * 0.2 [ns / step]
# lagtime(lag=300) = 60ns
lagslist = [100, 300, 800]# 50, 100, 200, 500, 1000]
strides  = [  5,  10,  30]#  4,   3,   3,   2,    2]
lags = dict(short = zip(strides, lagslist),
            #medium = [10*i for i in range(1,10)],
            #long =   [100*i for i in range(1,10)]
           )

ticass = dict()
## doesn't work since tica function isn't at top of pyemma module
#keys, vals = zip(*lags.items())
#data=[(inp,val) for val in np.concatenate(vals)]
#ticasresult = analtools.multiproc(pyemma.coordinates.tica, data)

for k,v in lags.items():
    ticass[k]          = dict()
    ticas              = list()
    ticass[k]['ticas'] = ticas
    for stride,lag in v:
        printout("CALCULATING TICA:  lag- {0}  stride-{1}".format(lag, stride))
        ticas.append(pyemma.coordinates.tica(inp, stride=stride, lag=lag, dim=n_TICs,))


printout("GETTING OUTPUT")
for lagsname, ticas in ticass.items():
    ticas.update({"ystacked_ticas": list()})
    #[ticas['ystacked_ticas'].append(tica)
    [ticas['ystacked_ticas'].append(tica.get_output())
    # this one cats the trajectories to make like 1 long one
    #[ticas['ystacked_ticas'].append(np.concatenate(tica.get_output()))
     for tica in ticas['ticas']]


printout("WE HAVE {} TOTAL TRAJECTORIES".format(len(ticas['ystacked_ticas'][0])))



printout("HISTOGRAMMING TICA TRAJ")
# TODO this needs to be moved below to avoid using huge ystacked_ticas
#      instead, create a sparse traj and ouput everything from that guy
scorestraj = ticass['short']['ystacked_ticas'][lagidx]
tica01_hist, tica0_coords, tica1_coords = np.histogram2d(
                                                         scorestraj[0][:,0],
                                                         scorestraj[0][:,1],
                                                         bins=1024)
                                                         #np.concatenate(scorestraj)[:,0],
                                                         #np.concatenate(scorestraj)[:,1],
                                                         #bins=64)


printout("FINDING BIN ENERGIES")
energy_shift = np.log(max([max(row) for row in tica01_hist]))
bin_energies = np.array(energy_shift + -np.log(tica01_hist))


tica0_coords[-1] += 1
tica0_coords[-2] = tica0_coords[-1] - 0.1
tica1_coords[-1] += 1
tica1_coords[-2] = tica1_coords[-1] - 0.1


printout("PLACING EACH FRAME")
framebins_0 = list()
for _traj in scorestraj:
    traj = _traj[::1112]
    print(len(traj))
    _framebins = list()
    for frame in traj:
        framebin = 0
        for tc in tica0_coords:
            if frame[0] > tc:
                framebin += 1
            else:
                break
        _framebins.append(framebin)
    framebins_0.append(np.array(_framebins))

    
framebins_1 = list()
for _traj in scorestraj:
    traj = _traj[::1112]
    print(len(traj))
    _framebins = list()
    for frame in traj:
        framebin = 0
        for tc in tica1_coords:
            if frame[1] > tc:
                framebin += 1
            else:
                break
        _framebins.append(framebin)
    framebins_1.append(np.array(_framebins))


printout("MAKING ENERGIES TRAJ")
framebins = np.dstack((framebins_0, framebins_1))
energies_traj = [[bin_energies[b1-1,b2-1]
                  if not np.isinf(bin_energies[b1-1,b2-1])
                  else energy_shift
                  for b1, b2 in fbs]
                 for fbs in framebins]


printout("WRITING ENERGY & TIC TRAJ")
trajdata = list()
for i, nrgtraj in enumerate(energies_traj):
    block = np.insert(scorestraj[i][::1112][:,:2], 0, nrgtraj, axis=1)
    trajdata.append(block)


header = '   '.join(['Energy', 'TIC0', 'TIC1', 'TIC2', 'TIC3'])
for i,td in enumerate(trajdata):
    analtools.dump_data(workshopfolder+'/traj_{0}.energy_tics.dat'.format(i),
                        td, 
                        format_template=rowstring(td.shape[1]),
                        header=header) 


printout("SAVING NRG TRAJ IN TICs PLOT")
fig, axes = plt.subplots(max(2,len(tictraj)), 1, figsize=(15,35), sharex=True)
for i,td in enumerate(trajdata):
    XY = td[:,1:3].T
    ########
    #CMAP ENERGY#
    #Z  = td[:,0]
    #label='energy'
    #cmap = plt.cm.plasma_r
    #CMAP ENERGY#
    ########
    ########
    #CMAP MICROSTATES#
    sortedmicrostates = tt['sortedmicrostates'][0]
    colorcoord, microstatesorderarray = np.array(zip(*sortedmicrostates))
    dt = microstates.dtrajs[0][::1112]
    Z = colorcoord[microstatesorderarray[dt].astype(int)]
    label = 'microstates'
    cmap = plt.cm.gist_ncar
    #CMAP MICROSTATES#
    ########
    ax = axes[i]
    s = ax.scatter(*XY, s=1, c=Z, cmap=cmap)
    ax.set_title('Lag: {0} Trajectory {1} Projected into TICs'.format(lags[lagsname][j],i+1))
    ax.set_xlabel('Coordinate in TIC1')
    ax.set_ylabel('Coordinate in TIC2')
    cb = fig.colorbar(s, ax=ax)
    cb.set_label('Frame Energy')

#plt.savefig(workshopfolder+'/traj_{0}_TICs.png'.format(label), dpi=600)
plt.savefig(workshopfolder+'/traj_{0}_TICs.png'.format(label), dpi=300)
plt.close()

####for i,td in enumerate(trajdata):
####    XY = td[:,1:3].T
####    sortedmicrostates = tt['sortedmicrostates'][0]
####    microstatesorderarray = np.array(zip(*sortedmicrostates)[1])
####    dt = microstates.dtrajs[0][::1112]
####    Z = dt#microstatesorderarray[dt]
####    label = 'microstates'
####    #s = ax.scatter(*XY, s=1, c=Z, cmap=plt.cm.plasma_r)
####    #cb = fig.colorbar(s, ax=ax)
####    #cb.set_label('Frame Energy')
####    plt.plot(*XY, marker='o', alpha=0)
####    for i, (x,y) in enumerate(XY.T):
####        plt.text(x,y,Z[i],fontsize=1)
####    #plt.set_title('Lag: {0} Trajectory {1} Projected into TICs'.format(lags[lagsname][j],i+1))
####    #plt.set_xlabel('Coordinate in TIC1')
####    #plt.set_ylabel('Coordinate in TIC2')
####
####plt.savefig(workshopfolder+'/traj_{0}_TICs.png'.format(label), dpi=800)
####plt.close()





for lagsname,tt in ticass.items():
    for tica in tt['ticas']:
        printout("WRITING TICA MATRIX:  lag- {0}  stride-{1}".format(tica.lag, tica.stride))
        #tica.eigenvalues
        #tica.eigenvectors
        #rowstring(tica.eigenvalues.shape[0])
        #tica.eigenvectors[0]
        with open(workshopfolder+'/eigendata.lag{0}.dat'.format(tica.lag), 'w') as out:
            out.write(rowstring(tica.eigenvalues.shape[0]).format(*tica.eigenvalues)+'\n\n')
            for vec in tica.eigenvectors:
                out.write(rowstring(vec.shape[0]).format(*vec)+'\n')


printout("CLUSTERING NOW")
for lagsname, tt in ticass.items():
    microstates       = list()
    tt['microstates'] = microstates
    printout("CLUSTERING:  k- {0}  stride-{1}".format(kk, stride))
    for tica in tt['ticas']:
        printout("CLUSTERING FROM TICAS:  lag- {0}  stride-{1}".format(tica.lag, tica.stride))
        microstates.append(coor.cluster_kmeans(tica, stride=stride, k=kk, max_iter=1000))


printout("ESTIMATING MARKOV MODEL")
for lagsname, tt in ticass.items():
    msms       = list()
    tt['msms'] = msms
    for i,microstates in enumerate(tt['microstates']):
        lag = lags[lagsname][i][1]
        printout("MSM with:  lag- {0}".format(lag))
        msms.append(pyemma.msm.estimate_markov_model(microstates.dtrajs, lag=lag))


#printout("COUNTING MICROSTATES DATA POPULATIONS")
#microstates=ticass[lagsname]['microstates'][i]
#sampledmicrostates = np.zeros(microstates.n_clusters)
#for dtraj in microstates.dtrajs:
#    for cluster in dtraj:
#        sampledmicrostates[cluster] += 1







for lagsname, tt in ticass.items():
    analtools.plot_msm_timescales(tt['msms'], plt)

plt.savefig(workshopfolder+'/modeltimescales.png')
plt.close()




for lagsname, tt in ticass.items():
    lagslist = [50, 100, 250, 500, 1000]
    for i,microstates in enumerate(tt['microstates']):
        analtools.plot_msm_its(microstates, plt, lags[lagsname][i], lagslist)
        plt.savefig(workshopfolder+'/impliedtimescales.{0}.png'.format(i))
        plt.close()








# TODO FIXME old iteration formatting starting here
# TODO FIXME old iteration formatting starting here
# TODO FIXME old iteration formatting starting here


printout("SORTING MICROSTATES ALING TIC1")
for lagsname, tt in ticass.items():
    sortedmicrostates       = list()
    tt['sortedmicrostates'] = sortedmicrostates
    for microstates in tt['microstates']:
        sortedmicrostates.append(sorted(zip(microstates.clustercenters[:,0],
                               range(microstates.clustercenters.shape[0])),
                           key=lambda x: x[0]
                          ))


printout("WRITING MICROSTATE CLUSTER TRAJECTORY")
for lagnsname,tt in ticass.items():
    for i,microstates in enumerate(tt['microstates']):
        sortedmicrostates = tt['sortedmicrostates'][i]
        filename = workshopfolder+'/traj_{0}.{1}.microclusters.dat'.format(i,lagsname)
        microstatesorderarray = np.array(zip(*sortedmicrostates)[1])
        with open(filename, 'w') as out:
            dt = microstates.dtrajs[0][::1112]
            out.write(rowstring(dt.shape[0]).format(*microstatesorderarray[dt]))



#printout("SAVING EIGENVALUES PLOT")
#lagsname = 'short'
#fig, axes = plt.subplots(1,3, figsize=(20,10))
#for i, (lagsname, ticas) in enumerate(ticass.items()):
#    #analtools.multiproc('plot_eigenvalues', ticas, axes[i])
#    analtools.plot_eigenvalues(ticas['ticas'], plt, axes[i],)
#    plt.savefig(workshopfolder+'/eigenvalues.png', dpi=600)
#    plt.close()



printout("SAVING TRAJ IN TICs PLOT")
for lagsname,tt in ticass.items():
    for j,tica in enumerate(tt['ticas'][:1]):
        tictraj = tica.get_output()
        fig, axes = plt.subplots(max(2,len(tictraj)), 1, figsize=(15,35), sharex=True)
        for i,tic in enumerate(tictraj):
            XY = tic[:,:2].T[:,:5001]
            ax = axes[i]
            s = ax.scatter(*XY, s=1, c=range(XY.shape[1]), cmap=plt.cm.plasma)
            ax.set_title('Lag: {0} Trajectory {1} Projected into TICs'.format(lags[lagsname][j],i+1))
            ax.set_xlabel('Coordinate in TIC1')
            ax.set_ylabel('Coordinate in TIC2')
            cb = fig.colorbar(s, ax=ax)
            cb.set_label('Frame Number')

plt.savefig(workshopfolder+'/traj_inTICs.png', dpi=600)
plt.close()



###printout("SAVING FEATURE CORRELATIONS PLOT")
#### TODO SET the ylimits to avoid extra space on top
#### TODO WHY features at the bottom without correlations
###for lagsname,tt in ticass.items():
###    for j,tica in enumerate(ticas['ticas']):
###        fig, ax = plt.subplots(1,1, figsize=(30,180))
###        n_tics = tica.feature_TIC_correlation.shape[1]
###        reordering = np.argsort(abs(tica.feature_TIC_correlation.T[0]))
###        for i,corr in enumerate(tica.feature_TIC_correlation.T[:,reordering]):
###            plt.plot(
###            abs(corr),
###            range(len(corr)),
###            linestyle = '',
###            marker    = 'o',
###            markersize= 5,
###            color     = plt.cm.gnuplot2_r(21*(n_tics-i)),
###            label     = 'TIC{0}'.format(i)
###            )
###        [ax.axhline(0.5+offset, c=str(0.85+0.08*(offset%2))) for offset in range(-1,len(feat_labels))]
###        ax.legend(loc='upper right')
###        ax.set_title("Lag {0} TIC Correlations with Features".format(lags[lagsname][j]))
###        ax.set_ylabel("Feature")
###        ax.set_xlabel("Absolute Correlation")
###        ax.set_yticks(range(len(feat_labels)))
###        ax.set_yticklabels(np.array(feat_labels)[reordering], fontsize=10)
###        plt.savefig(workshopfolder+'/feature_correlations.png', dpi=50)
###        plt.savefig(workshopfolder+'/feature_correlations.highres.png', dpi=200)
###        plt.close()



#def grouper(n, iterable):
#    it = iter(iterable)
#    while True:
#       chunk = tuple(itertools.islice(it, n))
#       if not chunk:
#           return None
#       yield chunk
#
#datasets = dict(anton=grouper(increment, iter(np.concatenate(ticass['short']['microstates'][lagidx].dtrajs))))
#visitations = dict()
#[visitations.update({k: set()}) for k in datasets]
#increment = 1000
#done = False
#while not done:
#    for nm,mstraj in datasets.items():
#        chunk = next(mstraj)
#        if chunk:
#            visitations[nm].append(unique
#        else:
#            




printout("CALCULATING BayesianMSM")
dtraj = ticass[lagsname]['microstates'][0].dtrajs
lag = lags[lagsname][0][1]
bmsm=pyemma.msm.bayesian_markov_model(dtraj, lag=lag)


printout("PREPPING reordered matrix")
cmap = "winter"
microstatesorder = list(zip(*sortedmicrostates))[1]
usematrix = bmsm.transition_matrix
newT = analtools.reorder_matrix(np.array(usematrix, copy=True),
                                microstatesorder, toreorder=True)



printout("PLOTTING Transition Matrix")
analtools.plot_heatmap_logvals(np.array(usematrix, copy=True),
                               plt, logmin=0.000001, cmap=cmap,
                               ylabels=list(range(microstates.n_clusters)))

plt.savefig(workshopfolder+'/unordered_T.png', dpi=600)
plt.close()



printout("PLOTTING reordered Transition Matrix")
analtools.plot_heatmap_logvals(newT, plt, logmin=0.000001,
                     cmap=cmap, ylabels=microstatesorder)

plt.savefig(workshopfolder+'/ordered_T.png', dpi=600)
plt.close()



#printout("PLOTTING sorted reordered Transition Matrix")
#analtools.plot_heatmap_logvals(newT, plt, logmin=0.000001,
#                     cmap=cmap, ylabels=microstatesorder, sort=True)

#plt.savefig(workshopfolder+'/ordered_resorted_T.png', dpi=600)
#plt.close()



printout("PLOTTING Stationary Distribution")
plt.plot(bmsm.stationary_distribution)
plt.savefig(workshopfolder+'/stationary_distribution.png', dpi=600)
plt.close()




printout("CALCULATING PCCA Coarse MSM")
# How to choose n_PCCA?
bpcca = bmsm.pcca(n_metastable_groups)

printout("PLOTTING Coarse Transition Matrix")
plt.imshow(bpcca.coarse_grained_transition_matrix)
plt.colorbar()
plt.savefig(workshopfolder+'/pcca_T.png', dpi=600)
plt.close()

macrostatemapping = list()
for i,micros in enumerate(bpcca.metastable_sets):
     for j in micros:
            macrostatemapping.append(tuple([j,i]))

macrostatemapping = np.array(zip(*sorted(macrostatemapping))[1])


printout("WRITING MACROCLUSTER Trajectory")
for i, dtraj in enumerate(ticass['short']['microstates'][lagidx].dtrajs):
    filename = workshopfolder+'/traj_{0}.macroclusters.dat'.format(i)
    with open(filename, 'w') as out:
        dt = dtraj[::1112]
        out.write(rowstring(dt.shape[0]).format(*macrostatemapping[dt]))

pccaorder = np.concatenate(bpcca.metastable_sets)
usematrix = bmsm.transition_matrix
printout("PREPPING PCCAordered Transition Matrix")
newT2 = analtools.reorder_matrix(np.array(usematrix, copy=True),
                                 pccaorder, toreorder=True)


printout("PLOTTING PCCAordered Transition Matrix")
analtools.plot_heatmap_logvals(newT2, plt, logmin=0.000001,
                     cmap='plasma', ylabels=list(pccaorder),
                              )#sort=True)

plt.savefig(workshopfolder+'/pccaordered_T.png', dpi=600)
plt.close()



printout("PLOTTING Timescales")
timescales = bmsm.timescales()
plt.plot(timescales[:75], marker='.')
plt.savefig(workshopfolder+'/timescales.png', dpi=600)
plt.close()


printout("PLOTTING PCCA Stationary Distribution")
plt.plot(bpcca.coarse_grained_stationary_probability)
plt.savefig(workshopfolder+'/stationary_distribution.macrostates.png', dpi=600)
plt.close()



printout("SAVING MSM Transition Matrix")
np.savetxt(workshopfolder+'/transitionmatrix.microstates.dat', newT2)
np.savetxt(workshopfolder+'/stationarydistribution.microstates.dat', analtools.reorder_matrix(bmsm.stationary_distribution,pccaorder,True))

printout("SAVING PCCA Transition Matrix")
np.savetxt(workshopfolder+'/transitionmatrix.macrostates.dat', bpcca.coarse_grained_transition_matrix)
printout("SAVING Stationary Distribution")
np.savetxt(workshopfolder+'/stationarydistribution.macrostates.dat', bpcca.coarse_grained_stationary_probability)


P_nonrev = bpcca.coarse_grained_transition_matrix
w_nonrev = bpcca.coarse_grained_stationary_probability
# make reversible
C = np.dot(np.diag(w_nonrev), P_nonrev)
Csym = C + C.T
P = Csym / np.sum(Csym,axis=1)[:,np.newaxis]
# CHEATING HERE!!
Pan = np.array([pa/np.sum(pa) for pa in abs(P)])
#Pan[ Pan < 0.001 ] = 0
#Pan = np.array([pa/np.sum(pa) for pa in Pan])

pyemma.plots.plot_markov_model(Pan, state_sizes=w_nonrev)
plt.savefig("maybemarkovmodel_2.png")
plt.close()


msss = bpcca.metastable_sets
spss = np.array([sum(bmsm.stationary_distribution[mss])**0.25 for mss in msss])
spss /= np.sum(spss)
np.savetxt(workshopfolder+'/stationarydistribution.macrostates.dat', bpcca.coarse_grained_stationary_probability)


tps = np.array([[np.sum(bpcca.transition_matrix[i][:,j])
                 if i[0] not in j else np.sum(bpcca.transition_matrix[i][:,j])**0.5
                 for i in msss] for j in msss])


for tp in tps:
    tp /= np.sum(tp)

pyemma.plots.plot_markov_model(tps, state_sizes=spss)
plt.savefig(workshopfolder+'/network.png', dpi=600)
plt.close()

np.savetxt(workshopfolder+'/stationarydistribution.macrostates.dat', spss)
np.savetxt(workshopfolder+'/transitionmatrix.macrostates.dat', tps)

