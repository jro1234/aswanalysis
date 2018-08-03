#from __future__ import print_function
import mdtraj
import subprocess, shlex
import os
from glob import glob

import pyemma.coordinates as coor
import pyemma.msm as msm
import numpy as np


# TODO use logger 


fmt_spacedlist = lambda l: ' '.join(['{}'.format(ll) for ll in l])

min_full_datas = lambda dts: min([np.sum(np.sum(ds['datashape']))/ds['stride'] for nm,ds in dts.items()])
onlydatas = lambda datasets: {k:v for k,v in datasets.items() if k is not 'analysis'}



def squish_tica_inputfiles(datasets):
    ticainputs = list()
    # The 'dim_reduction' dict must only have parameters for
    # the function call, store 'input_order' separately
    input_order = list()
    #datasets['analysis']['dim_reduction']['input_order'] = list()
    for nm,pars in datasets.items():
        ticainputs.extend(pars.get('tica_inputs', []))
        input_order.append(nm)
    return ticainputs, input_order


def prepare_tica_inputs(datasets, topfile, get_output=False):
    assert isinstance(get_output, bool)
    # FIXME This isn't very general... use adaptivemd.analysis.apply_feat_part?
    # PyEMMA data prep
    feat = coor.featurizer(topfile)
    feat.add_inverse_distances(feat.select_Ca())
    #feat.add_distances(feat.select_Ca())
    ticainputs,input_order = squish_tica_inputfiles(datasets)
    #datasets['analysis']['dim_reduction']['input_order'].append(nm)
    tica_inp     = coor.source(ticainputs, feat, chunksize=5000)
    if get_output:
        tica_inp = tica_inp.get_output()
    return tica_inp, input_order


def analyze(d, a):
    #process_tools.create_tica_inputs(d, namekey)
    tica_inp, input_order = prepare_tica_inputs(d, a['topfile'])
    tica_func    = build_call(a['dim_reduction'], coor, tica_inp)
    tica_results = tica_func()
    cluster      = build_call(a['clustering'],    coor, tica_results.get_output())
    clustering   = cluster()
    msm_inp      = clustering.dtrajs
    emm          = build_call(a['msm'],           msm,  msm_inp)
    mm           = emm()
    bind_dtraj_results(d, input_order, clustering.dtrajs)
    return tica_results, clustering, mm


def prepare_data(dataset, analysis, topdir):
    # NOW starting to process data
    create_trajlist     (dataset, topdir)
    assign_stride       (dataset)# add namekey
    determine_epochsize (dataset)
    determine_datashape (dataset, analysis)


def bind_dtraj_results(datasets, input_order, dtrajs):
    # This tries to reference dtrajs correctly from dataset dicts
    n_trajs_bound = 0
    for i,nm in enumerate(input_order):
        pars = datasets[nm]
        if pars.get('trajfiles', []):
            these          = len(pars['trajfiles'])
            pars['dtrajs'] = dtrajs[n_trajs_bound:n_trajs_bound+these]
            n_trajs_bound += these


def calc_eq_explored(dataset, mm):
    # We're going to get the total equilibrium population explored
    # for each dataset with consistent epochsizes across datasets
    # ie more parallel <--> smaller chunks per traj
    # i ~ number of epochs
    nfc = n_frames_compare(dataset)
    #for nm,pars in process_tools.onlydatas(datasets).items():
    for nm,pars in dataset.items():
        pars['n_clusters']  = list()
        pars['eq_explored'] = list()
        pars['clusters']    = list()
        # _array rows smoosh workloads to get data in row-pll count order
        _array = [[] for _ in range(pars['n_trajs'])]
        [_array[i].extend(np.concatenate(pars['dtrajs'][i::pars['n_trajs']])) for i in range(pars['n_trajs'])]
        # FIXME if dtrajs are inhomogenous then need to end with sequence
        # of last cluster to avoid extra clusters and make array
        dta = np.array(_array)
        i   =  0
        # while total_counted_frames < frames_to_compare:
        #while i*pars['epochsize']*pars['n_trajs'] < 1000:
        while i*pars['epochsize']/pars['stride']*pars['n_trajs'] < nfc:
            try:
                pars["clusters"].append(np.unique(
                        # this only looks at unique cluster IDs in the epoch
                        dta[:,i*pars['epochsize']/pars['stride']:(i+1)*pars['epochsize']/pars['stride']]
                        # this looks at unique cluster IDs up to & incl this epoch
                        #dta[:,:(i+1)*pars['epochsize']]
                ))
                # this collects unique cluster IDs found in all epochs
                exploredclusters = np.unique(np.concatenate(pars["clusters"]))
                pars["n_clusters"] .append(len(exploredclusters))
                pars["eq_explored"].append(np.sum(mm.stationary_distribution[exploredclusters]))
            except IndexError:
                print "Can't get uneven final chunk"
                print (i+1)*pars['epochsize']*j*pars['n_trajs'], "<", nfc
            i += 1

def n_frames_compare(datasets):
    mfd = min_full_datas(onlydatas(datasets))
    return mfd

def assign_trajs_to_clusters(datasets, MORE):
    inp_files = min_full_datas(datasets)
    for nm, pars in datasets.items():
        #if nm is not master_dataset:
        #if nm is not 'single trajectory':
        print 'Projecting {} into clustercenters'.format(nm)
        feat = coor.featurizer(pars['topfile'])
        feat.add_distances(feat.select_Ca())
        inp = coor.source(inp_files[nm], feat, chunksize=1000)
        dtrajs = coor.assign_to_centers(
                inp,
                centers=all_clust.clustercenters,
                stride=pars['stride'],
                metric='minRMSD',
                chunksize=1000
                )
        pars['feat']   = feat
        pars['inp']    = inp
        pars['dtrajs'] = dtrajs
        np.save(
                topdir + nm + '.npy',
                disc.dtrajs
               )



def write_pdb_topfile(filename, trajectory, atomselection='all', framenumber=0):
    trajectory.atom_slice(atomselection)[framenumber].save(filename)


def assign_stride(datasets):
    '''
    This function determines stride through datasets
    to match frame rate in simulation time:
     eg 
        dataset 1: 200 ps frames  --> stride 1
        dataset 2:  50 ps frames  --> stride 4
    '''
    longest = max([pars['timestep'] for nm,pars in datasets.items() if 'directory' in pars])
    for nm, pars in onlydatas(datasets).items():
        pars['stride'] = int(longest/pars['timestep'])



def create_trajlist(datasets, topdir='.'):
    '''
    Argument is a dict of `dict` dataset definitions.
    Each definition must have must have fields:
      - directory - relative home directory of trajectory files
      - filename  - pattern to match to pull the trajectory files
      - trajlist  - `list` for storing paths to trajectory files
    '''
    for nm, pars in onlydatas(datasets).items():
        filenames = os.path.join(topdir, pars['directory'], pars['filename'])
        first = pars.get('first', 0)
        last  = pars.get('last', 99999999)
        pars['trajfiles'] = sorted(glob(filenames))[first:last]


# TODO FIXME this probably misses trajectories if they are under
#            a dataset directory but no further subdirectory...
#             - ie the '**' implies a folder must exist
def create_tica_inputs(datasets, namekey, topdir='.'):
    for nm, pars in datasets.items():
        if 'directory' not in pars:
            continue
        middle_part = ''
        if pars['directory']:
            middle_part = os.path.join(pars['directory'], '**')
        filenames = os.path.join(topdir, middle_part, '*'+namekey+'*.dcd')
        print "USING this search key: ", filenames
        pars['tica_inputs'] = sorted(glob(filenames))
        print "found these guyts: ", pars['tica_inputs']


def build_call(calldict, parent, *args):
    '''
    >>> This is a dangerous but convenient function
    Build a call signature for function/method of parent
    given by calldict top key. 
    len(calldict) === 1
    Parent can be module or object with method of name

    Arguments
    ---------
    '''
    if len(calldict) != 1:
        return None
    fname = calldict.keys()[0]
    f = getattr(parent, fname)
    #print f
    kwargs = calldict[fname]
    #print args
    #print kwargs
    return lambda: f(*args, **kwargs)


# THIS IS SUPER SLOW
# AND UNNECESSARILY DUPLICATES THE STRIDE 1 FILES
# Can't figure out how to access DCD header via mdtraj...
# ... no other way to get ahold of length without full read
#  - a better chunksize for big-node processing would be like 100,000
# TODO   datashape isn't used anymore --> except for total frame counting
#      FIXME   just make feature trajs? -> ?^^^?
def determine_datashape(datasets, analysis, topdir='.', prepare=True, preptag='--', chunksize=4000, maxchunks=False, atomselection="all"):
    '''
    This function provides some way of creating (as `list`)
    a shape for each dataset that is used to navigate equal
    size epochs. Creates a sub-`list` for each workload,
    with each entry corresponding to a trajectory length
    from that workload. The trajs in a workload are considered
    degenerate, but the workloads are in a pretermined order
    given by their order-of-appearance in the trajlist.

    If `prepare` is `True`, save trajectories in correct
    stride for the analysis--> add to `tica_inputs`

     - Without correct n_trajs (workload size), datashape will mix workloads
     - A (final) workload of different size will be missed entirely via the `j` indexing strategy

    '''
    all_topofile           = analysis['topfile']
    analysis_atomselection = analysis['atomselection']
    if all_topofile:
        all_topofile = os.path.join(topdir, all_topofile)
    if analysis_atomselection:
        atomselection = analysis_atomselection
    for nm,pars in onlydatas(datasets).items():
        if 'tica_inputs' not in pars:
            pars.update({'tica_inputs': list()})
        datashape = pars['datashape']
        # TODO replace `j` with workload
        j    = -1
        done = False
        if prepare:
            f_topo_atomselection = os.path.join(topdir, pars['directory'], pars['topfile'].split('.pdb')[0]+'.features.idx')
            print "\n --- Now using this feature topology file ", f_topo_atomselection
        while not done:
            j += 1
            try:
                #print "on j=",j
                print "\nOn set ", nm
                lengths = list()
                for i in range(pars['n_trajs']):
                    tfn = pars['trajfiles'][ i + j*pars['n_trajs'] ]  # traj filename
                    tfi = mdtraj.iterload(tfn, top=os.path.join(pars['directory'],pars['topfile']), chunk=chunksize)
                    tot = 0  # total frames in traj
                    tfs = '' # trajectory file atomselected
                    adx = None
                    if prepare:
                        # # TODO implement offset for proper stitching
                        # #      when mixing chunks and stride
                        # offset = 0
                        # tfs = "{0}--stride-{1}{{0}}.dcd".format(
                        #           tfn.split('dcd')[0], pars['stride'])
                        # TODO replace "stride" with preptag
                        tfs = "{0}--stride-{1}.dcd".format(
                                tfn.split('.dcd')[0], pars['stride'])
                    #print "\nlooking at this traj ", tfn
                    #print "using this traj loader ", tfi
                    for i,segment in enumerate(tfi):
                        print segment
                        if i == 0:
                            adx = segment.topology.select(atomselection)
                        if all_topofile:
                            if not os.path.exists(all_topofile):
                                print "Going to write this topo",segment.atom_slice(adx).topology
                                segment.atom_slice(adx)[0].save(all_topofile)
                        if maxchunks:
                            if i >= maxchunks:
                                break
                        tot += segment.n_frames
                       # if tfs:
                       #     #segment[offset::pars['stride']].save(tfs)
                       #     #offset = ( segment.n_frames + offset ) % stride + 1
                       #     segment[::pars['stride']].save(tfs.format(i))
                    lengths.append(tot)
                    print "added this many frames traj to datashape ", str(tot)
                    #mdtraj.join([tfs.format(k) for k in range(i)])
                    # dependent on `prepare` and file existence
                    if prepare:
                        if prepare and not os.path.exists(f_topo_atomselection):
                            print "Creating new atom index file: ", f_topo_atomselection
                            print "to print these indices: {}".format(adx)
                            with open(f_topo_atomselection, 'w') as f:
                                f.write(fmt_spacedlist(adx))
                        else:
                            print "Found existing atoms index file: ", f_topo_atomselection
                        if not os.path.exists(tfs):
                            print "using this file to store output atomselection ", tfs
                            subprocess.call(shlex.split("mdconvert -o {output_str} -a {selection} -s {stride} {input_str} ".format(output_str=tfs, stride=pars['stride'], input_str=tfn, selection=f_topo_atomselection)))
                        else:
                            print "Found already existing output file: ", tfs
                        pars['tica_inputs'].append(tfs)
                datashape.append(lengths)
            except IndexError:
                print "Done reading from set ", nm, "\n"
                done = True


def determine_epochsize(datasets):
    '''
    The minimum matching epochs are calculated
    considering the level of parallelism and in the
    native stride of the dataset, then we calculate
    the number of frames from each dataset in the
    native strides.
    '''
    epochsize = 0
    for nm,pars in onlydatas(datasets).items():
        e = pars['n_trajs'] * pars['timestep'] * pars['stride']
        if e > epochsize:
            epochsize = e
    for nm,pars in onlydatas(datasets).items():
        pars['epochsize'] = int(
            epochsize / pars['n_trajs'] / pars['timestep'])# / pars['stride'])
    return epochsize



# TODO Upgrade to a datasets class and use this to shape the __init__
def create_parameterset(timestep, topfile, directory, filename, n_trajs=1):
    parameterdict = {
                     "timestep"  : timestep,
                     "topfile"   : topfile,
                     "directory" : directory,
                     "filename"  : filename,
                     "trajfiles" : list(),
                     "n_trajs"   : n_trajs,
                     "stride"    : None,
                     "epochsize" : None,
                     "totaldata" : None,
                     "n_clusters": list(),
                     "datashape" : list(),
                    }


