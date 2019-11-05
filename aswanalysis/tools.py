
import subprocess, shlex
import sys
import os
from glob import glob

import mdtraj
import pyemma
pyemma.config.show_progress_bars=False
import pyemma.coordinates as coor
import pyemma.msm as msm
import numpy as np

from pprint import pformat

# TODO add adaptivemd imports
# from adaptivemd.analysis.pyemma import apply_feat_part


# TODO use logger 
def get_logger(logname, logfolder=''):
    import logging

    if logname in {"builtins","__main__"}:
        import datetime
        logname = "session.{}".format(
             datetime.datetime.now().strftime('%m-%d-%y.%H:%M:%S'))

    loglevel  = logging.INFO
    logfile   = os.path.join(logfolder, logname+'.log')
    formatter = logging.Formatter(
        '%(asctime)s :: %(name)s :: %(levelname)s' + \
        ' :: %(lineno)d |||   %(message)s')
    logging   . basicConfig(level=loglevel)#, format=formatter)
    logger    = logging.getLogger(logname)
    ch        = logging.StreamHandler()#sys.stdout)
    ch        . setLevel(loglevel)
    ch        . setFormatter(formatter)
    fh        = logging.FileHandler(logfile)
    fh        . setLevel(loglevel)
    fh        . setFormatter(formatter)
    logger    . addHandler(ch)
    logger    . addHandler(fh)
    logger    . propagate = False
    return logger


listit = lambda x: x if isinstance(x, (list,tuple,set)) else [x]
format_spacedlist = lambda l: ' '.join(['{}'.format(ll) for ll in l])
#min_full_datas = lambda dts: min([np.sum(np.sum([dss for dss in ds['datashape']])/ds['stride'])
#                                  for nm,ds in dts.items()])
min_full_datas = lambda dts: min([int(np.sum([[(dsss-1)/ds['stride'] for dsss in dss]
                                  for dss in ds['datashape']])) for nm,ds in dts.items()])
onlydatas = lambda datasets: {k:v for k,v in datasets.items() if k is not 'analysis'}


def squish_tica_inputfiles(datasets, feat):
    ticainputs = list()
    # The 'dim_reduction' dict must only have parameters for
    # the function call, store 'input_order' separately
    input_order = list()
    #datasets['analysis']['dim_reduction']['input_order'] = list()
    for nm,pars in datasets.items():
        ticainputs.extend(pars.get('tica_inputs', []))
        input_order.append(nm)
        pars['feat']   = feat
    return ticainputs, input_order


# FIXME belongs as first-class function in adaptivemd
def apply_feat_part(featurizer, parts):
    if isinstance(parts, dict):

        items = list(parts.items())
        if len(items) == 1:
            func, attributes = items[0]
            kwargs = dict()

        elif len(items) == 2:
            if items[0][0] == 'kwargs':
                func, attributes = items[1]
                key, kwargs = items[0]

            elif items[1][0] == 'kwargs':
                func, attributes = items[0]
                key, kwargs = items[1]

            for k,v in kwargs.items():
                if isinstance(v, dict):

                    _func, _attr = list(v.items())[0]
                    _f = getattr(featurizer, _func)
                    if _attr is None:
                        idc = _f()

                    elif isinstance(_attr, (list, tuple)):
                        idc = _f(*apply_feat_part(featurizer,
                                 _attr))

                    kwargs[k] = idc

        assert isinstance(kwargs, dict)
        f = getattr(featurizer, func)

        if attributes is None:
            return f(**kwargs)

        elif isinstance(attributes, (list, tuple)):
            return f(*apply_feat_part(featurizer, attributes),
                     **kwargs)
        else:
            return f(apply_feat_part(featurizer, attributes),
                     **kwargs)

    elif isinstance(parts, (list, tuple)):
        return [apply_feat_part(featurizer, q)
                for q in parts]
    else:
        return parts

# TODO topology decorator
#@md_topology
def prepare_tica_inputs(datasets, topfile, features=None, selection=None, chunksize=10000, singletraj=False):
    #print("topfile: ", topfile)
    #print("selection", selection)
    if isinstance(topfile, mdtraj.Topology):
        topology = topfile
    elif os.path.exists(topfile):
        topology = mdtraj.load(topfile).topology

    assert isinstance(topology, mdtraj.Topology)

    if selection:
        topology = topology.subset(topology.select(selection_string=selection))

    #if isinstance(features, featurizer):
    #    feat = features
    feat = coor.featurizer(topology)
    if not features: # then use inverse Ca distances
        # PyEMMA equivalent: `feat.add_inverse_distances(feat.select_backbone())`
        features = {'add_inverse_distances': { 'select_Ca': None }}

    apply_feat_part(feat, features)
    ticainputs, input_order = squish_tica_inputfiles(datasets, feat)
    if singletraj:
        ticainputs = [ticainputs]
    #print("remove this & below!")
    #print(ticainputs)
    #datasets['analysis']['dim_reduction']['input_order'].append(nm)
    tica_inp = coor.source(ticainputs, feat, chunksize=chunksize)
    return tica_inp, input_order


def default_msm_analyze(d, a):
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
    # TODO rewrite `prepare_data` argument structure
    #      to accommodate call signature changes
    #      to downstream functions
    # NOW starting to process data
    create_trajlist     (dataset)
    assign_stride       (dataset)# add namekey
    determine_epochsize (dataset)
    determine_datashape (dataset, analysis)


def bind_ttraj_results(datasets, input_order, ttrajs):
    # This tries to reference dtrajs correctly from dataset dicts
    n_trajs_bound = 0
    for i,nm in enumerate(input_order):
        pars = datasets[nm]
        if pars.get('trajfiles', []):
            these          = len(pars['trajfiles'])
            pars['ttrajs'] = ttrajs[n_trajs_bound:n_trajs_bound+these]
            n_trajs_bound += these


def bind_dtraj_results(datasets, input_order, dtrajs):
    # This tries to reference dtrajs correctly from dataset dicts
    n_trajs_bound = 0
    for i,nm in enumerate(input_order):
        pars = datasets[nm]
        if pars.get('trajfiles', []):
            these          = len(pars['trajfiles'])
            pars['dtrajs'] = dtrajs[n_trajs_bound:n_trajs_bound+these]
            n_trajs_bound += these


def md_traj(trajfiles, topfile, stride=1):
    trajfiles = listit(trajfiles)

    assert os.path.exists(topfile)
    assert all([os.path.exists(tf) for tf in trajfiles])

    traj = mdtraj.load(trajfiles, top=topfile)
    #traj = mdtraj.load(trajfiles[0], top=topfile)
    #for tf in trajfiles[1:]:
    #    traj += mdtraj.load(tf, top=topfile)
    return traj[::stride]


# FIXME TODO use loading decorator to replace md_traj
#            and annotate traj loading function
def get_traj(func, trajfiles, topfile=None, stride=1):
    # TODO this and rmsd--> take existing traj as argument
    if isinstance(trajfiles,list) and topology:
        traj = md_traj(trajfiles, topfile, stride)
    elif isinstance(trajfiles,mdtraj.Trajectory):
        pass
        
    func()



def prep_rmsd_calc(traj, ref=None, selection='', precentered=False, parallel=True):
    if not ref:
        ref = traj
    elif isinstance(ref, str):
        ref = mdtraj.load(ref)
    kwargs_rmsd = dict()
    kwargs_rmsd['precentered'] = precentered
    kwargs_rmsd['parallel']    = parallel
    if selection:
        kwargs_rmsd['ref_atom_indices'] = ref.topology.select(selection)
        kwargs_rmsd['atom_indices']     = traj.topology.select(selection)
    if precentered:
        traj.center_coordinates()
    return kwargs_rmsd


def compute_dssp(trajfiles, topfile=None, stride=1, simplified=True):
    # TODO this and rmsd--> take existing traj as argument
    traj = md_traj(trajfiles, topfile, stride)
    dssp = mdtraj.compute_dssp(traj, simplified)
    return dssp


def dssp_per_residue(dssp_traj):
    dssp_residues = [{"H":0, "E":0, "C":0} for _ in dssp_traj.shape[0]]
    for dssp_frame in dssp_traj:
        for i, dssp_res in enumerate(dssp_frame):
            dssp_residues[i][dssp_res] += 1

    return dsp_residues


def rmsd_rolling(trajfiles, topfile, lag=1, **kwargs):
    # see TODO's for `rmsd_traj` and apply as applicable
    # TODO --> is there a faster way for this to work?
    #          ref seems always a single frame in MDTraj.rmsd
    #          maybe just multiprocess calls of next_rmsd
    #          and add to a shared dict, then return as list
    #           shared_dict :: key<idx> : value<rmsd>
    '''
    RMSD is calculated on a rolling basis from `lag` number of frames
    back. This function might be useful for checking the 'quality' or
    'sanity' of your simulation data, ie that there isn't a huge jump
    in the structure between saved frames.

    Notes: 
    - `precentered` to this function is an operation to be done, which
    is passed as `precentered` to the MDTraj function.
    - centering and RMSD calculation, not loading, timed if `timeit=True`
    - see MDTraj notes on centering
      http://mdtraj.org/development/api/generated/mdtraj.rmsd.html?highlight=precentered#
    '''
    traj        = md_traj(trajfiles, topfile)
    kwargs_rmsd = prep_rmsd_calc(traj, **kwargs)
    # TODO wretchedly slow from changing reference and re-emitting mdtraj call
    next_rmsd   = lambda i: mdtraj.rmsd(traj[i], traj[i+lag], **kwargs_rmsd)
    rmsds       = np.concatenate([next_rmsd(i) for i in range(len(traj) - lag)])
    return rmsds
    

#def rmsd_traj(trajfile, topfile, reference=None, selection='', precenter=False, timeit=False, parallel=True):
    # TODO  option of separate reference (PDB) file incase of compare
    #       slightly different topologies, then topfile is only used
    #       as the trajfile topology source.
def rmsd_traj(trajfiles, topfile=None, **kwargs):
    # TODO  look at timing of with/without precentering and eliminate
    #       option if it doesn't decrease the calculation time
    '''
    RMSD of frames in trajfile to reference structure in topfile.
    They (currently) need to have identical topology as topfile is
    both the reference structure and topology source. The logic here
    requiresd that `topfile` is a PDB filepath, not trajectory.

    Notes: 
    - `precenter` to this function is an operation to be done, which
    will result in `precentered` to the MDTraj function as True.
    - centering and RMSD calculation, not loading, timed if `timeit=True`
    - see MDTraj notes on centering
      http://mdtraj.org/development/api/generated/mdtraj.rmsd.html?highlight=precentered#
    '''
    ref         = mdtraj.load(topfile)
    traj        = md_traj(trajfiles, topfile)
    kwargs_rmsd = prep_rmsd_calc(traj, **kwargs)
    rmsds       = mdtraj.rmsd(traj, ref, **kwargs_rmsd)
    return rmsds


# TODO isolate property sifting function and give rmsd as property
#      get function to `microstates_rmsd` as a wrapper
def microstates_rmsd(dtrajs, trajfiles, topfile, return_traj=False):
    '''
    This function collects the observed RMSD distribution for each microstate
    from the given trajectory dataset and microstate:trajectory mapping. 
    '''
    rmsds = [rmsd_traj(trajfile, topfile) for trajfile in listit(trajfiles)]
    assert len(dtrajs) == len(rmsds)
    micrormsds = [[rmsds[i][np.where(dtrajs[i]==mdx)] for i in range(len(rmsds))] for mdx in np.unique(np.concatenate(dtrajs))]
    if return_traj:
        return micrormsds, rmsds
    else:
        return micrormsds


def save_microstate_structures(filename_tpt, mm, trajfiles, topfile=None, n_frames=10):
    sampled_frames = mm.sample_by_state(n_frames)
    outfiles = [filename_tpt.format(state) for state in mm.active_set]
    assert len(sampled_frames) == len(outfiles)
    if not isinstance(trajfiles, list):
        trajfiles = [trajfiles]
    if all(isinstance(trajfile, mdtraj.Trajectory) for trajfile in trajfiles):
        trajs = trajfiles
    else:
        if not topfile:
            raise Exception("Need to give topfile to load new trajectories")
        trajs = [mdtraj.load(trajfile, top=topfile) for trajfile in trajfiles]
    for otf, sampling in zip(outfiles, sampled_frames):
        for i in range(1+max(sampling[:,0])):
            fp = otf.split('.')
            fp[-2] += '-traj_{}'.format(i)
            trajs[i][sampling[np.where(sampling[:,0]==i)[0]][:,1]].save('.'.join(fp))
    return sampled_frames


def compare_sets(dataset, n_frames):
    if isinstance(n_frames,int) and n_frames > 0:
        nfc = n_frames
    else:
        nfc = n_frames_compare(dataset)


class MongoInstance(object): # TODO replace path insert formatting with os.path.join
    def __init__(self, dbpath, dbport=27017):
        super(MongoInstance, self).__init__()
        assert isinstance(dbpath, str) # TODO check is path-like, parent dir exists
        #assert isinstance(dbport, validport) # TODO isint and betwen X and Y
        self.dbpath = dbpath
        self._dblog_file = None
        self.dbport = dbport
    @property
    def pid(self):
        return self.mongo_proc.pid
    @property
    def mongo_proc(self):
        return self._mongo_proc
    @property
    def dblogfile(self):
        return self._dblog_file
    def open_mongodb(self, remove_socket=False, remove_locks=False):
        if remove_socket:
            self._remove_socket_file()
        if remove_locks:
            self._remove_lock_files()
        self._write_config_file()
        self._dblog_file = open("{0}/db.log".format(self.dbpath), 'w')
        mongo_launcher = "mongod --dbpath {0}/db --config {0}/db.cfg"
        self._mongo_proc = subprocess.Popen(shlex.split(mongo_launcher.format(self.dbpath)), stdout=self._dblog_file, stderr=self._dblog_file)
    def stop_mongodb(self):
        self._mongo_proc.kill()
        self._mongo_proc.wait()
        self._dblog_file.close()
        self._remove_socket_file()
        self._remove_lock_files()
    def _remove_lock_files(self):
        mongo_lock_file = '{0}/db/mongod.lock'.format(self.dbpath)
        wt_lock_file = '{0}/db/WiredTiger.lock'.format(self.dbpath)
        if os.path.exists(mongo_lock_file):
            os.remove(mongo_lock_file)
        if os.path.exists(wt_lock_file):
            os.remove(wt_lock_file)
    def _remove_socket_file(self):
        socket_file = '{0}/socket/mongodb-{1}.sock'.format(self.dbpath, self.dbport)
        if os.path.exists(socket_file):
            os.remove(socket_file)
    def _write_config_file(self):
        config_string = "net:\n   unixDomainSocket:\n      pathPrefix: {0}/socket\n   bindIp: 0.0.0.0\n   port:   {1}\n".format(self.dbpath, self.dbport)
        config_file = "{0}/db.cfg".format(self.dbpath)
        with open(config_file,"w") as cfg:
            cfg.write(config_string)


# TODO remove or specify purpose of epoch_coeff
#def calc_eq_explored(dataset, mm, epoch_coeff=1, n_frames=0):
def calc_eq_explored(dataset, mm, n_frames=0):
    # We're going to get the total equilibrium population explored
    # for each dataset with consistent epochsizes across datasets
    # ie more parallel <--> smaller chunks per traj
    #  - if a dataset gives initial frame,
    #    this code assumes a single traj from this frame
    #    and that there is more data here than in other datasets
    if isinstance(n_frames,int) and n_frames > 0:
        nfc = n_frames
    else:
        nfc = n_frames_compare(dataset)
    #print("Comparing {} frames at the master stride interval".format(nfc))
    for nm,pars in dataset.items():
        initial = pars.get('initial', 0)
        if initial:
            traj_initial, frame_initial = initial # maybe not zero
            # composite_trajs coordinates the indices used
            # to build the dtja/_array vars below with actual data
            composite_trajs = [traj_initial] # using individual segment
            n_trajs = 1
            frame_final = int(frame_initial + nfc)
            epochsize = int(pars['epochsize']/pars['stride']*pars['n_trajs'])
            #epochsize = pars['epochsize']/pars['stride']*epoch_coeff
        else:
            frame_initial = initial # == 0
            composite_trajs = range(pars['n_trajs'])
            n_trajs = len(pars['dtrajs']) # not same as pars['n_trajs']
            frame_final = -1
            epochsize = int(pars['epochsize']/pars['stride'])
        print(
  "Epoch size of {0} frames for dataset {1}".format(epochsize, nm))
        print(
  "Data from \"{0}\" workflow will start at frame {1}".format(nm, frame_initial))
        pars['eq_explored'] = list()
        pars['clusters']    = list()
        pars['clust_explored']  = list()
        # the dtj_array initialization logic forces single traj processing
        # if an initial frame was given
        dtj_array = [[]] if initial else [[] for _ in composite_trajs]
        # this line stacks the trajs from each workload in iteration order
        #     pars['dtrajs']         dtj_array
        #      [ traj1,             [ traj1traj3,
        #        traj2,      -->      traj2traj4,
        #        traj3,             ]
        #        traj4,             2 rounds of 2 replicates
        #      ]
        # or selects correct single trajectory if using data starting from
        # a single frame
        [dtj_array[i].extend(np.concatenate(
         pars['dtrajs'][tj_i:tj_i+n_trajs][
         ::len(composite_trajs)]))
         for i,tj_i in enumerate(composite_trajs)]
        # this line instantiates numpy arrray, trims
        # only if using a single trajectory segment
        # TODO     more general method (see below fixme)
        # FIXME if dtrajs are inhomogenous then need
        #       to avoid extra clusters
        #       - maybe... end with sequence of
        #         last cluster to make array
        dtja = np.array([np.array(_arr[frame_initial:frame_final])
                         for _arr in dtj_array])#, axis=0)
        del dtj_array
        print("SHAPE: ", dtja.shape)
        epoch  = 0 # incremented to give first frame of each epoch
        pars['clusters'].append(
              np.unique(dtja[:,0]))
        pars['clust_explored'].append(len(np.unique(
              np.concatenate(pars["clusters"]))))
        pars['eq_explored'].append(np.sum(
              mm.stationary_distribution[np.unique(
              np.concatenate(pars["clusters"]))]))
        while epoch*len(composite_trajs) < nfc:
            try:
                #print(epoch, type(epoch))
                #print(epochsize, type(epochsize))
                pars["clusters"].append(np.unique(
                        # this only looks at unique
                        # cluster IDs in the epoch
                        dtja[:,epoch:epoch+epochsize]
                ))
                # this collects unique cluster IDs
                # found in all epochs
                exploredclusters = np.unique(
                      np.concatenate(pars["clusters"]))
                pars["clust_explored"].append(
                      len(exploredclusters))
                pars["eq_explored"].append(np.sum(
                      mm.stationary_distribution[exploredclusters]))
            except IndexError: # FIXME pretty sure this isn't useful as is
                print("Can't get uneven final chunk")
                print(epoch * pars['n_trajs'], "<", nfc)
            epoch += epochsize


def n_frames_compare(datasets):
    mfd = min_full_datas(onlydatas(datasets))
    return mfd


def _check_square_same(P, Q):
    assert len(P.shape) == 2
    assert P.shape[0]   == P.shape[1]
    assert P.shape == Q.shape
    return True


def calc_Tmatrix_ratio(P, Q):
    assert _check_square_same(P,Q)
    I,J = [range(_) for _ in P.shape]
    D   = list()
    for i in I:
        D.append(list())
        for j in J:
            pij = P[i,j]
            qij = Q[i,j]
            if pij == 0 or qij == 0:
                if pij == qij:
                    dij = 1e-10 # 0 for log colormap
                elif pij > qij:
                    dij = 100
                else:
                    dij = 10000
            else:
                dij = np.log10( pij / qij )
            D[-1].append(dij)
    return np.array(D)


def calc_Tmatrix_difference(P, Q):
    assert _check_square_same(P,Q)
    I,J = [range(_) for _ in P.shape]
    return np.array([[ P[i,j] - Q[i,j] for j in J] for i in I])


def expand_matrix(M, existing, size):
    '''
    M is treated as a square matrix
    '''
    
    return M


# FIXME no prior assumed before C-->T calculation
#   - probably need to take C matrices and calc T
#     (ie P,Q) and then the relative entropy
def calc_rel_entropy(u, P, Q, Q_states=None, reduce=False):
    '''
    Q can be expanded to match dims of P if not all states
    are represented in Q, or else P contracted to compare
    only the entries existing in both matrices

    u :: array of reference equilibrium distribution
    P :: reference distribution transition matrix
    Q :: test distribution transition matrix
    Q_states :: array of the state indices in Q
         if Q is smaller than P
    '''
    if Q_states:
        if reduce:
            # rows,cols only from states in Q
            P = P[Q_states,Q_states] # FIXME if assert fails
        else:
            Q = expand_matrix(Q, Q_states, size)
    assert _check_square_same(P,Q)
    D   = [0 for i in range(P.shape)]
    I,J = [range(_) for _ in P.shape]
    for i in I:
        for j in J:
            pij = P[i,j]
            qij = Q[i,j]
            dij = u[i] * pij * np.log( pij / qij )
            D[i] += dij
    return D


def assign_trajs_to_clusters(datasets, MORE):
    inp_files = min_full_datas(datasets)
    for nm, pars in datasets.items():
        #if nm is not master_dataset:
        #if nm is not 'single trajectory':
        print('Projecting {} into clustercenters'.format(nm))
        feat = pars['feat']
        inp = coor.source(inp_files[nm], feat, chunksize=1000)
        dtrajs = coor.assign_to_centers(
                inp,
                centers=all_clust.clustercenters,
                stride=pars['stride'],
                metric='minRMSD',
                chunksize=1000
                )
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



def create_trajlist(datasets, topdir='.', sortkey=None, filterkey=None, minimum_filesize=1):
    '''
    `datasets` argument is a `dict` of `dict` dataset definitions who all
    share a common home directory given by `topdir`. Subdirectories are
    determined via the dataset configurations. `sortkey` is passed as the
    `key` argument to the python `sorted` builtin if you need a special way
    to sort your trajectory files by name. `minimum_filesize` can
    be passed (in bytes) to prevent some small files from being included.
    An addition function `filterkey` can be provided to do additional
    screening of the trajectory file list.

    Each dataset sub-`dict` must have must have fields:
      - directory - relative home directory of trajectory files
      - filename  - pattern `glob` matches to pull the trajectory files
    Result is stored in each dataset `dict` as a new / replaced field:
      - trajlist  - `list` for storing paths to trajectory files
    '''
    assert isinstance(minimum_filesize, int)
    assert os.path.exists(topdir)

    if sortkey:
        assert callable(sortkey)

    if filterkey:
        assert callable(filterkey)

    for nm, pars in onlydatas(datasets).items():

        filenames = os.path.join(topdir, pars['directory'], pars['filename'])
        first     = pars.get('first', 0)
        last      = pars.get('last', None)
        trajfiles = glob(filenames)
        if filterkey:
            trajfiles = filter(filterkey, trajfiles)

        pars['trajfiles'] = list(filter(
                lambda fname: os.stat(fname).st_size > minimum_filesize,
                sorted(trajfiles, key=sortkey)[first:last]
        ))


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
        print("USING this search key: ", filenames)
        pars['tica_inputs'] = sorted(glob(filenames))
        print("found these guyts: ", pars['tica_inputs'])


def build_call(calldict, parent, *args):
    # TODO register and offer report of function call
    '''
    This is a dangerous but convenient function
    Build a call signature for function/method of parent
    given by calldict top key. 
    len(calldict) === 1
    Parent can be module or object with method of name

    Arguments
    ---------
    '''
    if len(calldict) != 1:
        return None
    fname = list(calldict.keys())[0]
    f = getattr(parent, fname)
    kwargs = calldict[fname]
    return lambda: f(*args, **kwargs)


# THIS IS SUPER SLOW
# AND UNNECESSARILY DUPLICATES THE STRIDE 1 FILES
# Can't figure out how to access DCD header via mdtraj...
# ... no other way to get ahold of length without full read
#  - a better chunksize for big-node processing would be like 100,000
# TODO   datashape isn't used anymore --> except for total frame counting
#      FIXME   just make feature trajs? -> ?^^^?
def determine_datashape(datasets, analysis=None, topdir='.',
               prepare=True, selectiontag='-stride-', chunksize=10000,
               maxchunks=False, atomselection="all"):
    # FIXME this function is doing too much
    #        - seperate traj reduction, navigating topo
    #          and index files, getting datashape
    '''
    This function provides some way of creating (as `list`)
    a shape for each dataset that is used to navigate equal
    size epochs. Creates a sub-`list` for each workload,
    with each entry corresponding to a trajectory length
    from that workload. The trajs in a workload are considered
    degenerate/replicates, but the workloads are in a predetermined order
    given by their order-of-appearance in the trajlist.

    If `prepare` is `True`, save trajectories in correct
    stride for the analysis--> add to `tica_inputs`

    Fields required in each dataset sub-`dict`:
     - Without correct n_trajs (workload size), datashape will mix workloads
     - A (final) workload of different size will be missed entirely via
       the `j` indexing strategy
    '''
    if analysis: # configuration should be compatible for all datasets
        all_topofile           = analysis['topfile']
        analysis_atomselection = analysis.get('atomselection', 'all')

    else: # FIXME using first dataset's topfile as all topfile
        source = datasets[list(datasets)[0]]
        all_topofile = os.path.join(source['directory'], source['topfile'])

    all_topofile = str(selectiontag+"-selection").join(
            [os.path.join(topdir, all_topofile).split('.pdb')[0], '.pdb'])

    if atomselection != 'all':
        analysis_atomselection = atomselection

    for nm,pars in onlydatas(datasets).items():
        fntag = "{0}-{1}".format(selectiontag, pars['stride'])

        if 'tica_inputs' not in pars:
            pars['tica_inputs'] = list()
        if 'datashape' not in pars:
            pars['datashape'] = list()

        datashape = pars['datashape']
        # TODO replace `j` with workload
        j    = -1
        done = False
        if prepare:
            f_topo_atomselection = os.path.join(
                    topdir,
                    pars['directory'],
                    pars['topfile'].split('.pdb')[0]+fntag+'.atoms.idx')
            #print("\n --- Now using this atom index file ", f_topo_atomselection)
            pars['topfile_atoms'] = all_topofile

        while not done:
            j += 1
            try:
                #print "on j=",j
                #print "\nOn set ", nm
                lengths = list()
                for i in range(pars['n_trajs']):
                    tfn = pars['trajfiles'][ i + j*pars['n_trajs'] ]  # traj filename
                    #print("LOADING TRAJ FROM:")
                    #print("  --> TRAJFILE: ", tfn)
                    #print("  --> TOPFILE : ", os.path.join(pars['directory'],pars['topfile']))
                    tot = 0  # total frames in traj
                    tfs = '' # trajectory file atomselected
                    adx = None
                    if prepare:
                        # # TODO implement offset for proper stitching
                        # #      when mixing chunks and stride
                        # offset = 0
                        # tfs = "{0}--stride-{1}{{0}}.dcd".format(
                        #           tfn.split('dcd')[0], pars['stride'])
                        tfs = "{0}-{1}.dcd".format(
                                tfn.split('.dcd')[0], fntag)

                        if os.path.exists(tfs): # FIXME assuming prepare complete with 1 check
                            #print("Found existing trajectory with atomselection: ", tfs)
                            #pars['tica_inputs'].append(tfs)
                            #continue
                            pass

                    tfi = mdtraj.iterload(
                            tfn,
                            top=os.path.join(pars['directory'],pars['topfile']),
                            chunk=chunksize)

                    #print "\nlooking at this traj ", tfn
                    #print "using this traj loader ", tfi
                    for i,segment in enumerate(tfi):
                        #print segment
                        if i == 0:
                            adx = segment.topology.select(analysis_atomselection)

                        if all_topofile:
                            if not os.path.exists(all_topofile):
                                #print("Going to write this topo",segment.atom_slice(adx).topology)
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
                    #print("added this many frames traj to datashape ", str(tot))
                    #mdtraj.join([tfs.format(k) for k in range(i)])
                    # dependent on `prepare` and file existence
                    if prepare:
                        if prepare and not os.path.exists(f_topo_atomselection):
                            #print("Creating new atom index file: ", f_topo_atomselection)
                            #print("to print these indices: {}".format(adx))
                            with open(f_topo_atomselection, 'w') as f:
                                f.write(format_spacedlist(adx))

                        else:
                            #print("Found existing atoms index file: ", f_topo_atomselection)
                            pass

                        if not os.path.exists(tfs):
                            #print("using this file to store output atomselection ", tfs)
                            subprocess.call(shlex.split(
                                  "mdconvert -o {output_str} -a {selection} -s {stride} {input_str} ".format(
                                  output_str=tfs, stride=pars['stride'], input_str=tfn, selection=f_topo_atomselection)))

                        else: # TODO this all messed up, & we shouldn't get here
                            #print("Found already existing output file: ", tfs)
                            pass

                        pars['tica_inputs'].append(tfs)

                datashape.append(lengths)

            except IndexError:
                #print("Done reading from set ", nm, "\n")
                done = True


def determine_epochsize(datasets):
    # FIXME this function is printing something somehow, where? -> remove
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
def prime_dataset(dataset):
    added_parameters    = {"n_clusters",}
    optional_parameters = {"first","last","initial"}
    preset_parameters   = {"timestep","topfile","directory","filename","n_trajs"}
    required_parameters = {
           "trajfiles" : list(),
           "stride"    : None,
           "epochsize" : None,
           "totaldata" : None,
           "datashape" : list(),
    }
    for rqd,vi in required_parameters.items():
        if not dataset.get("rqd", 0):
            dataset[rqd] = vi
    if not all([preset in dataset for preset in preset_parameters]):
        print("WARNING: Not all parameters requiring preset values are found in dataset definition.")


normalize = lambda vec: vec/np.linalg.norm(vec)

###  TODO get_(cg)msm method for msm objects
###       then use this function to get distribution
def calc_xma_distribution(model, num_macrostates=25, reversible=True):
    import msmtools
    def MinMaxScale(X, min=-1, max=1):
        X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
        X_scaled = X_std * (max - min) + min
        return X_scaled
    data = model
    c         = data['msm']['C']
    counts    = np.array(np.sum(c, axis=1), dtype=int)
    array_ok  = msmtools.estimation.largest_connected_set(c)
    num_macrostates = min(num_macrostates, array_ok.shape[0]/3)
    connected = msmtools.estimation.is_connected(c[array_ok,:][:,array_ok])
    disconnected_microstates = [i for i in range(c.shape[0]) if i not in array_ok]
    p = msmtools.estimation.transition_matrix(c[array_ok,:][:,array_ok], reversible=reversible)
    current_MSM_obj = pyemma.msm.markov_model(p)
    current_MSM_obj.pcca(num_macrostates)
    #macrostate_assignments = { k:v for k,v in enumerate(current_MSM_obj.metastable_sets) }
    macrostate_assignment_of_visited_microstates = current_MSM_obj.metastable_assignments
    corrected_num_macrostates = num_macrostates + len(disconnected_microstates)
    corrected_ma_assignment = np.zeros(c.shape[0])
    corrected_ma_assignment[array_ok] = macrostate_assignment_of_visited_microstates
    for n,i in enumerate(disconnected_microstates):
        corrected_ma_assignment[i]=n+num_macrostates
    macrostate_counts = np.array([np.sum(counts[corrected_ma_assignment == macrostate_label]) for macrostate_label in range(corrected_num_macrostates)])
    ma_counted = macrostate_counts
    ma_counted = ma_counted + (np.sum(ma_counted)/float(ma_counted.shape[0]))**0.5
    ma_pdist = 1.0 / ma_counted
    ma_pdist[macrostate_counts==0] = 0 # Don't select unvisited states
    ma_pdist = normalize(ma_pdist)
    mi_pdist = np.zeros(corrected_ma_assignment.shape[0])
    for ma in range(corrected_num_macrostates):
        mi_pdist[corrected_ma_assignment==ma] = ma_pdist[ma]
    return mi_pdist

