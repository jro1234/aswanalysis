
import inspect

## No successful use yet
#from numba import jit

from count_changes import count_changes, dtraj_from_changepoints
import pandas
import numpy as np

import copy

import pyemma

def refined_regspace(data, mindist=1, min_counts=10, **kwargs):
    dmin = kwargs.pop("dmin")
    dmin *= mindist
    regspace = pyemma.coordinates.cluster_regspace(data, dmin=dmin, **kwargs)
    
    #return regspace
    if isinstance(min_counts, int):
        min_counts = [min_counts]
    
    for min_count in min_counts:
        keep_mask = np.concatenate(np.argwhere(
            np.unique(regspace.dtrajs, return_counts=True)[1] > min_count))
        new_centers = regspace.clustercenters[keep_mask]
        #print("Removing %d of %d centers" % (
        #    regspace.clustercenters.shape[0] - new_centers.shape[0], regspace.clustercenters.shape[0]
        #))
        regspace = pyemma.coordinates.assign_to_centers(data, new_centers, return_dtrajs=False)
    
    return regspace


nworkloads     = lambda nt,ttj: int(len(ttj)/nt)
workloadbounds = lambda nt,ttj: [nt*i for i in range(1+nworkloads(nt,ttj))]

def partition_workflow(w, ttrajs=None):
    if ttrajs is None:
        if "ttrajs" in w:
            ttrajs = w["ttrajs"]
        else:
            return list()
        
    bounds = workloadbounds(w["n_trajs"], ttrajs)
    return zip(bounds[:-1],bounds[1:])


def in_notebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


def copy_include(d, include=None):
    if include:
        return copy.deepcopy({k:v for k,v in d.items() if k in include})
    else:
        return copy.deepcopy({k:v for k,v in d.items()})

    
def get_matching_input(setups, key=None, val=None):
    if isinstance(setups, list):
        for setup in setups:
            if setup["kwargs"][key] == val:
                return setup["result"]
    # For consistent use by scan_single_parameter
    # when the input isn't from a list
    else:
        return setups


def iter_models(all_models):
    for features in all_models:
        for nm, dataset in all_models[features].items():
            yield features, nm, dataset

        
def scan_single_parameter(
    calculation, models, parameter, fixed_pars=dict(),
    wkf_pars=dict(), wkf_include=lambda w: w,
    inp_source=None, inp_key=None, inp_val=None
):

    for feat,nm,dataset in iter_models(models):
        if not wkf_include(nm):
            continue

        dataset[calculation] = news = list()
        inp = get_matching_input(dataset[inp_source], inp_key, inp_val)

        pnm, pvals = list(parameter.items())[0]
        for p in pvals:
            kwargs = {pnm: p}
            kwargs.update(fixed_pars)
            kwargs.update(wkf_pars(nm))
            
            newd = dict(result=False, input=inp, kwargs=kwargs)

            news.append(newd)


#class ColorCycler(object):
#    # TODO color mod application generalization
#    _apply = lambda f,c: [sum(r) for r in zip(c, (-0.05, -0.05, -0.05, 0))]
#    _darken  = lambda c: [sum(r) for r in zip(c, (-0.05, -0.05, -0.05, 0))]
#    _redder  = lambda f,c: [sum(r) for r in zip(c, (f*0.05, 0, 0, 0))]
#    _lesred  = lambda f,c: [sum(r) for r in zip(c, (f*-0.05, 0, 0, 0))]
#    _greener = lambda f,c: [sum(r) for r in zip(c, (0, f*0.05, 0, 0))]
#    _colors = [
#        _darken(_lesred(3, colorcycle[1])),
#        _darken(_redder(3, colorcycle[3])),
#        _darken(_greener(3, colorcycle[4])),
#        _darken(colorcycle[5]),
#        _darken(_darken(_darken(colorcycle[6]))),
#    ]
#    def __len__(self):
#        return len(self._colors)

#    def __init__(self):
#        super(ColorCycler, self).__init__()
#        self._counter = -1

#    def __next__(self):
#        self._counter += 1
#        return self._colors[self._counter % len(self)]



def mix_parameters(pars_dicts):
    import itertools

    pars_lists = dict()
    #for par, values_lists in pars_dicts.items():
    #    if not isinstance(values_lists[0], list):
    #        pars_dicts[par] = [[v] for v in values_lists]

    pars_names, _pars_lists = list(zip(*pars_dicts.items()))
    pars_lists = {k : list() for k in pars_names}

    for product in itertools.product(*_pars_lists):
        [pars_lists[k].append(v) for k,v in zip(pars_names, product)]
        
    return pars_lists


#@jit
def count_transitions_runningtotal(
    transitions, increments, increment_length):
    return [
        len(list(filter(
            lambda x: x<(i*increment_length), transitions
        ))) for i in range(1,1+increments)
    ]


#@jit doesn't work with highly embedded lambda whatever
def observed_transition_rate(
    transitions, increments, increment_length):
    return [
        float(len(list(filter(
            lambda x: x<(i*increment_length), transitions
        )))) / float(i*increment_length)
        for i in range(1, increments+1)
    ]


# FIXME this only does forwards and backwards
#       to/from state 0
def calculate_observed_rates(
    transition_data, incr_length, n_incr, step_per_ns):
    '''Return values in MHz
    '''

    observed_rates = dict()

    for key, trannies in transition_data.items():

        trannies = trannies[0]
        forwards, backwards = trannies[0][0], trannies[1][0]

        #if key == "OpenMM_interfaces-[0.13, 0.3]_min_residence-10":
        #    print(key)
        #    print(key)
        #    print(forwards)
        #    print(backwards)

        #if isinstance(forwards[0], list):
        #    forwards = forwards

        datakey = key.split("_")[0]
        _n_incr = n_incr[datakey]

        observed_rates.update({key : {
            "forwards" :
                pandas.Series([
                    r * step_per_ns[datakey] * 1e3# ns/s / MHz/Hz
                    for r in observed_transition_rate(
                        forwards, _n_incr, incr_length)
                ]),
            "backwards" :
                pandas.Series([
                    r * step_per_ns[datakey] * 1e3# ns/s / MHz/Hz
                    for r in observed_transition_rate(
                        backwards, _n_incr, incr_length)
                ]),
        }})
        #if key == "OpenMM_interfaces-[0.13, 0.3]_min_residence-10":
        #    print(["%f, "%f for f in observed_rates[key]["forwards"]])
        #    print(["%f, "%f for f in observed_rates[key]["backwards"]])
        #print(len(count_transitions_runningtotal(trannies, nincr, incr_length, incr_length)))
    return observed_rates


def _orig_only_forwards_from_state0_calculate_observed_rates(
#def calculate_observed_rates(
    transition_data, incr_length, n_incr, step_per_ns):

    observed_rates = dict()

    for key, trannies in transition_data.items():
        trannies = trannies[0][0]
        if trannies and isinstance(trannies[0], list):
            trannies = trannies[0]

        datakey = key.split("_")[0]
        _n_incr = n_incr[datakey]

        observed_rates.update({key : 
            pandas.Series([
                r / step_per_ns[datakey] * 1e9
                for r in observed_transition_rate(
                    trannies, _n_incr, incr_length)
            ])}
        )
        #print(len(count_transitions_runningtotal(trannies, nincr, incr_length, incr_length)))
    return observed_rates

