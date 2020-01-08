
#test_timeseries = [0,0,1,1,1,2,3,4,5,5,5,5,5,5,4,3,2,1,0,1,2,3,2,1,0,0,0,1,2,3,4,5,5,5]
#interfaces = [1,4]
#count_changes(test_timeseries, interfaces)

interfaces = [1,4,6,8]
test_timeseries = [0,0,1,1,1,2,3,4,5,5,5,5,6,6,6,7,8,9,9,9,9,9,9,8,7,8,9,9,9,9,8,7,7,6,6,6,5,5,5,4,3,2,1,2,2,3,4,5,5,5,5,5,4,3,2,1,0,1,2,3,2,1,0,0,0,1,2,3,4,5,5,5]


def count_changes(timeseries, interfaces, min_residence=1, **cc_kwargs):
    #if len(interfaces) == 2:
    if False:
        result = count_changes_2state(timeseries, interfaces, min_residence)
    else:
        result = count_changes_Nstate(timeseries, interfaces, min_residence, **cc_kwargs)

    return result


def dtraj_from_changepoints(transitions, full_length):
    import numpy as np
    dtraj       = np.zeros(0, dtype=np.int32)
    for transs in transitions:
        t_fwd,t_bwd = transs
        assert len(t_fwd)==len(t_bwd)==1
        assert len(t_fwd[0])-len(t_bwd[0])<=1
    
        if t_fwd[0][0] < t_bwd[0][0]:
            state0 = 0
        else:
            state0 = 1
    
        if len(t_fwd) > len(t_bwd):
            t_bwd.append(None)
        
        elif len(t_bwd) > len(t_fwd):
            t_fwd.append(None)
    
        if t_fwd[0][-1] > t_bwd[0][-1]:
            stateN = 1
        else:
            stateN = 0
    
        for moves in zip(t_fwd[0],t_bwd[0]):
            if state0 == 0:
                newN = np.zeros(moves[0]-dtraj.shape[0], dtype=np.int32)
                newNN = np.ones(moves[1]-moves[0], dtype=np.int32)
            
            elif state0 == 1:
                newN = np.ones(moves[1]-dtraj.shape[0], dtype=np.int32)
                newNN = np.zeros(moves[0]-moves[1], dtype=np.int32)
            
            dtraj = np.concatenate((
                dtraj,
                newN,
                newNN,
            ))
        
        else:
            try:
                if stateN == 0:
                    newN = np.zeros(full_length-dtraj.shape[0], dtype=np.int32)
                elif stateN == 1:
                    newN = np.ones(full_length-dtraj.shape[0], dtype=np.int32)
        
            except Exception as e:
                print(dtraj.shape, full_length)
                raise e
                
            dtraj = np.concatenate((
                dtraj,
                newN,
            ))
            
    return dtraj


def count_changes_Nstate(timeseries, interfaces, min_residence=1, use_initial_entry=True):
    '''
    region    state or interface or transition placement
    1         dimensions in state space
    N         states
    2(N-1)    interfaces
    N-1       forward changes
    N-1       backward changes

    concept            interface            assigned
                       parity               region idx

       [ state 1 ]                          0
    --- interface ---  o     ^  arrow in    0
       < xregion >           |  backward    1
    --- interface ---  e     |  direction   1
       [ state 2 ]                          2
    --- interface ---  o                    2
       < xregion >                          3
    --- interface ---  e                    3
       [ state 3 ]                          4
    --- interface ---  o                    4
       < xregion >                          5
    --- interface ---  e                    5
       [ state 3 ]                          6
    --- interface ---  o                    6
       < xregion >                          7
    --- interface ---  e                    7

          ...

       < xregion >
    --- interface ---                       2N-1
       [ state N ]                          2N

    Returns
    -------
    (forward_changes, backward_changes)

    2-tuple with sub N-tuples containing counts of
    forward and backward changes into states {N}
    '''

    assert len(interfaces) % 2 == 0 # always even since N states is counting num
                                    #  & inside a 1D dimension (expecting RMSD)
                                    # --------------------------------------------
                                    # 2 interfaces per possible transition
                                    # 1 transition per state count increment
                                    # 2 states, 1 transition minimum

    assert min_residence >= 1  # residence time
                               # before transition "into" state

    n_interfaces       = len(interfaces)
    n_states           = int(n_interfaces / 2) + 1
    #------------------------------------------------#
    # INIT list of lists to store transition indices #
    #------------------------------------------------#
    forward_changes    = [list() for _ in range(n_states - 1)]
    backward_changes   = [list() for _ in range(n_states - 1)]

    def assign_current_index(timepoint):
        '''Get first interface index higher than timepoint value
        '''
        # assume state 1: below interface[0]
        i = 0
        try:
            while timepoint > interfaces[i]:

                # if greater: check next
                i += 1

            # after no satisfy while condition return last i
            else:
                return i

        # timepoint is beyond last interface
        except IndexError:
            return i

    #--------------=====-----------------#
    # TODO allow arbitrary start state   #
    # Initial counting state cannot be the transition region
    #   - loop breaks with first actual state
    ts             = iter(enumerate(timeseries))
    current_region = assign_current_index(next(ts)[1])

    # while in transition region
    #while current_region in xregion:
    while current_region % 2 == 1:
        current_region = assign_current_index(next(ts)[1])
    else:
        last_region = current_region
        last_state  = current_region
    # END TODO for arbitrary start state #
    #--------------=====-----------------#

    #----------=====-----------------#
    # Main loop to count transitions #
    #----------=====-----------------#
    residence_time = -1
    for stepnum, tp in ts:
        current_region = assign_current_index(tp)

        # SAME idea
        #if current_region in states:
        if current_region % 2 == 0:

            # If in new region, residence tracking
            if current_region != last_state:

                # First visit, count starts
                if current_region != last_region:
                    last_region = current_region
                    residence_time += 1

                else:
                    residence_time += 1

                    #---------=====--------------#
                    # OFFICIALLY entered a state #
                    #=========-----==============#
                    if residence_time >= min_residence:
                        # Forwards if new index is higher
                        forwards   = True if last_state < last_region else False
                        last_state = last_region
                        
                        if use_initial_entry:
                            registered_transition_step = stepnum - residence_time

                        else:
                            registered_transition_step = stepnum

                        residence_time = -1

                        if forwards:
                            forward_changes[current_region//2-1].append(registered_transition_step)
                        else:
                            backward_changes[current_region//2].append(registered_transition_step)

            # Residing in state
            else:
                pass

        # In transition region
        else:
            # From a quick exit below `min_residence`
            #   - reset the residence counter
            # TODO FIXME probably don't need reset `last_region`
            if last_region != last_state:
                residence_time = -1
                last_region = last_state

            # From well-established `last_state`
            else:
                pass

    return (forward_changes, backward_changes)



def count_changes_2state(timeseries, interfaces, min_residence=1):
    '''Count changes between 2 primary states separated by a transition region
    Each state is defined by an exit/entrance interface, when the timeseries
    successfully exits a state and enters the other, a change is counted.
    
       [ state 1 ]

    --- interface ---

       < xregion >

    --- interface ---

       [ state 2 ]

    Returns
    -------
    (forward_changes, backward_changes)

    2-tuple with counts of forward and backward changes
    '''

    forward_changes    = list() # from state1 to state2
    backward_changes   = list() # from state2 to state1
    
    assert len(interfaces) == 2 # to use as described
    assert min_residence >= 1 # how long to remain before "in" state

    state1 = "state1"
    state2 = "state2"
    states = [state1, state2]
    xregion = "xregion"

    def check_current(timepoint):
        if timepoint < interfaces[0]:
            return state1
        elif timepoint > interfaces[1]:
            return state2
        else:
            return xregion

    # Initial state cannot be the transition region
    #   - loop breaks with first actual state
    ts = iter(enumerate(timeseries))
    current_region = check_current(next(ts)[1])
    while current_region == xregion:
        current_region = check_current(next(ts)[1])

    else:
        # `_last_state` tracks residence w/out assigning state
        _last_state = current_region
        # `last_state` is assignment of last state
        last_state = current_region

    # Main loop to count transitions
    residence_time = -1
    for stepnum, tp in ts:
        current_region = check_current(tp)

        if current_region in states:

            # If in new region, start residence tracker
            if current_region != last_state:

                if current_region != _last_state:
                    _last_state = current_region
                    residence_time += 1

                # SAME as elif here
                #elif current_region == _last_state:
                else:
                    residence_time += 1

                    if residence_time >= min_residence:
                        last_state = _last_state
                        residence_time = -1

                        if current_region == state2:
                            forward_changes.append(stepnum)

                        else:
                            backward_changes.append(stepnum)

            # Residing in state
            else:
                pass

        # In transition region
        else:
            # From a quick exit below `min_residence`
            #   - reset the residence counter
            if _last_state != last_state:
                residence_time = -1
                _last_state = last_state

            # From well-established `last_state`
            else:
                pass

    return (forward_changes, backward_changes)
