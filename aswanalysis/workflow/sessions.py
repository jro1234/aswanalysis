
__all__ = [
    "read_sessions"
]

# TODO This stuff should all be done via class
#      have giant function instead. Break at
#      the organization points and create
#      attributes that store the bridging data

import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family']     = "Arial"
from matplotlib import pyplot as plt

#=====----------------------------------------=====#
# From outside packge                              #
#=====----------------------------------------=====#
import os
import sys
import yaml

from glob import glob
from pprint import pformat

import numpy as np

#=====----------------------------------------=====#
# From within package                              #
#=====----------------------------------------=====#
from . import timestamp
from aswanalysis import get_logger


#=====----------------------------------------=====#
# Small helper tools                               #
#=====----------------------------------------=====#
def flatten_nested_dict(d):
    new = dict()
    assert isinstance(d, dict)
    for k,v in d.items():
        if isinstance(v, dict):
            if k in new:
                print(
"WARNING: value for existing key {} being overwritten"
                )
            new.update(flatten_nested_dict(v))
        else:
            new.update({k:v})
    return new

def collect_XY_columns(data, key):
    assert isinstance(data, dict)
    X = sorted(data)
    Y = [data[x][key] for x in X]
    return np.array(X, Y)


_is_globbable = lambda x: x.find('*') >= 0

_create_sequence  = lambda y: [
    single for multi in [
        f if i==0 else list(reversed(f))
        for i,f in enumerate(zip(
        *[s for x in y for s in x.values()]))
    ] for single in multi
]

#_wkf_data_summary_format = ("{0: <7} {1: <12} ",
#"{2: <12} {3: <12} {4: <12} {5: <12} {6: <12} ",
#"{7: <12} {8: <12} {9: <12} {10: <12} {11: <12} ".format(
#"N Tasks", "W Dur", "T Dur,avg", "T Dur,std",
#"T Ini,avg", "T Ini,std",)

#=====----------------------------------------=====#
# Main Component from this module                  #
#=====----------------------------------------=====#
def read_sessions(sessions_home, cfg, loglevel="INFO"):
    logger = get_logger(__name__, "WARNING")
    logger.setLevel(loglevel)
    logger.info("LOGLEVEL set to: {}".format(logger.level))

    #-----------------------------------------------------------#
    # Second thing first, read and set configuration
    config_location = cfg
    logger.info("Reading config from file: %s" % config_location)

    # FIXME with the later proper structure
    with open(config_location, 'r') as f_config:
        analyze_configuration = yaml.safe_load(f_config)

    #workload_filenames = analyze_configuration['workload_filenames']

    #del#workload_folder     = analyze_configuration['workload']['folder']
    workload_filename   = analyze_configuration['workload']['filenames'][0]
    workload_components = analyze_configuration['workload']['timestamps']

    tasks_filenames     = analyze_configuration['task']['filenames']
    tasks_folder        = analyze_configuration['task']['folder']
    task_components     = analyze_configuration['task']['timestamps']
    error_components    = analyze_configuration['task'].get('errors', [])

    #workload_keys     = _unpack_stampkeys(workload_components)
    #task_keys         = _unpack_stampkeys(task_components)
    workload_sequence = _create_sequence(workload_components)
    task_sequence     = _create_sequence(task_components)
    error_connection  = _create_sequence([error_components[0]])
    error_execution   = _create_sequence([error_components[1]])
    error_sequence    = error_connection + error_execution

    #-----------------------------------------------------------#
    # First thing second, handle arguments (use config a bit)
    session_directories = list()

    if os.path.isfile(os.path.join(sessions_home, workload_filename)):
        # Looking to process SINGLE directories
        session_directories.append(sessions_home)

    elif os.path.isdir(sessions_home):
        # Looking to process MANY   directories
        for d in os.listdir(sessions_home):
            if os.path.isfile(os.path.join(sessions_home, d, workload_filename)):
                session_directories.append(os.path.join(sessions_home, d))

        if len(session_directories) == 0:
            logger.info("session_directories must be single folder or top of")
            logger.info("a set of folders containing data for analysis.")
            logger.info("The given value '%s' does not lead to any folders" % 
                sessions_home)
            logger.info("or subfolders with the workflow file '%s' specified" %
                workflow_filename)
            logger.info("by the configuration, Exiting")
            sys.exit(1)

    else:
        logger.info("session_directories must be single folder or top of")
        logger.info("a set of folders containing data for analysis.")
        logger.info("The given value '%s' does not meet the criteria" % 
            sessions_home)
        logger.info("Exiting")
        sys.exit(1)

    #-----------------------------------------------------------#
    # Second thing again, configuring processing setup
    #  - making human readable instead of more general processing
# TODO queryable object
#class timeline(object):
    # TODO incorporate all this into config file
    workloadStart = workload_sequence[0]
    workloadBoot  = workload_sequence[1]
    workloadLive  = workload_sequence[2]
    workloadDone  = workload_sequence[3]
    workloadStop  = workload_sequence[4]

    taskExecStart = task_sequence[0]
    taskExecStop  = task_sequence[-1]
    taskWorkStart = task_sequence[1]
    taskWorkStop  = task_sequence[-2]
    taskTaskStart = task_sequence[2]
    taskTaskStop  = task_sequence[-3]

    interval_timestamp_keys = dict()
    interval_timestamp_keys['workload'] = workload = [workloadStart, workloadStop]
    interval_timestamp_keys['execlive'] = execlive = [taskExecStart, taskExecStop]
    interval_timestamp_keys['worklive'] = worklive = [taskWorkStart, taskWorkStop]
    interval_timestamp_keys['execboot'] = execboot = [workloadStart, taskExecStart]
    interval_timestamp_keys['execinit'] = execinit = [taskExecStart, taskWorkStart]
    interval_timestamp_keys['taskinit'] = taskinit = [taskWorkStart, taskTaskStart]
    interval_timestamp_keys['taskmain'] = taskmain = [taskTaskStart, taskTaskStop]
    interval_timestamp_keys['taskend']  = taskend  = [taskTaskStop, taskWorkStop]
    interval_timestamp_keys['execend']  = execend  = [taskWorkStop, taskExecStop]
    interval_timestamp_keys['workend']  = workend  = [taskExecStop, workloadStop]

    interval_plotlabel_keys = dict()
    interval_plotlabel_keys['workload'] = "Wtotal"
    interval_plotlabel_keys['execlive'] = "Elive"
    interval_plotlabel_keys['worklive'] = "Wlive"
    interval_plotlabel_keys['execboot'] = "Winit"
    interval_plotlabel_keys['execinit'] = "Einit"
    interval_plotlabel_keys['taskinit'] = "Tinit"
    interval_plotlabel_keys['taskmain'] = "Ttotal"
    interval_plotlabel_keys['taskend']  = "Tclose"
    interval_plotlabel_keys['execend']  = "Eclose"
    interval_plotlabel_keys['workend']  = "Wclose"

    n_files_per_task = len(tasks_filenames)

    timestamp_keys = {
        'workload' : workload_sequence,
        'task'     : task_sequence,
        'error'    : error_sequence,
    }

    assert all([_is_globbable(tfnm) for tfnm in tasks_filenames])

    #-----------------------------------------------------------#
    # Third, initialize data structures
    all_timestamps = dict()
    durations = dict()
    analysis = dict()

    #-----------------------------------------------------------#
    # Fourth, read in the workload profile data
    for session_directory in session_directories:

        all_timestamps[session_directory] = timestamp.get_session_timestamps(
            session_directory,
            workload_filename,
            tasks_folder,
            tasks_filenames,
            timestamp_keys,
            convert_to_quantity=True,
        )

    logger.info("All Timestamps:\n{}".format(pformat(all_timestamps)))

    #-----------------------------------------------------------#
    # Fifth, calculate execution durations
    for session_directory, timestamps in all_timestamps.items():
        # TODO aslso assert there is only 1 timestamp
        #      for things like start, stop
        assert len(timestamps['workload']) == 1

        durations[session_directory] = durs = dict()
        [
         durs.update({interval:list()})
         for interval in interval_timestamp_keys
        ]

        for taskstamps in timestamps['task']:
            _stamps = {
                k:v for k,v in timestamps['workload'][0].items()
            }
            _stamps.update(taskstamps)
            for key, (start, stop) in interval_timestamp_keys.items():
                if key in durs and stop in _stamps and start in _stamps:
                    durs[key].append(_stamps[stop][0] - _stamps[start][0])

    #-----------------------------------------------------------#
    # Sixth, calculate duration statistics
    for session_directory, durs in durations.items():
        analysis[session_directory] = anls = dict()

        for interval,duration in durs.items():
            anls[interval] = (np.average(duration), np.std(duration))

    #-----------------------------------------------------------#
    # Seventh, save the analysis
    for session_directory in session_directories:

        durs = durations[session_directory]
        anls = analysis[session_directory]

        timestamps = all_timestamps[session_directory]

        output_profile_path = os.path.join(
            session_directory, "profile.txt")

        output_analysis_path = os.path.join(
            session_directory, "analysis.txt")

        output_timestamps_path = os.path.join(
            session_directory, "timestamps.txt")

        with open(output_profile_path, 'w') as f_out:
            f_out.write(pformat(durs)+'\n')

        with open(output_analysis_path, 'w') as f_out:
            f_out.write(pformat(anls)+'\n')

        with open(output_timestamps_path, 'w') as f_out:
            f_out.write(pformat(timestamps)+'\n')

    #-----------------------------------------------------------#
    # Eighth, FIXME redundant w/ all_timestamps
    #               but try to start at zero
    executors = dict()

    for session_directory, timestamps in all_timestamps.items():

        taskstamps     = timestamps["task"]
        workflowstamps = timestamps["workload"]

        workflow_start = workflowstamps[0][workloadStart][0]

        _executors = list()

        for ts in taskstamps:
            executor = dict()
            for tskey, tslist in ts.items():
                if tslist:
                    if not tskey in executor:
                        executor[tskey] = list() 

                    executor[tskey].extend([
                        timestamp.timestamp_to_workflowtime(t, workflow_start)
                        for t in tslist])

            _executors.append(executor)

        executors[session_directory] = _executors

    #-----------------------------------------------------------#
    # Ninth, classify executors by profile

    # Different Profile Types
    _nostart_      = lambda x: taskWorkStart not in x
    _completed_    = lambda x: all([x.get(tslabel, None) for tslabel in task_sequence])
    _errored_      = lambda x: any([x.get(tse, None) for tse in error_sequence])
    _noerror_      = lambda x: all([not bool(x.get(tse, None)) for tse in error_sequence])
    _disconnected_ = lambda x: bool(x.get(error_connection[0], None))
    _reconnected_  = lambda x: _disconnected_(x) and bool(x.get(error_connection[1],None))

    # Smoosh them for ordered querying
    #  - each sessions' executors will be organized into these groups
    executor_types_list = ["clean","noerror","reconn","disconn","error"]
    executors_by_type = dict()

    # Classify set of executors from each session
    for session_directory, _executors in executors.items():

        # Find members for each profile
        nostart_executors = list(filter(lambda x: _nostart_(x), _executors))

        errored_executors = sorted(
            list(filter(lambda x: _errored_(x), _executors)),
            key=lambda x: x[taskExecStart][0])

        othererror_executors = sorted(
            list(filter(lambda x: _errored_(x) and not _disconnected_(x), _executors)),
            key=lambda x: x[taskExecStart][0])

        disconnected_executors = sorted(
            list(filter(lambda x: _disconnected_(x) and not _reconnected_(x), _executors)),
            key=lambda x: x[taskExecStart][0])

        reconnected_executors = sorted(
            list(filter(lambda x: _reconnected_(x), _executors)),
            key=lambda x: x[taskExecStart][0])

        clean_executors = sorted(
            list(filter(lambda x: _completed_(x) and _noerror_(x), _executors)),
            key=lambda x: x[taskExecStart][0])

        noerror_executors = sorted(
            list(filter(lambda x: _noerror_(x) and not _completed_(x), _executors)),
            key=lambda x: x[taskExecStart][0])

        executors_by_type[session_directory] = {
            "clean"   : clean_executors,  # These 2 should be same
            "noerror" : noerror_executors,# guys, check why if not
            "reconn"  : reconnected_executors,
            "disconn" : disconnected_executors,
            "error"   : othererror_executors,
        }

    #-----------------------------------------------------------#
    # Tenth, plot stuff
    for session_directory, ebt in executors_by_type.items():

        workflow_start = all_timestamps[session_directory]["workload"][0][workloadStart][0]
        workflow_stop  = all_timestamps[session_directory]["workload"][0][workloadStop][0]
        workflow_duration = workflow_stop - workflow_start

        iter_executors = [
            ex for type_ in executor_types_list
            for ex in ebt[type_]
        ]

        n_tasks = len(iter_executors)

        #-----------------------------------------------------------#
        # Nth, preparing for plots
        executors_dispatched = np.concatenate(list(map(
            lambda ex: ex[taskExecStart],
            filter(lambda x: taskExecStart in x, iter_executors)
        )))
        executors_booted = np.concatenate(list(map(
            lambda ex: ex.get(taskWorkStart, [-100]),
            filter(lambda x: taskWorkStart in x, iter_executors)
        )))
        executors_closed = np.concatenate(list(map(
            lambda ex: ex.get(taskWorkStop, [-100]),
            filter(lambda x: taskWorkStop in x, iter_executors)
        )))
        executors_exited = np.concatenate(list(map(
            lambda ex: ex[taskExecStop],
            filter(lambda x: taskExecStop in x, iter_executors)
        )))
        tasks_started = np.array(np.concatenate(list(map(
            lambda ex: list(filter(lambda x: isinstance(x, float),
                ex.get(taskTaskStart, [-100.]))),
            filter(lambda x: taskTaskStart in x, iter_executors)
        ))), dtype=float)
        tasks_completed = np.array(np.concatenate(list(map(
            lambda ex: list(filter(lambda x: isinstance(x, float),
                ex.get(taskTaskStop, [-100.]))),
            filter(lambda x: taskTaskStop in x, iter_executors)
        ))), dtype=float)

        #-----------------------------------------------------------#
        # Nth, doing plots

        minsep = 10 # in seconds, just to help visually
        # Separating the dots if too close, leave be if weird thign
        # that boot is before dispatch (impossible to happen irl)
        pad_2lists = lambda l1,l2: list(map(
            lambda x,y: x if (x-y>minsep or x<y) else y+minsep, l1, l2))

        plot_filename = "timeseries.sort-start.{}.png".format(n_tasks)
        fig = plt.figure(figsize=(8,4))
        ax  = fig.add_subplot(111)

        ax.axvline(
            0, linestyle="--", linewidth=1, color="black",
            label="Workload Start", zorder=1)

        ax.axvline(
            workflow_duration, linestyle="-.", linewidth=1,
            color="black", label="Workload Stop",  zorder=1)

        ax.scatter(
            executors_dispatched, range(len(executors_dispatched)), s=15,
            label="Worker Dispatched", color='blueviolet', zorder=2)

        ax.scatter(
            pad_2lists(executors_booted, executors_dispatched),
            range(len(executors_booted)), s=15,
            label="Worker Booted", color='forestgreen', zorder=3)

        ax.scatter(
            executors_closed, range(len(executors_closed)), s=15,
            label="Worker Closed", color='darkorange', zorder=5)

        ax.scatter(
            pad_2lists(executors_exited, executors_closed),
            range(len(executors_exited)), s=15,
            label="Worker Exited", color='firebrick', zorder=4)

        #-----------------------------------------------------------#
        # MAGIC part to plot bars over each task duration
        #     NOTE bars are shortened to make gaps easier to see
        exid_order = list(range(n_tasks))
        startslabelledyet = 0
        for exid,executor in enumerate(iter_executors):
            try:
                starts = filter(lambda x: isinstance(x, float),
                    executor[taskTaskStart])
                stops  = filter(lambda x: isinstance(x, float),
                    executor[taskTaskStop])

                try:
                    mapped_tasks = list(zip(starts, stops))
                    X = [exid_order.index(exid), exid_order.index(exid)]

                    for mt in mapped_tasks:
                        if startslabelledyet:
                            ax.plot([mt[0]+3, mt[1]-3], X, linewidth=1.5, color='grey', zorder=2)
                        else:
                            startslabelledyet = 1
                            ax.plot([mt[0]+4, mt[1]-4], X, label="Task Runtime", linewidth=1.2, color='grey', zorder=2)

                except Exception as e:
                    print("See if this error preveents plottable data from being processed!")
                    print(e)
                    pass

            except KeyError:
                pass

        ax.set_xlim([-10, workflow_duration + 0.1 * workflow_duration])
        ax.set_ylabel("Worker Index", fontsize=14)
        ax.set_xlabel("Elapsed Time (seconds)", fontsize=14, fontname="Arial")
        ax.set_title("Workload: %s"%session_directory, fontsize=18, fontname="Arial")

        plt.tight_layout()
        legend = plt.legend(fancybox=True, framealpha=1, loc='center', fontsize=9)
        #legend.remove()
        #ax2 = ax.twinx()
        #ax2.add_artist(legend)
        for lh in legend.legendHandles:
            lh._sizes = [20]

    return session_directories, durations, analysis, all_timestamps, executors, executors_by_type


def plot_staircase():
    pass
def plot_taskhist():
    pass
def plot_normalizedlag():
    pass
