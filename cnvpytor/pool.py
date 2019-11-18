""" cnvpytor.pool

Multiprocessing pool functions
"""
import multiprocessing
import logging

_logger = logging.getLogger("cnvpytor.pool")


def _fun(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))


def parmap(f, x_arg, cores=multiprocessing.cpu_count(), info=True):
    """
    Calculates list(map(f,x_args)) using multiprocessing module.

    Parameters
    ----------
    f : callable
    x_arg : list of function arguments
    cores : maximal number of cores

    Returns
    -------
    results : list of values
        List of results.

    """
    cores = min(cores, multiprocessing.cpu_count())
    if info:
        _logger.info("Parallel processing using %d cores" % cores)
    else:
        _logger.debug("Parallel processing using %d cores" % cores)

    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=_fun, args=(f, q_in, q_out))
            for _ in range(cores)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(x_arg)]
    [q_in.put((None, None)) for _ in range(cores)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]
