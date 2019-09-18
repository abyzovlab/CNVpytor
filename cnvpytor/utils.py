""" cnvpytor.utils

Misc functions

"""
from __future__ import absolute_import, print_function, division
import numpy as np
from argparse import ArgumentTypeError
from scipy.stats import norm
from scipy.optimize import curve_fit
import logging

_logger = logging.getLogger("cnvpytor.utils")


def gc_at_compress(gc, at):
    """ Commpress GC/AT content of 100bins using fact that #at+#gc=100 in large majority of bins.
    If #at+#gc=100 it will store just one byte (#gc) otherwise it will store two bytes: (128+#gc and #at).
    First bit is used as flag: if it is set there are additional byte for #at.

    Parameters
    ----------
    gc : list of int
        Binned GC content (100bp bins).
    at : list of int
        Binned AT content (100bp bins).

    Returns
    -------
    gcat : numpy.ndarray
        Array contains compressed GC/AT content.

    """
    gcat = []
    for g, a in zip(gc, at):
        if a + g == 100:
            gcat.append(g)
        else:
            gcat.append(g + 128)
            gcat.append(a)
    return np.array(gcat, dtype="uint8")


def gc_at_decompress(gcat):
    """ Decompress GT/AC content - inverse function of gc_at_compress(gc, at).

    Parameters
    ----------
    gcat : numpy.ndarray
        Array contains compressed GC/AT content.

    Returns
    -------
    gc : list of int
        Binned GC content (100bp bins).
    at : list of int
        Binned AT content (100bp bins).

    """
    fil = filter(lambda x: (x[2] < 128), zip(list(gcat), list(gcat[1:]) + [0], [0] + list(gcat[:-1])))
    y = list(map(lambda x: (x[0] - 128, x[1]) if x[0] > 127 else (x[0], 100 - x[0]), fil))
    return [list(c) for c in zip(*y)]


def gcp_decompress(gcat):
    """ Decompress GT/AC content and calculate GC percentage.

    Parameters
    ----------
    gcat : numpy.ndarray
        Array contains compressed GC/AT content.

    Returns
    -------
    gcp : list of int
        GC percentage.

    """
    fil = filter(lambda x: (x[2] < 128), zip(list(gcat), list(gcat[1:]) + [0], [0] + list(gcat[:-1])))
    return list(map(lambda x: float(x[0] - 128) / (x[0] + x[1] - 128 + 1e-10) if x[0] > 127 else float(x[0]), fil))


bp_to_binary = {'A': 0, 'T': 3, 'G': 1, 'C': 2, '.': 4}
binary_to_bp = {0: 'A', 3: 'T', 1: 'G', 2: 'C', 4: '.'}


def snp_compress(pos, ref, alt, nref, nalt, gt, flag, qual):
    """ Commpress SNP information binary in four arrays:
            snp_pos - diferences in positions
            snp_counts - 32 bits: ref count (16 bits), alt count (16 bits)
            snp_desc -  16 bits: ref (3 bits), alt (3 bits), gt (3 bits), flag (2 bits)
            snp_qual - 8 bits
    """
    snp_pos = list(map(lambda i, j: i - j,
                       pos, [0] + pos[:-1]))
    snp_counts = list(map(lambda i, j: (i << 16) + j,
                          nref, nalt))
    snp_desc = list(map(lambda r, a, g, f: ((((((bp_to_binary[r] << 3) + bp_to_binary[a]) << 3) + g) << 2) + f),
                        ref, alt, gt, flag))
    snp_qual = qual
    return np.array(snp_pos, dtype="uint32"), np.array(snp_desc, dtype="uint16"), np.array(snp_counts,
                                                                                           dtype="uint32"), np.array(
        snp_qual, dtype="uint8")


def snp_decompress(snp_pos, snp_desc, snp_counts, snp_qual):
    """ Decommpress SNP information
    """
    pos = []
    s = 0
    for i in snp_pos:
        s += i
        pos.append(s)
    ref = []
    alt = []
    gt = []
    flag = []
    for i in snp_desc:
        ref.append(i >> 8)
        alt.append((i >> 5) & int('111', 2))
        gt.append((i >> 2) & int('111', 2))
        flag.append(i & int('11', 2))
    nref = []
    nalt = []
    for i in snp_counts:
        nref.append(i >> 16)
        nalt.append(i & int('ffff', 16))
    qual = list(snp_qual)

    return pos, ref, alt, nref, nalt, gt, flag, qual


def mask_compress(mask):
    """ Commpress strict mask P value intervals.
    Argument:
        mask - list of intervals [(start_1, end_1), (start_2, end_2),...]
    """
    pos = [p for interval in mask for p in interval]
    cmask = list(map(lambda i, j: i - j, pos, [0] + pos[:-1]))
    return np.array(cmask, dtype="uint32")


def mask_decompress(cmask):
    """ Decommpress strict mask P value intervals.
    """
    pos = []
    s = 0
    for i in cmask:
        s += i
        pos.append(s)
    return list(zip(pos[::2], pos[1::2]))


def rd_compress(rd_p, rd_u, data_type="uint16"):
    """ Compres RD signal
    """
    return rd_p.astype(data_type), (rd_p - rd_u).astype(data_type)


def rd_decompress(crd_p, crd_u):
    """ Decommpress SNP information
    """
    return np.array(crd_p), np.array(crd_p) + np.array(crd_u)


def binsize_type(x):
    x = int(x)
    if x % 100 != 0 or x <= 0:
        raise ArgumentTypeError("Bin size should be positive integer divisible by 100!")
    return x


def normal_overlap(m1, s1, m2, s2):
    if m1 > m2:
        m1, m2 = m2, m1
        s1, s2 = s2, s1
    a = 1. / (2. * s1 ** 2) - 1. / (2. * s2 ** 2)
    b = m2 / (s2 ** 2) - m1 / (s1 ** 2)
    c = m1 ** 2 / (2 * s1 ** 2) - m2 ** 2 / (2 * s2 ** 2) - np.log(s2 / s1)
    print(np.roots([a, b, c]))
    roots = sorted(np.roots([a, b, c]))
    if len(roots) == 0:
        return 1
    elif len(roots) == 1:
        r = roots[0]
        return norm.cdf(r, m2, s2) + (1. - norm.cdf(r, m1, s1))
    elif s1 < s2:
        return norm.cdf(roots[1], m2, s2) - norm.cdf(roots[0], m2, s2) + 1 - norm.cdf(roots[1], m1, s1) + norm.cdf(
            roots[0], m1, s1)
    return norm.cdf(roots[1], m1, s1) - norm.cdf(roots[0], m1, s1) + 1 - norm.cdf(roots[1], m2, s2) + norm.cdf(
        roots[0], m2, s2)


def normal_merge(m1, s1, m2, s2):
    return (m1 * s2 * s2 + m2 * s1 * s1) / (s1 * s1 + s2 * s2), np.sqrt(s1 * s1 * s2 * s2 / (s1 * s1 + s2 * s2))


def normal(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def fit_normal(x, y):
    """ Fit Gaussian
    """
    if sum(y) == 0:
        _logger.debug("Problem with fit: all data points have zero value. Return zeros instead fit parameters!")
        return [0, 0, 0], None
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    if sigma == 0:
        _logger.debug("Problem with fit: sigma equals zero. Return zeros instead fit parameters!")
        return [0, 0, 0], None
    area = sum(y[:-1] * (x[1:] - x[:-1]))
    try:
        popt, pcov = curve_fit(normal, x, y, p0=[area, mean, sigma])
        return popt, pcov
    except ValueError:
        _logger.warning("Problem with fit: optimal parameters not found. Using mean and std instead!")
        return [area, mean, sigma], None
    except RuntimeError:
        _logger.warning("Problem with fit: optimal parameters not found. Using mean and std instead!")
        return [area, mean, sigma], None


def calculate_gc_correction(his_rd_gc, mean, sigma, bin_size=1):
    """ Calculate GC correction from RD-GC histogram
    """
    max_bin = int(max(2 * mean, mean + 5 * sigma) / bin_size)
    his = his_rd_gc[1:max_bin, :]
    rd = np.repeat((np.arange(1, max_bin * bin_size, bin_size)[:max_bin - 1]).reshape((max_bin - 1, 1)), 101, axis=1)
    np.seterr(divide='ignore', invalid='ignore')
    gc_corr = np.sum(rd * his, axis=0) / np.sum(his, axis=0)
    gc_corr[np.isnan(gc_corr)] = 1.0
    gc_corr = gc_corr / (np.sum(gc_corr * np.sum(his, axis=0)) / np.sum(np.sum(his, axis=0)))
    return gc_corr
