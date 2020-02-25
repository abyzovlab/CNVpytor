""" cnvpytor.utils

Misc functions

"""
from __future__ import absolute_import, print_function, division
import numpy as np
from argparse import ArgumentTypeError
from scipy.special import erf
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
import logging
import readline
import requests

_logger = logging.getLogger("cnvpytor.utils")


def gc_at_compress(gc, at):
    """
    Commpress GC/AT content of 100bins using fact that #at+#gc=100 in large majority of bins.
    Transforms #at, #gc -> #at, 100-#at-#gc

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
    cgc = np.array(gc)
    cat = 100 - np.array(at) - cgc
    return np.concatenate((cgc, cat)).astype("uint8")


def gc_at_decompress(gcat):
    """
    Decompress GT/AC content - inverse function of gc_at_compress(gc, at).

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
    return list(gcat[:gcat.size // 2]), list(100 - gcat[:gcat.size // 2] - gcat[gcat.size // 2:])


def gcp_decompress(gcat, bin_ratio=1):
    """
    Decompress GT/AC content and calculate GC percentage.

    Parameters
    ----------
    gcat : numpy.ndarray
        Array contains compressed GC/AT content.

    Returns
    -------
    gcp : list of int
        GC percentage.

    """
    gc, at = gcat[:gcat.size // 2].astype("float"), (100 - gcat[:gcat.size // 2] - gcat[gcat.size // 2:]).astype(
        "float")
    if bin_ratio == 1:
        return 100. * gc / (at + gc + 1e-10)
    n = len(gc)
    his_gc = np.concatenate(
        (np.array(gc), np.zeros(bin_ratio - n + (n // bin_ratio * bin_ratio)))).astype("float")
    his_gc.resize((len(his_gc) // bin_ratio, bin_ratio))
    his_gc = his_gc.sum(axis=1)
    his_at = np.concatenate(
        (np.array(at), np.zeros(bin_ratio - n + (n // bin_ratio * bin_ratio)))).astype("float")
    his_at.resize((len(his_at) // bin_ratio, bin_ratio))
    his_at = his_at.sum(axis=1)

    return 100. * his_gc / (his_at + his_gc + 1e-10)


bp_to_binary = {'A': 0, 'T': 3, 'G': 1, 'C': 2, '.': 4}
binary_to_bp = {0: 'A', 3: 'T', 1: 'G', 2: 'C', 4: '.'}


def snp_compress(pos, ref, alt, nref, nalt, gt, flag, qual):
    """
    Commpress SNP information binary in four arrays:
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
    """
    Decommpress SNP information

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
        ref.append(binary_to_bp[i >> 8])
        alt.append(binary_to_bp[(i >> 5) & int('111', 2)])
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
    return np.array(crd_p), np.array(crd_p) - np.array(crd_u)


def segments_code(segments):
    """
    Convert segments to numpy array e.g. [[1,2],[3]] -> [1,2,MAX,3,MAX]
    Parameters
    ----------
    segments : list of list of int
        Segments e.g. [[1,2],[3]].

    Returns
    -------
    aseg : numpy.ndarra

    """
    max = 2 ** 32 - 1
    l = []
    for s in segments:
        for i in s:
            l.append(i)
        l.append(max)
    return np.array(l, dtype="uint32")


def segments_decode(aseg):
    """
    Decode segments.
    Parameters
    ----------
    aseg : numpy.ndarra of uint32

    Returns
    -------
    segments : list of list of int

    """
    max = 2 ** 32 - 1
    segments = []
    l = []
    for x in list(aseg):
        if x == max:
            segments.append(l)
            l = []
        else:
            l.append(x)
    return segments


def binsize_type(x):
    x = int(x)
    if x % 100 != 0 or x <= 0:
        raise ArgumentTypeError("Bin size should be positive integer divisible by 100!")
    return x


def binsize_format(x):
    if x >= 1000000:
        return str(x // 1000000) + "M"
    elif x >= 1000:
        return str(x // 1000) + "K"
    else:
        return str(x)


def normal_overlap(m1, s1, m2, s2):
    if m1 > m2:
        m1, m2 = m2, m1
        s1, s2 = s2, s1
    a = 1. / (2. * s1 ** 2) - 1. / (2. * s2 ** 2)
    b = m2 / (s2 ** 2) - m1 / (s1 ** 2)
    c = m1 ** 2 / (2 * s1 ** 2) - m2 ** 2 / (2 * s2 ** 2) - np.log(s2 / s1)
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
    if s1 == 0 and s2 == 0:
        return 0.5 * (m1 + m2), 0
    else:
        return (m1 * s2 * s2 + m2 * s1 * s1) / (s1 * s1 + s2 * s2), np.sqrt(s1 * s1 * s2 * s2 / (s1 * s1 + s2 * s2))


def normal(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) / np.sqrt(2 * np.pi) / sigma


def fit_normal(x, y):
    """ Fit Gaussian
    """
    if sum(y) == 0:
        _logger.debug("Problem with fit: all data points have zero value. Return zeros instead fit parameters!")
        return [0, 0, 0], None
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    area = sum(y[:-1] * (x[1:] - x[:-1]))
    _logger.debug("%f %f %f %d" % (area, mean, sigma, len(x)))
    if sigma == 0:
        _logger.debug("Problem with fit: sigma equals zero. Using mean and std instead fitting parameters!")
        return [area, mean, sigma], None

    if len(x) < 3:
        _logger.warning("Problem with fit: insufficient data points. Using mean and std instead fitting parameters!")
        return [area, mean, sigma], None
    try:
        popt, pcov = curve_fit(normal, x, y, p0=[area, mean, sigma])
        popt[2] = np.abs(popt[2])
        if popt[1] <= 0:
            _logger.warning("Problem with fit: negative mean. Using mean and std instead fitting parameters!")
            return [area, mean, sigma], None
        return popt, pcov
    except ValueError:
        _logger.warning("Problem with fit: Value Error. Using mean and std instead fitting parameters!")
        return [area, mean, sigma], None
    except RuntimeError:
        _logger.warning("Problem with fit: Runtime Error. Using mean and std instead fitting parameters!")
        return [area, mean, sigma], None


def t_test_1_sample(mean, m, s, n):
    if s == 0:
        s = 1
    t = (mean - m) / s * np.sqrt(n)
    p = 1.0 - stats.t.cdf(np.abs(t), df=n - 1)
    return p


def t_test_2_samples(m1, s1, n1, m2, s2, n2):
    if s1 == 0:
        s1 = 1
    if s2 == 0:
        s2 = 1
    t = (m1 - m2) / np.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2)
    df = (s1 ** 2 / n1 + s2 ** 2 / n2) ** 2 * (n1 - 1) * (n2 - 1) / (
            s1 ** 4 * (n2 - 1) / n1 ** 2 + s2 ** 4 * (n1 - 1) / n2 ** 2)
    p = 1.0 - stats.t.cdf(np.abs(t), df=int(df + 0.5))
    return p


def getEValue(mean, sigma, rd, start, end):
    aver = np.mean(rd[start:end])
    s = np.std(rd[start:end])
    if s == 0:
        s = sigma * aver / mean if sigma > 0 else 1
    return t_test_1_sample(mean, aver, s, end - start) / (end - start)


def gaussianEValue(mean, sigma, rd, start, end):
    aver = np.mean(rd[start:end])
    max = np.max(rd[start:end])
    min = np.min(rd[start:end])

    if aver < mean:
        x = (max - mean) / (sigma * np.sqrt(2.))
        return np.power(0.5 * (1. + erf(x)), end - start)
    x = (min - mean) / (sigma * np.sqrt(2.))
    return np.power(0.5 * (1. - erf(x)), end - start)


def adjustToEvalue(mean, sigma, rd, start, end, pval, max_steps=1000):
    val = getEValue(mean, sigma, rd, start, end)
    step = 0
    done = False
    while val > pval and not done and step < max_steps:
        done = True
        step += 1
        v1, v2, v3, v4 = 1e10, 1e10, 1e10, 1e10
        if start > 0:
            v1 = getEValue(mean, sigma, rd, start - 1, end)
        if end - start > 2:
            v2 = getEValue(mean, sigma, rd, start + 1, end)
            v3 = getEValue(mean, sigma, rd, start, end - 1)
        if end < len(rd):
            v4 = getEValue(mean, sigma, rd, start, end + 1)
        if min(v1, v2, v3, v4) < val:
            done = False
            if v1 == min(v1, v2, v3, v4):
                start -= 1
                val = v1
            elif v2 == min(v1, v2, v3, v4):
                start += 1
                val = v2
            elif v3 == min(v1, v2, v3, v4):
                end -= 1
                val = v3
            elif v4 == min(v1, v2, v3, v4):
                end += 1
                val = v4
    if val <= pval:
        return start, end
    return None


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


def beta(k, m, p, phased=False):
    """
    Returns likelihood beta function f(p) where 0 <= p <= 1.
    Function is not normalized.

    Parameters
    ----------
    k : int
        Number of haplotype 1 reads.
    m : int
        Number of haplotype 2 reads.
    p : float
        Parameter.
    phased : bool
        Likelihood will be symmetrised if not phased.

    Returns
    -------
    f : float
        Value od likelihood function,.

    """
    if k == m or phased:
        return 1.0 * p ** k * (1. - p) ** m
    else:
        return 1.0 * p ** k * (1. - p) ** m + 1.0 * p ** m * (1. - p) ** k


def log_beta(k, m, p, phased=False):
    """
    Returns logarithm of likelihood beta function f(p) where 0 <= p <= 1.
    Function is not normalized.

    Parameters
    ----------
    k : int
        Number of haplotype 1 reads.
    m : int
        Number of haplotype 2 reads.
    p : float
        Parameter.
    phased : bool
        Likelihood will be symmetrised if not phased.

    Returns
    -------
    f : float
        Value od likelihood function,.

    """
    if k == m or phased:
        return np.log(1.0 * p ** k * (1. - p) ** m)
    else:
        return np.log(1.0 * p ** k * (1. - p) ** m + 1.0 * p ** m * (1. - p) ** k)


def likelihood_overlap(lk1, lk2):
    """
    Returns overlap area of two likelihood functions.

    Parameters
    ----------
    lk1 : numpy.ndarray
        First likelihood function.
    lk2 : numpy.ndarray
        Second likelihood function.

    Returns
    -------
    overlap : float
        Overlap area.

    """
    return np.sum(np.min((lk1, lk2), axis=0))


def decode_position(s):
    """

    Parameters
    ----------
    s : str

    Returns
    -------
    str
    """
    return int(s.replace("K", "000").replace("k", "000").replace("M", "000000").replace("m", "000000"))


def decode_region(s):
    """

    Parameters
    ----------
    s : str

    Returns
    -------
    list of tuples

    """
    regs = s.split(",")
    ret = []
    for r in regs:
        chr_interval = r.split(":")
        if len(chr_interval) > 1:
            begin_end = chr_interval[1].split("-")
            ret.append((chr_interval[0], (decode_position(begin_end[0]), decode_position(begin_end[1]))))
        else:
            ret.append((chr_interval[0], (1, 1000000000)))
    return ret


def likelihood_baf_pval(likelihood):
    """
    Calculates baf level and p-value for given likelihood function.
    Parameters
    ----------
    likelihood

    Returns
    -------
    b : float
        BAF level (difference from 1/2)
    p : float
        p-value for event different than 1/2

    """
    res = likelihood.size
    max_lh = np.amax(likelihood)
    ix = np.where(likelihood == max_lh)[0][0]
    if ix > res // 2:
        ix = res - 1 - ix
    b = 1.0 * (res // 2 - ix) / (res + 1)

    ix1 = (res // 2 + ix) // 2
    ix2 = res - 1 - ix1
    p = np.sum(likelihood[ix1:ix2]) / np.sum(likelihood)
    return b, p

def likelihood_pixels_pval(likelihood):
    """
    Calculates maximum likelihood pixel positions and p-value.
    Parameters
    ----------
    likelihood

    Returns
    -------
    b : float
        BAF level (difference from 1/2)
    p : float
        p-value for event different than 1/2

    """
    res = likelihood.size
    max_lh = np.amax(likelihood)
    ix = np.where(likelihood == max_lh)[0][0]
    if ix > res // 2:
        ix = res - 1 - ix
    ix1 = (res // 2 + ix) // 2
    ix2 = res - 1 - ix1
    p = np.sum(likelihood[ix1:ix2]) / np.sum(likelihood)
    return ix, res - 1 - ix, p

def is_downloadable(url):
    """
    Returns does the url contain a downloadable resource
    Parameters
    ----------
    url : str
        Resource url

    Returns
    -------
    downloadable : bool
        True if the url contain a downloadable resource

    """
    h = requests.head(url, allow_redirects=True)
    header = h.headers
    content_type = header.get('content-type')
    if 'text' in content_type.lower():
        return False
    if 'html' in content_type.lower():
        return False
    return True


def download(url, file):
    """
    Download resource into file

    Parameters
    ----------
    url : str
        Resource url
    file : str
        Filename

    """
    r = requests.get(url, allow_redirects=True)
    with open(file, 'wb') as f:
        f.write(r.content)


class PromptCompleter:
    def __init__(self, command_tree):
        """
        Initialize completer class for GNU readline tab comleter.

        Parameters
        ----------
        command_tree : dictionary
            Example: command_tree = {"comannd1":None,"command2":{"subcommand2.1":None,subcommand2.1:None}}

        """
        self.command_tree = command_tree

    def _traverse(self, tokens, tree):
        if tree is None or len(tokens) == 0:
            return []
        if len(tokens) == 1:
            return [x + ' ' for x in tree if x.startswith(tokens[0])]
        else:
            if tokens[0] in tree.keys():
                return self._traverse(tokens[1:], tree[tokens[0]])
        return []

    def complete(self, text, state):
        try:
            tokens = readline.get_line_buffer().split()
            if not tokens or readline.get_line_buffer()[-1] == ' ':
                tokens.append("")
            results = self._traverse(tokens, self.command_tree) + [None]
            return results[state]
        except Exception as e:
            print(e)


class TerminalColor:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN2 = '\033[92m'
    GREEN = '\033[32m'
    YELLOW = '\033[93m'
    YELLOW2 = '\033[33m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def add_tabs(s, n=4):
    return "\n".join(list(map(lambda x: " " * n + x, s.split("\n"))))


def key_val_str(d, indent=4):
    s = ""
    for i in d:
        s += " " * indent + "* " + i + " - " + d[i] + "\n"
    return s


def help_format(topic="", p_desc="", p_usage="", p_type="", p_default="", p_affects="", p_example="", p_see=""):
    ret_str = "\n"
    if p_desc != "":
        ret_str += TerminalColor.BOLD + topic + "\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_desc) + TerminalColor.END + "\n\n"
    if p_usage != "":
        ret_str += TerminalColor.BOLD + "USAGE\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_usage) + TerminalColor.END + "\n\n"
    if p_type != "":
        ret_str += TerminalColor.BOLD + "TYPE\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_type) + TerminalColor.END + "\n\n"
    if p_default != "":
        ret_str += TerminalColor.BOLD + "DEFAULT\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_default) + TerminalColor.END + "\n\n"
    if p_affects != "":
        ret_str += TerminalColor.BOLD + "PLOTS AFFECTS\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_affects) + TerminalColor.END + "\n\n"
    if p_example != "":
        ret_str += TerminalColor.BOLD + "EXAMPLE(s)\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_example) + TerminalColor.END + "\n\n"
    if p_see != "":
        ret_str += TerminalColor.BOLD + "SEE ALSO\n" + TerminalColor.END
        ret_str += TerminalColor.DARKCYAN + add_tabs(p_see) + TerminalColor.END + "\n\n"
    return ret_str[:-1]
