#!/usr/bin/env python
"""


"""

import numpy as np
import os
import sys
import glob
import scipy.stats as st
import pylab as plt



def process_one_file(fname):
    """Process the result of az_scan_en_norm from az.tcl for a single file

    fname: the {base}.scan file written out by az_scan_en_norm in az.tcl
        We also search for {base}.sim in the same location base {base}.scan resulting from
        az_sim_dchi2 to get the data-based distribution of chi2 to use for calculating
        the probabilities instead of ftest.

    Write the resulting plot to a veusz (https://veusz.github.io/) file with the name
        {path_to_input_dir}/{input_file}.2d.plot
    """

    # do we have a simulation data file for the reference probability density?
    sim_data = None
    if fname[-4:] == 'scan':
        sfile = fname[:-4] + 'sim'
        if len(glob.glob(sfile)) != 0:
            sim_data = np.loadtxt(sfile)


    # read the scanning data of energy-normalization
    data  = np.loadtxt(fname)
    with open(fname, 'r') as fp:
        line = fp.readline()
    stat_0 = np.double(line.split()[2])
    dof_0  = np.int(line.split()[5])

    en_arr = data[:, 0]
    nm_arr = data[:, 1]
    stat   = data[:, 2]

    # turn the 1d arrays to 2d
    nen, nnm = len(np.unique(en_arr)), len(np.unique(nm_arr))
    en_arr = en_arr.reshape((nnm, nen))
    nm_arr = nm_arr.reshape((nnm, nen))
    en     = en_arr[0]
    nm     = nm_arr[:, 0]
    inm    = np.argsort(nm)
    nm     = nm[inm]
    stat   = stat.reshape((nnm, nen))[inm, :]


    ## statistics ##
    # use simulated statistics if available to obtain the significance
    if not sim_data is None:
        sim_dchi2 = (sim_data[:,-1]/2) / (sim_data[:,-2] / (dof_0 - 2))
        obs_dchi2 = ((stat_0 - stat)/2) / (stat / (dof_0 - 2))
        # approximate the simulated distribution with an f-distribution
        sim_pvalue = 1 - st.f(*(st.f.fit(sim_dchi2))).cdf(obs_dchi2)
        #sim_pvalue = np.argmin(np.abs(np.sort(sim_dchi2)[None,None,::-1] - 
        #                              obs_dchi2[:,:,None]), axis=2) / len(sim_dchi2)
        sim_nsigma = -st.norm.ppf(sim_pvalue/2)
        sim_nsigma = np.clip(sim_nsigma, 0, 10)
        sim_prob = 100 * (1 - sim_pvalue)


    # siginifance of features using ftest #
    # fstat = ((c1-c2)/(d1-d2)) / (c2/d2)
    # ftest = 1 - st.f.cdf(fstat, dfn=2, dfd=d2) # 2 for en/norm
    fstat  = ((stat_0 - stat) / 2) / (stat / (dof_0 - 2))
    pvalue = 1 - st.f.cdf(fstat, dfn=2, dfd=dof_0-2)
    prob   = 100*(1 - pvalue)
    nsigma = -st.norm.ppf(pvalue/2)


    from scipy.ndimage import zoom
    nn, xx, yy = zoom(nsigma, 1), zoom(en, 1), zoom(nm, 1)
    lv = [0.5, 1, 2, 3, 4, 5]
    plt.contourf(xx, yy, nn, levels=lv, cmap='GnBu')
    cs = plt.contour(xx, yy, nn, levels=lv, linewidths=0.5, colors='k')
    ylim = 1.2*np.abs(nm[np.any(nsigma > 0.5, 1)][[0,-1]]).max()
    plt.ylim([-ylim, ylim])

    from IPython import embed; embed();exit(0)

    # prepare the output; write to folder veusz in the same location as input file(s)
    basedir = '/'.join(['.'] + fname.split('/')[:-1])

    outfile = '{}/{}.2d.plot'.format(basedir, 
            '.'.join(os.path.basename(fname).split('.')[:-1]))
    # print a useful line that helps when reading the files to veusz
    with open(outfile, 'w') as fp:
        fp.write(('# pvalue nsigma prob ') + (
                    '\n' if sim_data is None else 'pvalue_sim nsigma_sim\n'
                ))
    write_2d_veusz(outfile, pvalue.T, en, nm, 1)
    write_2d_veusz(outfile, nsigma.T, en, nm, 1)
    write_2d_veusz(outfile, prob.T,   en, nm, 1)

    if not sim_data is None:
        write_2d_veusz(outfile, sim_pvalue.T, en, nm, 1)
        write_2d_veusz(outfile, sim_nsigma.T, en, nm, 1)
        write_2d_veusz(outfile, sim_prob.T,   en, nm, 1)


    if sim_data is None:
        return stat, stat_0, dof_0, pvalue, nsigma, nm, en
    else:
        return stat, stat_0, dof_0, [pvalue, sim_pvalue], [nsigma, sim_nsigma], nm, en



def combine_results(res):
    """Combined the statistics from several scans.
    
    res: a list of individual results returned by process_one_file

    """
    # make sure we have more than one file #
    if len(res) < 2:
        raise ValueError('combining data needs more than one file')

    if not np.all([r[0].shape==res[0][0].shape for r in res]):
        raise ValueError('array dimensions of input files do not match')

    # combine data, first by summing the statistics #
    stat   = np.array([r[0] for r in res]).sum(0)
    stat_0 = np.array([r[1] for r in res]).sum(0)
    dof_0  = np.array([r[2] for r in res]).sum(0)
    nm,en  = res[0][5], res[0][6]
    ipvalue = np.array([r[3] for r in res])
    ipvalue_sim = None
    if len(ipvalue.shape) == 4:
        ipvalue_sim = ipvalue[:,1]
        ipvalue = ipvalue[:,0]

    # use Fisher's method to combine pvalues (https://en.wikipedia.org/wiki/Fisher%27s_method)
    cpvalue = 1 - st.chi2.cdf(-2*np.sum(np.log(ipvalue), 0), df=2*len(ipvalue))
    cnsigma = -st.norm.ppf(cpvalue/2)
    cnsigma = np.clip(cnsigma, 0, 10)
    cprob = 100*(1 - cpvalue)
    if not ipvalue_sim is None:
        cpvalue_sim = 1 - st.chi2.cdf(-2*np.sum(np.log(ipvalue_sim), 0), df=2*len(ipvalue))
        cnsigma_sim = -st.norm.ppf(cpvalue_sim/2)
        cnsigma_sim = np.clip(cnsigma_sim, 0, 10)
        cprob_sim   = 100 * (1 - cpvalue_sim)


    # statistics #
    # siginifance of features using ftest #
    # fstat = ((c1-c2)/(d1-d2)) / (c2/d2)
    # ftest = 1 - st.f.cdf(fstat, dfn=2, dfd=d2) # 2 for en/norm
    fstat  = ((stat_0 - stat) / 2) / (stat / (dof_0 - 2))
    pvalue = 1 - st.f.cdf(fstat, dfn=2, dfd=dof_0-2)
    nsigma = -st.norm.ppf(pvalue/2)
    prob   = 100 * (1 - pvalue)

    outfile = 'combined_scan.2d.plot'
    with open(outfile, 'w') as fp:
        fp.write(('# pvalue nsigma prob cpvalue cnsigma cprob ') + (
                    '\n' if ipvalue_sim is None else 
                    'cpvalue_sim cnsigma_sim cprob_sim\n'
                ))
    write_2d_veusz(outfile, pvalue.T, en, nm, 1)
    write_2d_veusz(outfile, nsigma.T, en, nm, 1)
    write_2d_veusz(outfile, prob.T,   en, nm, 1)

    write_2d_veusz(outfile, cpvalue.T, en, nm, 1)
    write_2d_veusz(outfile, cnsigma.T, en, nm, 1)
    write_2d_veusz(outfile, cprob.T,   en, nm, 1)
    if not ipvalue_sim is None:
        write_2d_veusz(outfile, cpvalue_sim.T, en, nm, 1)
        write_2d_veusz(outfile, cnsigma_sim.T, en, nm, 1)
        write_2d_veusz(outfile, cprob_sim.T,   en, nm, 1)

    print('** result saved to {} **'.format(outfile))


def write_2d_veusz(fname, arr, xcent=None, ycent=None, append=False):
    """Write a 2d array to a file for veusz viewing
    
    Args:
        fname: name of file to write
        arr: array to write, shape (len(xcent), len(ycent))
        xcent, ycent: central points of axes.
        append: append to file? Default=False
    """

    thead = '\n\n'
    if xcent is not None and ycent is not None:
        assert( arr.shape==(len(xcent),len(ycent)))
        thead += (  'xcent ' + 
                    ' '.join(['{}'.format(x) for x in xcent]) + '\n')
        thead += (  'ycent ' +
                    ' '.join(['{}'.format(x) for x in ycent]) + '\n')
    arr = arr.T
    txt2d = '\n'.join([
            '{}'.format(' '.join(['{:3.3e}'.format(arr[i, j])
                        for j in range(len(arr[0]))]))
                for i in range(len(arr))
        ])
    with open( fname, 'a' if append else 'w' ) as fp:
        fp.write(thead+txt2d)


if __name__ == '__main__':
    """Process the files resulting from calling scan_en_norm in az.tcl
    If one file is given, process it.
    if more than one file are given, process each one, then produce the
    result from their combination.


    """

    if len(sys.argv) < 2:
        raise ValueError('I need at least a file name')


    # process individual files #
    res = [process_one_file(f) for f in sys.argv[1:]]

    # combine results #
    if len(res) > 1:
        combine_results(res)
    


