#!/usr/bin/env python


import numpy as np
import os
import sys
import glob
import aztools as az
import scipy.stats as st
import pylab as plt


def process_one_file(fname):
    """Process the result of scan_en_norm from az.tcl for a single file
    If fname has a form {base}.scan, we search for {base}.sim resulting from
    az_sim_dchi2 to get the data-based distribution of chi2 to use instead of ftest
    """

    sim_data = None
    if fname[-4:] == 'scan':
        sfile = fname[:-4] + 'sim'
        if len(glob.glob(sfile)) != 0:
            sim_data = np.loadtxt(sfile)


    
    data  = np.loadtxt(fname)
    with open(fname, 'r') as fp:
        line = fp.readline()
    stat_0 = np.double(line.split()[2])
    dof_0  = np.int(line.split()[5])

    en_arr = data[:, 0]
    nm_arr = data[:, 1]
    stat   = data[:, 2]


    nen, nnm = len(np.unique(en_arr)), len(np.unique(nm_arr))
    en_arr = en_arr.reshape((nnm, nen))
    nm_arr = nm_arr.reshape((nnm, nen))
    en     = en_arr[0]
    nm     = nm_arr[:, 0]
    inm    = np.argsort(nm)
    nm     = nm[inm]
    stat   = stat.reshape((nnm, nen))[inm, :]
    ipos   = nm > 0
    ineg   = nm < 0


    # statistics #
    if not sim_data is None:
        sim_dchi2 = (sim_data[:,-1]/2) / (sim_data[:,-2] / (dof_0 - 2))
        obs_dchi2 = ((stat_0 - stat)/2) / (stat / (dof_0 - 2))
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


    basedir = '/'.join(['.'] + fname.split('/')[:-1]) + '/veusz'
    os.system('mkdir -p {}'.format(basedir))

    outfile = '{}/{}.2d.plot'.format(basedir, 
            '.'.join(os.path.basename(fname).split('.')[:-1]))
    with open(outfile, 'w') as fp:
        fp.write(('# 1,2,3sigma: 2.296 6.18 11.829\n'
                  '# stat pvalue nsigma prob stat_pos pvalue_pos nsigma_pos prob_pos '
                    'stat_neg pvalue_neg nsigma_neg prob_neg ') + (
                    '\n' if sim_data is None else 
                    'pvalue_sim nsigma_sim pvalue_pos_sim nsigma_po_sim pvalue_neg_sim pvalue_neg_sim\n'
                ))
    az.misc.write_2d_veusz(outfile, stat.T,   xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, pvalue.T, xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, nsigma.T, xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, prob.T,   xcent=en, ycent=nm, append=1)

    # positive norm #
    az.misc.write_2d_veusz(outfile, stat[ipos].T,   xcent=en, ycent=nm[ipos], append=1)
    az.misc.write_2d_veusz(outfile, pvalue[ipos].T, xcent=en, ycent=nm[ipos], append=1)
    az.misc.write_2d_veusz(outfile, nsigma[ipos].T, xcent=en, ycent=nm[ipos], append=1) 
    az.misc.write_2d_veusz(outfile, prob[ipos].T,  xcent=en, ycent=nm[ipos], append=1) 

    # negative norm #
    az.misc.write_2d_veusz(outfile, stat[ineg].T,   xcent=en, ycent=-nm[ineg], append=1)
    az.misc.write_2d_veusz(outfile, pvalue[ineg].T, xcent=en, ycent=-nm[ineg], append=1)
    az.misc.write_2d_veusz(outfile, nsigma[ineg].T, xcent=en, ycent=-nm[ineg], append=1) 
    az.misc.write_2d_veusz(outfile, prob[ineg].T,   xcent=en, ycent=-nm[ineg], append=1) 

    if not sim_data is None:
        az.misc.write_2d_veusz(outfile, sim_pvalue.T, xcent=en, ycent=nm, append=1)
        az.misc.write_2d_veusz(outfile, sim_nsigma.T, xcent=en, ycent=nm, append=1)
        az.misc.write_2d_veusz(outfile, sim_prob.T,   xcent=en, ycent=nm, append=1)

        # positive norm #
        az.misc.write_2d_veusz(outfile, sim_pvalue[ipos].T, xcent=en, ycent=nm[ipos], append=1)
        az.misc.write_2d_veusz(outfile, sim_nsigma[ipos].T, xcent=en, ycent=nm[ipos], append=1) 
        az.misc.write_2d_veusz(outfile, sim_prob[ipos].T, xcent=en, ycent=nm[ipos], append=1) 

        # negative norm #
        az.misc.write_2d_veusz(outfile, sim_pvalue[ineg].T, xcent=en, ycent=-nm[ineg], append=1)
        az.misc.write_2d_veusz(outfile, sim_nsigma[ineg].T, xcent=en, ycent=-nm[ineg], append=1) 
        az.misc.write_2d_veusz(outfile, sim_prob[ineg].T, xcent=en, ycent=-nm[ineg], append=1) 

    if sim_data is None:
        return stat, stat_0, dof_0, pvalue, nsigma, nm, en, ipos, ineg
    else:
        return stat, stat_0, dof_0, [pvalue, sim_pvalue], [nsigma, sim_nsigma], nm, en, ipos, ineg

    # add constrains from the simulations


def combine_results(res):
    """Combined the statistics from several scans.
    
    res: a list of individual results returned by process_one_file

    """
    # make sure we have more than one file #
    if len(res) < 2:
        raise ValueError('combining data needs more than one file')

    if not np.all([r[0].shape==res[0][0].shape for r in res]):
        raise ValueError('files properties do not match')

    # combine data #
    stat   = np.array([r[0] for r in res]).sum(0)
    stat_0 = np.array([r[1] for r in res]).sum(0)
    dof_0  = np.array([r[2] for r in res]).sum(0)
    nm,en  = res[0][5], res[0][6]
    ipos   = res[0][7]
    ineg   = res[0][8]
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
        fp.write(('# stat pvalue nsigma prob stat_pos pvalue_pos nsigma_pos prob_pos '
                    'stat_neg pvalue_neg nsigma_neg prob_neg cpvalue cnsigma cprob ') + (
                    '\n' if ipvalue_sim is None else 
                    'cpvalue_sim cnsigma_sim cprob_sim\n'
                ))
    az.misc.write_2d_veusz(outfile, stat.T,   xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, pvalue.T, xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, nsigma.T, xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, prob.T,   xcent=en, ycent=nm, append=1)

    # positive norm #
    az.misc.write_2d_veusz(outfile, stat[ipos].T,   xcent=en, ycent=nm[ipos], append=1)
    az.misc.write_2d_veusz(outfile, pvalue[ipos].T, xcent=en, ycent=nm[ipos], append=1)
    az.misc.write_2d_veusz(outfile, nsigma[ipos].T, xcent=en, ycent=nm[ipos], append=1) 
    az.misc.write_2d_veusz(outfile, prob[ipos].T,   xcent=en, ycent=nm[ipos], append=1) 

    # negative norm #
    az.misc.write_2d_veusz(outfile, stat[ineg].T,   xcent=en, ycent=-nm[ineg], append=1)
    az.misc.write_2d_veusz(outfile, pvalue[ineg].T, xcent=en, ycent=-nm[ineg], append=1)
    az.misc.write_2d_veusz(outfile, nsigma[ineg].T, xcent=en, ycent=-nm[ineg], append=1) 
    az.misc.write_2d_veusz(outfile, prob[ineg].T,   xcent=en, ycent=-nm[ineg], append=1) 

    az.misc.write_2d_veusz(outfile, cpvalue.T, xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, cnsigma.T, xcent=en, ycent=nm, append=1)
    az.misc.write_2d_veusz(outfile, cprob.T,   xcent=en, ycent=nm, append=1)
    if not ipvalue_sim is None:
        az.misc.write_2d_veusz(outfile, cpvalue_sim.T, xcent=en, ycent=nm, append=1)
        az.misc.write_2d_veusz(outfile, cnsigma_sim.T, xcent=en, ycent=nm, append=1)
        az.misc.write_2d_veusz(outfile, cprob_sim.T,   xcent=en, ycent=nm, append=1)

    print('** result saved to {} **'.format(outfile))


if __name__ == '__main__':
    """Process the files resulting from calling scan_en_norm in az.tcl
    If one file is given, process it.
    if more than one file are given, process each one, then produce the
    result from the combination.


    """

    if len(sys.argv) < 2:
        raise ValueError('I need at least a file name')


    # process individual files #
    res = [process_one_file(f) for f in sys.argv[1:]]

    # combine results #
    if len(res) > 1:
        combine_results(res)
    


