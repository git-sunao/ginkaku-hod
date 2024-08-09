import numpy as np
import matplotlib.pyplot as plt
from dark_emulator import model_hod
import meascorr
import pymultinest
import time
import json
import os
import argparse
from mpi4py import MPI
here = os.path.dirname(os.path.abspath(__file__))

class Sampler:
    def __init__(self, zbin, use_model=True):
        # Read the data
        ng     = np.loadtxt(os.path.join(here,'../data/ng_sig_z{:03d}.dat'.format(zbin)))
        cov_ng = np.loadtxt(os.path.join(here,'../data/ng_cov_10p_z{:03d}_z{:03d}.dat'.format(zbin, zbin)))
        rp, wp = np.loadtxt(os.path.join(here,'../data/wp_sig_z{:03d}.dat'.format(zbin)), unpack=True)
        cov_wp = np.loadtxt(os.path.join(here,'../data/wp_cov_z{:03d}_z{:03d}.dat'.format(zbin, zbin)))
        icov_wp= np.linalg.inv(cov_wp)
        zl_rep = [0.26, 0.51, 0.63][zbin]
        pimax  = 100.0
        dlnR   = np.log(rp[1]/rp[0])
        self.data = ng, cov_ng, rp, wp, cov_wp, icov_wp, zl_rep, pimax, dlnR

        if not use_model: return None
        # Instantiate the objects
        ## HOD model
        hod = model_hod.darkemu_x_hod()
        ## Cosmological parameters (this is fixed to Planck 2018)
        cparam = np.array([0.02225,  0.1198 ,  0.6844 ,  3.094  ,  0.9645 , -1.])
        hod.set_cosmology(cparam)
        ## Measurement correction
        mc = meascorr.wp_meascorr_class({'Om':0.233, 'wde':-1.0})
        mc.set_param({'Om':1.0-cparam[2], 'wde':cparam[5]})
        self.models = hod, mc

    def sample_with_multinest(self, prior, names_to_sample, output):
        # unpack the data
        ng, cov_ng, rp, wp, cov_wp, icov_wp, zl_rep, pimax, dlnR = self.data

        def mn_prior(cube, ndims, nparams):
            for i, name in enumerate(names_to_sample):
                vmin, vmax = prior[name]
                cube[i] = vmin + (vmax-vmin)*cube[i]

        Neval = 0
        def mn_likelihood(cube, ndim, nparams):
            t0 = time.time()
            # get the prediction
            wp_pred, ng_pred = self.get_prediction(cube, names_to_sample)

            # compute chi2
            chi2_wp = (wp-wp_pred)@icov_wp@(wp-wp_pred)
            chi2_ng = (ng-ng_pred)**2/cov_ng # Note: ng is a scalar
            chi2_tot= chi2_wp + chi2_ng

            # return the likelihood
            loglike = -0.5*chi2_tot

            nonlocal Neval
            print('LIKE evaluation, Neval: {} time: {} loglike: {}'.format(Neval, time.time()-t0, loglike))
            Neval += 1

            return loglike

        # MCMC sampling
        pymultinest.run(mn_likelihood, mn_prior, len(names_to_sample), 
            outputfiles_basename=output, resume=False, verbose=True)
        json.dump(names_to_sample, open(output + 'params.json', 'w'))

    def get_prediction(self, params, names):
        # unpack the data
        ng, cov_ng, rp, wp, cov_wp, icov_wp, zl_rep, pimax, dlnR = self.data

        # models
        hod, mc = self.models
        
        # HOD parameters
        gparam = gparam_fiducial.copy()
        for i, name in enumerate(names):
            gparam[name] = params[i]
        
        # Set the parameters
        hod.set_galaxy(gparam)

        # get model prediction
        f_pimax, f_rp = mc.get_corrs(zl_rep)
        wp_pred = hod.get_wp(f_rp*rp, zl_rep, pimax=pimax*f_pimax, rsd=True, dlnrp=dlnR) # Mpc/h
        f_ng = mc.get_ng_corr(zl_rep)
        ng_pred = f_ng * hod.get_ng(zl_rep) # (Mpc/h)^-3

        return wp_pred, ng_pred

    def load_chain(self, fname_chain_post_equal):
        samples = np.loadtxt(fname_chain_post_equal)
        names_to_sample = json.load(open(fname_chain_post_equal.replace('post_equal_weights.dat', 'params.json')))
        return samples, names_to_sample

    def derived_predictions(self, fname_chain_post_equal):
        samples, names_to_sample = self.load_chain(fname_chain_post_equal)

        print('The samples size: ', samples.shape[0])

        # parallelize the computation
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        sub_samples = samples[rank::size]

        sub_samples_sig = []
        sub_samples_hod = []
        for sample in sub_samples:
            # signal
            wp, ng = self.get_prediction(sample, names_to_sample)
            sub_samples_sig.append(np.append(wp, ng))
            # hod
            Mh = self.models[0].Mh
            Nc = self.models[0].Ncen
            Ns = self.models[0].Nsat
            sub_samples_hod.append(np.append(Nc, Ns))

        # gather the results
        samples_sig = comm.gather(sub_samples_sig, root=0)
        samples_hod = comm.gather(sub_samples_hod, root=0)
        if rank == 0:
            # signal
            samples_sig = np.concatenate(samples_sig)
            output = fname_chain_post_equal.replace('.dat', '-derived_signal.dat')
            np.savetxt(output, samples_sig)
            # hod
            samples_hod = np.concatenate(samples_hod)
            output = fname_chain_post_equal.replace('.dat', '-derived_hod.dat')
            np.savetxt(output, samples_hod)
            # halo mass bins
            np.savetxt(fname_chain_post_equal.replace('.dat', '-Mh.dat'), Mh)

    def get_bestfit_predictions(self, fname_chain_post_equal):
        samples, names_to_sample = self.load_chain(fname_chain_post_equal)
        sample_beftfit = samples[np.argmax(samples[:,-1]), :]
        wp, ng = self.get_prediction(sample_beftfit, names_to_sample)
        Mh, Nc, Ns = self.models[0].Mh, self.models[0].Ncen, self.models[0].Nsat
        return wp, ng, Mh, Nc, Ns

prior   =  {'logMmin': [12.0, 14.5], 
            'sigma_sq': [0.01, 1.0], 
            'logM1':[12.0, 16.0], 
            'alpha':[0.5, 3.0], 
            'kappa':[0.01, 3.0], 
            'poff':[0.0, 1.0], 
            'Roff':[0.01, 1.0], 
            'alpha_inc':[0, 10], 
            'logM_inc':[12.0, 15.3]}

# HOD parameters (This is updated during MCMC sampling)
gparam_fiducial={"logMmin": 13.13, 
                "sigma_sq":0.22, 
                "logM1": 14.21, 
                "alpha": 1.13, 
                "kappa": 1.25, # HOD parameters
                "poff": 0.2, # off-centering parameters p_off is the fraction of off-centered galaxies. 
                "Roff": 0.1, # Roff is the typical off-centered scale with respect to R200m.
                "sat_dist_type": "NFW",
                "alpha_inc": 0.44, 
                "logM_inc": 13.57} # incompleteness parameters. For details, see More et al. (2015)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('zbin', type=int, help='Redshift bin')
    parser.add_argument('output', type=str, help='Output file name')
    args = parser.parse_args()

    sampler = Sampler(args.zbin)
    output = os.path.join(here, '../chains/{}-z{}-'.format(args.output, args.zbin))
    # sampler.sample_with_multinest(prior, ['logMmin', 'sigma_sq', 'logM1', 'alpha', 'kappa'], output)
    sampler.derived_predictions(output+'post_equal_weights.dat')
