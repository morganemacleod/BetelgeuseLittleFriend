import os

import matplotlib
import numpy as np
import pylab as pl
import pandas as pd
from scipy import optimize
import corner
import random
import radvel
import radvel.likelihood
from radvel.plot import orbit_plots, mcmc_plots
import argparse
from astropy.time import Time

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',help='rv filename to process',type=str)
    parser.add_argument('--pmin',default=1600.0,help='pmin [d]',type=float)
    parser.add_argument('--pmax',default=2400.0,help='pmax [d]',type=float)
    parser.add_argument('--fix_gp_length',action=argparse.BooleanOptionalAction)
    parser.add_argument('--no_gp',action=argparse.BooleanOptionalAction)
    parser.add_argument('--circular',action=argparse.BooleanOptionalAction)
    parser.add_argument('--drop_sample',action=argparse.BooleanOptionalAction)
    parser.add_argument('--no_jit',action=argparse.BooleanOptionalAction)
    parser.add_argument('--gp_length',default=5000.0,help='GP length [d]',type=float)
    parser.add_argument('--nwalkers',default=20,type=int)
    parser.add_argument('--nensembles',default=8,type=int)
    parser.add_argument('--nrun',default=1000,type=int)
    parser.add_argument('--emax',default=0.6,type=float)
    parser.add_argument('--drop_option',default=1,type=int)
    parser.add_argument('--niter',default=1,type=int)
    parser.add_argument('--gp_pow',default=1.0,type=float)
    parser.add_argument('--gp_kernel',default="Exp",type=str)
    parser.add_argument('--gp_amp',default=1.0,type=float)
    parser.add_argument('--fix_gp_amp',action=argparse.BooleanOptionalAction)
    parser.add_argument('--do_mcmc',action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    print("args = ",args,"\n")
    
    # SET PARAMS
    myKernel = args.gp_kernel
    hnames = radvel.gp.KERNELS[myKernel]
    print(myKernel,hnames)

    eccentric = np.where(args.circular, False, True)
    use_GP = np.where(args.no_gp, False, True)
    fix_amplitude=np.where(args.fix_gp_amp,True,False)
    
    tc1_guess = 2459945
    vfactor = 1e3 #km/s ->m/s

    def initialize_model_GP(eccentric=True):
        params = radvel.Parameters(1,basis='per tc secosw sesinw logk') 
        params['per1'] = radvel.Parameter(value=(args.pmax + args.pmin)/2, mcmcscale=10 )
        params['tc1'] = radvel.Parameter(value=tc1_guess,mcmcscale=30 )
        params['secosw1'] = radvel.Parameter(value=0.0,vary=eccentric )
        params['sesinw1'] = radvel.Parameter(value=0.0,vary=eccentric )
        params['logk1'] = radvel.Parameter(value=np.log(2*vfactor)  )
        params['gamma'] = radvel.Parameter(value=21*vfactor, vary=True)
        if args.no_jit:
            params['jit'] = radvel.Parameter(value=0,vary=False)
        else:
            params['jit'] = radvel.Parameter(value=1.0*vfactor)
        params['dvdt'] = radvel.Parameter(value=0,vary=False)
        params['curv'] = radvel.Parameter(value=0,vary=False)
        if use_GP:
            if myKernel=='SqExp':
                params['gp_length'] = radvel.Parameter(value=args.gp_length*1.01,vary=True)
                params['gp_amp'] = radvel.Parameter(value=2.5*vfactor,vary=True)
            if myKernel=='RatQuad':
                params['gp_length'] = radvel.Parameter(value=args.gp_length*1.01,vary=True)
                params['gp_amp'] = radvel.Parameter(value=1*vfactor,vary=True)
                params['gp_pow'] = radvel.Parameter(value=args.gp_pow,vary=False)
            if myKernel=='ExpPow':
                params['gp_length'] = radvel.Parameter(value=args.gp_length*1.01,vary=True)
                params['gp_amp'] = radvel.Parameter(value=args.gp_amp*vfactor,vary=(not fix_amplitude) )
                params['gp_pow'] = radvel.Parameter(value=args.gp_pow,vary=False)
            if myKernel=='Exp':
                params['gp_length'] = radvel.Parameter(value=args.gp_length*1.01,vary=True)
                params['gp_amp'] = radvel.Parameter(value=args.gp_amp*vfactor,vary=True )
 
        mod = radvel.RVModel(params)
        return mod

    # Read data


    for ii in range(args.niter):
    
        # drop sample
        fn = args.filename 
        rv = pd.read_csv(fn)
        if args.drop_sample:
            if args.drop_option==1:
                ind_keep = sorted( random.sample(range(len(rv)), int(0.8*len(rv))) )
            elif args.drop_option==2:
                ind_keep = (rv.JD < Time(1935,format='jyear').jd) | (rv.JD > Time(1980,format='jyear').jd)
            elif args.drop_option==3:
                ind_keep = (rv.JD < Time(2008,format='jyear').jd)
            elif args.drop_option==4:
                ind_keep = (rv.JD < Time(1960,format='jyear').jd)
            elif args.drop_option==5:
                ind_keep = (rv.JD > Time(1960,format='jyear').jd)
            elif args.drop_option==6:
                ind_keep = (rv.JD < Time(2014,format='jyear').jd)
            else:
                print("DROP OPTION NOT RECOGNIZED")
    
            t = np.array(rv.JD[ind_keep])
            vel = np.array(rv.RV[ind_keep])*vfactor 
            errvel = np.array(rv.RV_err[ind_keep])*vfactor 
        else:
            t = np.array(rv.JD)
            vel = np.array(rv.RV)*vfactor 
            errvel = np.array(rv.RV_err)*vfactor 
    

        print("drop sampling: ",args.drop_sample)
        print("using this length tables:" ,len(t),len(vel), len(errvel) )


        ti = np.linspace(t[0]-365,t[-1]+365,1000)


        # setup model
        gpmod = initialize_model_GP(eccentric=eccentric)
        if use_GP:
            gplike = radvel.likelihood.GPLikelihood(gpmod, t, vel, errvel,
                                                hnames=hnames,kernel_name=myKernel)
        else:
            gplike = radvel.likelihood.RVLikelihood(gpmod,t,vel,errvel)

        # likelihood
        print("-------------- likelihood ---------------")
        print(gplike)
        print("-----------------------------------------")

        # posterior
        gppost = radvel.posterior.Posterior(gplike)
        gppost.priors += [radvel.prior.Gaussian('gamma',21.1*vfactor,3*vfactor)]
        gppost.priors += [radvel.prior.HardBounds( 'per1',args.pmin,args.pmax)]
        gppost.priors += [radvel.prior.HardBounds( 'tc1',tc1_guess-1500,tc1_guess+1500)]
        if args.no_jit == None:
            gppost.priors += [radvel.prior.Jeffreys( 'jit', 0.01*vfactor, 3*vfactor)]
            gppost.priors += [radvel.prior.HardBounds( 'logk1', np.log(0.1*vfactor), np.log(10*vfactor) ) ]
        if eccentric:
            gppost.priors += [radvel.prior.EccentricityPrior(1, upperlims=args.emax )]
        if use_GP:
            if args.fix_gp_amp:
                gppost.priors += [radvel.prior.Gaussian('gp_amp',args.gp_amp*vfactor,0.01*args.gp_amp*vfactor)]
            else:
                gppost.priors += [radvel.prior.Jeffreys('gp_amp', 0.01*vfactor, 20.*vfactor)]
            if args.fix_gp_length:
                gppost.priors += [radvel.prior.Gaussian('gp_length',args.gp_length,1)]
            else:
                gppost.priors += [radvel.prior.HardBounds('gp_length',args.gp_length,30000)]
        

        # optimize
        res  = optimize.minimize(
            gppost.neglogprob_array,     # objective function is negative log likelihood
            gppost.get_vary_params(),    # initial variable parameters
            method='Powell', # Powell also works
            options={'disp':True}
        )


        print("-------------- posterior ---------------")
        print(gppost)
        print('BIC=',gppost.bic())
        print("-----------------------------------------")

        #plot_results(gplike)    
        rv_resid = (gplike.y-gplike.predict(gplike.x)[0]-gplike.params['gamma'].value-gplike.model(gplike.x))/vfactor,
        rv_resid_std = np.sqrt(gplike.yerr**2 + gpmod.params['jit'].value**2 + gplike.predict(gplike.x)[1]**2)/vfactor

        mybins = np.linspace(-5,5,21)
        hv,hb = np.histogram(rv_resid[0]/rv_resid_std,bins=mybins,density=False)
    
        from astropy.modeling import models, fitting

        gauss_fit = fitting.LevMarLSQFitter()

        gm  = gauss_fit(models.Gaussian1D(), (hb[1:]+hb[0:-1])/2,hv ) 
        print('\n gm= \n',gm)


    
        # MCMC
        if args.do_mcmc:
            gpchains = radvel.mcmc(gppost,nwalkers=args.nwalkers,nrun=args.nrun,ensembles=args.nensembles,
                                   maxGR=1.01,minAfactor=40,maxArchange=0.03)
            gppost.chains=gpchains
        
        
        
        
        # WRITE POSTERIOR TO FILE
        
        gppost.argparse = args

        
        drop_str = str(np.where(args.drop_sample,"drop"+str(args.drop_option)+"_","_"))
        circ_str = str(np.where(args.circular,"circ_","e"+str(args.emax)+"_" ))
        jit_str = str(np.where(args.no_jit,"_nojit_","_"))
        gp_str = str(args.gp_length)    
        iter_str = str(ii)
        save_fn = "BG_radvel_MLE_"+myKernel+gp_str+fn[2:-4]+"_p"+str(args.pmin)+"_"+str(args.pmax)+jit_str+circ_str+drop_str+"_iter_"+iter_str+".pkl"
        print("SAVING TO:",save_fn)
        gppost.writeto(save_fn)


