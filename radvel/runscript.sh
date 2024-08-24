python betelgeuse_radvel_mcmc.py BG_rv10.csv --pmin 2000 --pmax 2200 --gp_length 300 --fix_gp_length --do_mcmc --nrun 15000 --emax 0.6 --niter 1 --gp_kernel Exp
python betelgeuse_radvel_mcmc.py BG_rv10.csv --pmin 2000 --pmax 2200 --gp_length 300 --fix_gp_length --do_mcmc --nrun 15000 --drop_sample --drop_option 6 --circular --niter 1 --gp_kernel Exp
