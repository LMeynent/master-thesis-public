import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings
import stan
import arviz as az

from datetime import datetime
from argparse import ArgumentParser

from scipy.stats import gamma, weibull_min, lognorm
import scipy.integrate as integrate

def discretise_pdf(pdf, n=100):
    res = np.zeros(n)
    res[0] = quad(pdf, 0, 1.5)[0]
    for i in range(1, n):
        res[i] = quad(pdf, i+0.5, i+1.5)[0]
    return res / np.sum(res)

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

FILES = {
        'Hovedstaden': 'hope_data/Rt-CAPITAL_chains=5_foldername=2023-04-13_samples=1000_thinning=1_tree=8_warmup=1000.csv',
        'Midtjylland': 'hope_data/Rt-CENTRAL_chains=5_foldername=2023-04-13_samples=1000_thinning=1_tree=8_warmup=1000.csv',
        'Nordjylland': 'hope_data/Rt-NORTH_chains=5_foldername=2023-04-13_samples=1000_thinning=1_tree=8_warmup=1000.csv',
        'Sjælland': 'hope_data/Rt-ZEALAND_chains=5_foldername=2023-04-13_samples=1000_thinning=1_tree=8_warmup=1000.csv',
        'Syddanmark': 'hope_data/Rt-SOUTH_chains=5_foldername=2023-04-13_samples=1000_thinning=1_tree=8_warmup=1000.csv'
}

with open('./combined_model.stan', 'r') as ifh:
    STAN_CODE = ifh.read()

def gamma_from_mean_sd(mean, sd):
    return (mean**2)/(sd**2), (sd**2)/mean

def gamma_from_mean_cv(mean, cv):
    k = 1/np.sqrt(cv)
    return k, mean/k
    
if __name__ == '__main__':
    
    parser = ArgumentParser(description='MCMC-based epidemiologic modelling using hospitalisation data')
    parser.add_argument('-c', '--chains', type=int, default=4)
    parser.add_argument('-s', '--samples', type=int, default=1000)
    parser.add_argument('-t', '--maxtreedepth', type=int, default=5)
    parser.add_argument('-l', '--long', action='store_true')
    
    args = parser.parse_args()
    
    n_chains = args.chains
    n_samples = args.samples
    max_tree_depth = args.maxtreedepth
    long = args.long

    W = 40   # convolution window size for g and pi

    df = pd.read_csv('../hope_data/epidata.csv')
    df['date'] = pd.to_datetime(df['date'])

    begin = '2020-03-01' if long else '2020-08-01'
    end   = '2021-12-01' if long else '2021-02-22'
    
    H = df.groupby(['date', 'region'], as_index=False)['hospit'].agg('mean')\
          .pivot(index='date', columns='region', values='hospit')\
          .loc[begin:end]
    
    D = df.groupby(['date', 'region'], as_index=False)['deaths'].agg('mean')\
          .pivot(index='date', columns='region', values='deaths')\
          .loc[begin:end]
    
    # g ~ Weibull(shape=2.826, scale=5.665)
    g = discretise_pdf(lambda x: weibull_min.pdf(x, c=2.826, scale=5.665), n=40)

    # i2o ~ LogNormal(μ = 1.621, σ = 0.418)
    i2o = discretise_pdf(lambda x: lognorm.pdf(x, s=0.418, scale=np.exp(1.621)), n=40)

    # o2h ~ Weibull(shape=0.845, scale=5.506)
    o2h = discretise_pdf(lambda x: weibull_min.pdf(x, c=0.845, scale=5.506), n=40)
    
    # o2d ~ Gamma(shape=4.938, scale=3.605)
    o2d = discretise_pdf(lambda x: gamma.pdf(x, a=4.938, scale=3.605), n=40)
    
    Pi = np.convolve(i2o, o2h)[:W]
    Tau = np.convolve(i2o, o2d)[:W]

    Pi /= np.sum(Pi)
    Tau /= np.sum(Tau)
    
    STAN_DATA = {
        'd': H.shape[0],
        'r': H.shape[1],
        'w': W,
        'H': H.values.astype('int'),
        'D': D.values.astype('int'),
        'I0': np.ones((H.shape[1]), dtype='int') if long else df[df['date'] == '2020-08-01'].groupby('region')['cases'].mean().values,
        'Pi': Pi,
        'Tau': Tau,
        'g': g
    }

    posterior = stan.build(STAN_CODE, data=STAN_DATA, random_seed=2020)
    
    fit = posterior.sample(num_chains=n_chains, 
                           num_samples=n_samples,
                           max_depth=max_tree_depth,
                           delta=0.99
                          )
    
    stan_data = az.from_pystan(
        posterior=fit,
        observed_data=["H", "D"]
    )
    stan_data.to_netcdf('{}_combined_mcmc_stan.netcdf'.format(
        datetime.now().strftime("%y%m%d_%H%M")
    ))
