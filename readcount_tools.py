import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd
import anndata
from scipy import sparse

SMALL_SIZE = 7
MEDIUM_SIZE = 7
BIGGER_SIZE = 7
LINEWIDTH=1
POINTSIZE=1
POINTSIZE_HIGHLIGHT=3
POINTSIZE_SMALL=0.5
TICKLENGTH=3
LEGEND_FONTSIZE=5.5
PAGEWIDTH_IN = 6.25

SPINEWIDTH=0.5
LETTER_LOC_X = -0.30
LETTER_LOC_Y = 0.95

PAPER_CONTEXT = {'font.size':SMALL_SIZE,
'axes.titlesize':SMALL_SIZE,
'axes.labelsize':MEDIUM_SIZE,
'xtick.labelsize':SMALL_SIZE,
'ytick.labelsize':SMALL_SIZE,
'legend.fontsize':SMALL_SIZE,
'figure.titlesize':BIGGER_SIZE,
'xtick.major.width':SPINEWIDTH,   #tick thickness
'ytick.major.width':SPINEWIDTH,
'xtick.major.size':TICKLENGTH,     #tick length
'ytick.major.size':TICKLENGTH,
"figure.facecolor":(1.0, 1.0, 1.0, 1.0),
"axes.facecolor":(1.0, 1.0, 1.0, 1.0)}

def add_largedot_legend(ax,loc,kwargs={},fix_alpha=False):
    #allow kwargs to overwrite default
    frameon=True
    if 'frameon' in kwargs.keys():
        frameon = kwargs.pop('frameon')
        
    lgnd = ax.legend(loc=loc,frameon=frameon,**kwargs)
    for l in lgnd.legendHandles:
        l._sizes = [30]
        if fix_alpha:
            l.set_alpha(1)

def compute_gene_stats(ad,suffix=''):
    denseX = ad.X.A
    ad.var[f'gene_var{suffix}'] = np.var(denseX,axis=0)
    ad.var[f'gene_mean{suffix}'] = np.mean(denseX,axis=0)
    ad.var[f'gene_FF{suffix}'] = ad.var[f'gene_var{suffix}']/ad.var[f'gene_mean{suffix}']
    ad.var[f'gene_fraction_zeros{suffix}'] = 1 - np.count_nonzero(denseX,axis=0)/denseX.shape[0]

def compute_kde(x,y,xscale,yscale,weights=None):
    xtrans = np.log10(x) if xscale == 'log' else x
    ytrans = np.log10(y) if yscale == 'log' else y
    data = np.vstack((xtrans,ytrans))
    kde = gaussian_kde(data,weights=weights)
    logdensity = kde.logpdf(data)
    return logdensity

def compute_marginals(counts):
    '''compute depths per cell (ns) and relative expression fractions per gene (ps)'''
    ns = np.sum(counts,axis=1)
    ps = np.sum(counts,axis=0)
    ps = ps / np.sum(ps)    
    return np.squeeze(np.array(ns)), np.squeeze(np.array(ps))

def get_tag(alpha,theta,clipping=True):
    if clipping:
        return f'pr_theta{theta}_alpha{alpha:.1f}'
    else:
        return f'pr_theta{theta}_alpha{alpha:.1f}_unclipped'
    
def pearson_residuals_compound(counts, theta, alpha, clipping=True):
    '''
    Computes analytical residuals for NB model with a fixed theta
    
    `theta=np.Inf` corresponds to Poisson
    
    `clipping=True` will clip outlier residuals to sqrt(N)
    
    `alpha=1` corresponds to a non-compound NB/Poisson
    Other `alpha` values correspond to a compound model with the following specs:
    
    number of unique molecules/UMIs:
        `k ~ NB(mu, theta)`
        
    amplified/sequenced counts of the i-th molecule/UMI :
        `Z_i ~ unknown distribution`
        
    total counts per gene:
        `X ~ sum_{i=1..k} (Z_i)`
        
    From the general variance of the compund process X
    
        `var[X] = E[X] * (E[Z] + var[Z]/E[Z])`
        
    this leads to the following for the NB compound we have here:
    
        `var[X] = E[X] * alpha + E[X]^2/theta`
        
    where
    
        `alpha = E[Z] + ( Var[Z]/E[Z] )`
    '''
    
    counts_sum0 = np.sum(counts, axis=0, keepdims=True)
    counts_sum1 = np.sum(counts, axis=1, keepdims=True)
    counts_sum  = np.sum(counts)

    #get residuals
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(alpha*mu + mu**2/theta)


    #clip to sqrt(n)
    if clipping:
        n = counts.shape[0]
        z[z >  np.sqrt(n)] =  np.sqrt(n)
        z[z < -np.sqrt(n)] = -np.sqrt(n)
    
    return z

def plot_mean_vs(adatas,axes, xkey = 'mean', ykey = 'ff', yscale='log', colorby='',vmin=None,vmax=None,
                    kde=True, autotitle=True, kde_cmap=None,s=3 ):
   
    for i,(ad,ax) in enumerate(zip(adatas,axes.flatten())):

        plt.sca(ax)        
        x=ad.var[xkey]
        y=ad.var[ykey]
        if kde:
            logdensity = compute_kde(x,y,xscale='log',yscale=yscale)
            ax.scatter(x,y,s=s,linewidth=0,c=logdensity,cmap=kde_cmap,rasterized=True)
        elif colorby:
            ax.scatter(x,y,s=s,linewidth=0,c=ad.var[colorby],cmap=kde_cmap,rasterized=True,vmin=vmin,vmax=vmax)
        else:
            ax.scatter(x,y,s=s,linewidth=0,rasterized=True)
        ax.set_xscale('log')
        ax.set_yscale(yscale)
        if autotitle:
            plt.title(ad.uns['protocol'])
                
                
    sns.despine()
    
    return axes


def broken_zeta(a1=1.4,
                a2=8.0,
                breakpoint=100,
                size=10000,
                z_max=100000,
                seed=42,
                return_p=False):    
    '''
    Samples from a broken zeta as follows:
    
        z ~ BrokenZeta(a1,a2,b)
    
        p(z) ~ z**-a1                   if z<b
        p(z) ~ b**(-a1+a2) * z**-a2     otherwise
    
    z > z_max are ignored when computing the PMF and thus cannot occur by design.
    '''
    
    np.random.seed(seed)
    
    z = np.arange(1,z_max)
    p = np.zeros(z.shape)
    
    below_b_idx = z<breakpoint
    above_b_idx = ~below_b_idx
    p[below_b_idx] = z[below_b_idx]**-a1
    p[above_b_idx] = breakpoint**(-a1+a2) * z[above_b_idx]**-a2
    p = p/np.sum(p)
    if return_p:
        return np.random.choice(z, size=size, p=p), p
    else:
        return np.random.choice(z, size=size, p=p)

    
def zeta_params_to_str(zeta_params):
    return f'a1{zeta_params["a1"]}_a2{zeta_params["a2"]}_b{zeta_params["breakpoint"]}_zmax{zeta_params["z_max"]}'

def molsim_params_to_str(molsim_params):
    return f'theta{molsim_params["theta_molecules"]}_n{molsim_params["depth"]}'

def params_to_pretty_str(molsim_params,zeta_params,amplification_stats):
    pretty_str = fr'molecules simulated with $\theta={molsim_params["theta_molecules"]}$, '\
                 fr'$n={molsim_params["depth"]}$' '\n'\
                 fr'amplified by BrokenZeta(a1={zeta_params["a1"]}, a2={zeta_params["a2"]}, '\
                 fr'b={zeta_params["breakpoint"]}, zmax={zeta_params["z_max"]})' '\n' \
                 fr'leading to E[Z]={amplification_stats["mean"]:.0f} and FF[Z]={amplification_stats["ff"]:.0f}'
    return pretty_str
    
def simulate_readcounts(molsim_params,zeta_params,amplification_seed=42,tag='untagged',color='tab:blue'):
    
    molecules_sim,ps_molsim_input, ns_observed, ps_observed = simulate_molecules(**molsim_params)
    
    amplification_params = dict(zeta_params=zeta_params,
                            seed=amplification_seed)

    if zeta_params['constant']:
        #UMI case
        readcounts_sim = molecules_sim.copy()
        amplification_stats = dict(mean=1,median=1,var=0,ff=0,alpha=1,max=1)
    else:
        #non UMI case
        readcounts_sim,amplification_stats = simulate_amplification(molecules_sim,
                                        **amplification_params)
    
    ad = anndata.AnnData(X=sparse.csc_matrix(readcounts_sim),layers=dict(molecules=molecules_sim))
    
    ad.var['ps_molsim_input']=ps_molsim_input
    ad.var['ps_observed']=ps_observed
    ad.obs['ns_observed']=ns_observed

    ad.uns["clustername"] = '_'.join (('simulated_readcounts', tag, molsim_params_to_str(molsim_params), zeta_params_to_str(zeta_params)))
    ad.uns['clustercolor'] = color
    ad.uns['simulation_info_pretty'] = params_to_pretty_str(molsim_params,zeta_params,amplification_stats)
    ad.var['gene_mean_withinCluster'] = np.mean(ad.X.A,axis=0)
    ad.var['genes'] = [f'simulated_gene_{i}' for i in range(ad.shape[1])]
    ad.var.set_index('genes',inplace=True,drop=False)
    ad.var.index.name = 'gene_name' #needed to be able to "drop=False" when saving to h5ad
    
    ad.uns['zeta_params'] = zeta_params
    ad.uns['molecule_sim_simple_params'] = molsim_params
    ad.uns['amplification_params'] = amplification_params
    ad.uns['amplification_stats'] = amplification_stats
    params_sim = dict(**molsim_params,**zeta_params,**amplification_stats,zeta_seed=amplification_seed)
    
    return ad,params_sim

def simulate_molecules(n_cells, ps_input, depth=100000, theta_molecules=10, seed=42):
    #simulate marginals
    ns_constant = depth*np.ones(n_cells)    
    #simulate molecules
    mu_molecules = ns_constant[:, np.newaxis] @ ps_input[np.newaxis]
    p_molecules = theta_molecules / (theta_molecules + mu_molecules)
    np.random.seed(seed)
    if np.isinf(theta_molecules):
        #poisson case
        molecules_simulated = np.random.poisson(mu_molecules)
    else:
        #NB case
        molecules_simulated = np.random.negative_binomial(theta_molecules, p_molecules)
    
    ns_observed,ps_observed =compute_marginals(molecules_simulated)
    zero_genes_idx = ps_observed==0
    print(f'removing {sum(zero_genes_idx)} all-zero genes after simulation')
    molecules_simulated=molecules_simulated[:,~zero_genes_idx]
    ps_input_after_filter = ps_input[~zero_genes_idx]
    ps_observed_after_filter = ps_observed[~zero_genes_idx]
    
    return molecules_simulated, ps_input_after_filter, ns_observed, ps_observed_after_filter


def simulate_amplification(molecules,
                           zeta_params=dict(a1=1.4,
                                                 a2=8.0,
                                                 breakpoint=100,
                                                 z_max=100000),
                           seed=42,):
    '''
    Simulates amplification of input molecules / UMIs by BrokenZeta compound model
    
        z~BrokenZeta()
        readcounts = sum(z)
        
    `return_z_stats`
    '''
    
    molecules=molecules.astype(int)
    readcounts=np.zeros(molecules.shape)
    
    #work only on nonzero (gene-x-cell)-observations 
    molecules_nonzero_idx = molecules>0
    molecules_nonzero_counts = molecules[molecules_nonzero_idx]
    split_idx = np.cumsum(molecules_nonzero_counts).astype(int)
    
    #one large sample instead of separate ones is more efficient (one sample per molecule observed)
    if 'constant' in zeta_params.keys():
        constant_flag = zeta_params.pop('constant')
        
    zs = broken_zeta(**zeta_params,size=int(sum(molecules_nonzero_counts)),seed=seed)
    
    zeta_params['constant'] = constant_flag

    mean=np.mean(zs)
    maxx=np.max(zs)
    median=np.median(zs)
    var=np.var(zs)
    ff=var/mean
    empirical_alpha=mean+ff
    
    print('Broken Zeta amplification with',zeta_params)
    print('Effectively amplifying with Zs that have mean=%.1f, median=%.1f, var=%.1f, FF=%.1f, leading to alpha=%.1f' % (mean,median,var,ff,empirical_alpha))
    
    #splitting up into separate groups of samples for each (gene-x-cell)-observation
    zs_per_cell_x_gene=np.split(zs,split_idx[:-1])
    
    #summing reads for each (gene-x-cell)-observation
    zs_sums =[sum(z) for z in zs_per_cell_x_gene]
    
    #mapping back to count matrix
    readcounts[molecules_nonzero_idx]=zs_sums   
    

    return readcounts, dict(mean=mean,median=median,var=var,ff=ff,alpha=empirical_alpha,max=maxx)



