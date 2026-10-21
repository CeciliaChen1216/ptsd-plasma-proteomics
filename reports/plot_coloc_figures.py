"""Render publication Figures S2/S9 from frozen tables; no association fitting.

Usage: python reports/plot_coloc_figures.py --output /path/to/figures
Dependencies: pandas, numpy, matplotlib. R equivalents are supplied alongside.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=Path('reports_coloc'))
args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
root=Path(__file__).resolve().parents[1]
snap=root/'results_snapshot'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})

e=pd.read_csv(snap/'sensitivity_review_20261019/evidence_matrix.csv',keep_default_na=False)
genes=['AKT3','CGREF1','FES','KHK','SNX18','CD40','FURIN','UBE2L6','SIRPA','CD101']
domains=e.domain.drop_duplicates().tolist()
labels=['Primary\nMR','Steiger\ndirectionality','MR-PRESSO\nassessment','Reverse-MR\nassessment','Multi-SNP\nIVW','eQTL-SMR\n+ direction','Colocalization','Brain pQTL\n/ PWAS']
colors={'S':'#2C7FB8','C':'#D95F0E','I':'#BDBDBD','N':'#F0F0F0'}
fig,ax=plt.subplots(figsize=(12,7.5));fig.subplots_adjust(left=.10,right=.985,top=.84,bottom=.11)
for row in e.itertuples(index=False):
    x=domains.index(row.domain);y=genes.index(row.protein)
    ax.add_patch(Rectangle((x-.48,y-.48),.96,.96,facecolor=colors[row.state],edgecolor='white'))
    if row.annotation:ax.text(x,y,row.annotation,ha='center',va='center',fontsize=9,color='white' if row.state in ('S','C') else '#283442')
ax.set_xlim(-.5,7.5);ax.set_ylim(9.5,-.5)
ax.set_xticks(range(8),labels,fontsize=9);ax.xaxis.tick_top()
ax.set_yticks(range(10),genes,fontweight='bold');ax.tick_params(length=0,pad=10)
for s in ax.spines.values():s.set_visible(False)
fig.suptitle('Evidence matrix for ten PTSD-prioritized proteins',x=.10,y=.965,ha='left',fontsize=15,fontweight='bold')
ax.legend(handles=[Patch(facecolor=colors[k],label=t) for k,t in zip('SCIN',['Supportive','Cautionary','Inconclusive','Not evaluable'])],
    loc='upper center',bbox_to_anchor=(.5,-.045),ncol=4,frameon=False,fontsize=10)
for ext in ('png','pdf'):fig.savefig(args.output/f'Supplementary_Figure_S2.{ext}',dpi=400)
plt.close(fig)

p=snap/'coloc_full_region_20261020'
su=pd.read_csv(p/'coloc_SuSiE_prior_grid.csv');ab=pd.read_csv(p/'coloc_ABF_prior_grid.csv');q=pd.read_csv(p/'summary_counts.csv').set_index('protein')
mx=su.groupby(['protein','N_mode','p12'],as_index=False)['PP.H4.abf'].max()
order=['FES','FURIN','CGREF1','KHK','UBE2L6','CD40','AKT3','SIRPA','SNX18','CD101']
fig,axes=plt.subplots(5,2,figsize=(10.2,12),sharex=True,sharey=True)
fig.subplots_adjust(left=.09,right=.98,bottom=.07,top=.86,hspace=.60,wspace=.18)
for ax,g in zip(axes.flat,order):
    a=ab[ab.protein==g].sort_values('p12')
    ax.plot(a.p12,a['PP.H4.abf'],color='#C46B19',marker='^',markersize=4,lw=1.5)
    for mode,col,ls,marker in [('NTOT','#2166AC','-','o'),('NEFF','#00897B','--','o')]:
        z=mx[(mx.protein==g)&(mx.N_mode==mode)].sort_values('p12')
        if len(z):ax.plot(z.p12,z['PP.H4.abf'],color=col,ls=ls,marker=marker,markersize=5 if mode=='NEFF' else 3.5,markerfacecolor='none' if mode=='NEFF' else col,lw=1.2)
    ax.axvline(1e-5,color='#AAAAAA',ls=':',lw=.8);ax.axhline(.8,color='#AAAAAA',ls='--',lw=.8)
    ax.set_xscale('log');ax.set_xlim(8e-8,1.3e-4);ax.set_ylim(-.025,1.04)
    ax.set_xticks([1e-7,1e-6,1e-5,1e-4],['10⁻⁷','10⁻⁶','10⁻⁵','10⁻⁴']);ax.minorticks_off()
    ax.set_yticks([0,.2,.4,.6,.8,1]);ax.tick_params(labelbottom=True,labelsize=9)
    ax.set_title(f'{g}  ({int(q.loc[g,"final_variants"]):,} SNPs)',fontweight='bold',fontsize=11,pad=23)
    if g not in set(su.protein):
        text='SuSiE: protein fit did not converge' if g=='SIRPA' else 'SuSiE: no PTSD credible set'
        ax.text(.5,1.045,text,transform=ax.transAxes,ha='center',va='bottom',fontsize=8,color='#555555')
    if ax in axes[:,0]:ax.set_ylabel('Posterior H4',fontsize=10)
fig.suptitle('Dense regional colocalization and prior sensitivity',fontsize=15,fontweight='bold',y=.975)
fig.legend(handles=[Line2D([0],[0],color='#2166AC',marker='o',label='SuSiE NTOT'),Line2D([0],[0],color='#00897B',ls='--',marker='o',markerfacecolor='none',label='SuSiE NEFF'),Line2D([0],[0],color='#C46B19',marker='^',label='Single-variant ABF')],loc='upper center',bbox_to_anchor=(.5,.953),ncol=3,frameon=False)
fig.text(.5,.025,'Shared-variant prior p₁₂ (p₁ = p₂ = 10⁻⁴)',ha='center',fontsize=11)
fig.text(.5,.006,'SuSiE curves show the maximum across signal pairs. Missing curves are unavailable, not zero.',ha='center',fontsize=8.5)
for ext in ('png','pdf'):fig.savefig(args.output/f'Supplementary_Figure_S9.{ext}',dpi=400)
plt.close(fig)
print('Rendered S2 and S9 from the frozen Oct20 tables.')
