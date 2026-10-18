#!/usr/bin/env python3
"""Review plot with explicit unsupported radii and control-only error comparison."""
import argparse
import csv
import json
from pathlib import Path
import platform

import numpy as np
from matplotlib import pyplot as plt
from download_pilot import file_sha256 as sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    with (a.analysis/'recovery_groups.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for h,color in [(100,'#286a9d'),(250,'#c06b20')]:
        for mode,style in [('raw','o-'),('corrected','s--')]:
            group=sorted([r for r in rows if r['cohort']=='area_qualified' and r['kron_factor']=='3.5'
                          and float(r['scale_kpc'])==h and float(r['amplitude_at_100kpc'])==3
                          and r['extraction']==mode],key=lambda r:float(r['r_min_kpc']))
            x=np.array([(float(r['r_min_kpc'])+float(r['r_max_kpc']))/2 for r in group])
            y=np.array([1+float(r['median_fractional_bias']) for r in group])
            n=np.array([int(r['n_eligible']) for r in group]);good=n>=10
            axes[0].plot(x,y,style,color=color,label=f'h={h} kpc; '+('raw' if mode=='raw' else 'outer plane'))
            axes[0].scatter(x[~good],y[~good],s=70,facecolors='white',edgecolors=color,zorder=5)
            if h==100 and mode=='raw':
                for xx,nn in zip(x,n):axes[0].text(xx,1.19,f'n={nn}',ha='center',fontsize=8)
    axes[0].axhline(1,color='k',lw=.7);axes[0].axhspan(.9,1.1,alpha=.08,color='gray')
    axes[0].axvspan(300,500,color='gray',alpha=.15);axes[0].text(400,.76,'n < 10\nexploratory',ha='center',fontsize=9)
    axes[0].set_ylim(.74,1.22);axes[0].set_xlim(50,500)
    axes[0].set_xlabel('Radius (kpc)');axes[0].set_ylabel('Re-detected increment / geometric-ring truth')
    axes[0].legend(fontsize=8,loc='lower left')
    with np.load(a.analysis/'covariances.npz') as z:
        for name,label in [('all_k3.5_raw','5 degree cells'),('all_k3.5_raw_cell2.5','2.5 degree cells'),
                           ('all_k3.5_raw_cell10.0','10 degree cells'),('all_k3.5_raw_shifted','shifted 5 degree cells')]:
            ratio=np.sqrt(np.diag(z[name+'__spatial_bootstrap_covariance'])/np.diag(z[name+'__independent_stack_covariance']))
            axes[1].plot(np.arange(6),ratio,'o-',label=label)
    axes[1].axhline(1,color='k',lw=.7);axes[1].set_xlabel('Radial bin (50-500 kpc)')
    axes[1].set_ylabel('Spatial-bootstrap SE / independent-field SE');axes[1].legend(fontsize=8)
    axes[1].set_title('Control-only estimate: 40 complete fields')
    fig.suptitle('Coadd injection recovery and control covariance\nPaired cluster-control covariance: zero complete six-bin pairs',fontsize=11)
    fig.tight_layout()
    path=a.analysis/'prestack_validation_review.png'
    with path.open('xb') as stream:fig.savefig(stream,format='png',dpi=150)
    plt.close(fig)
    inputs=[a.analysis/'recovery_groups.tsv',a.analysis/'covariances.npz',Path(__file__)]
    with (a.analysis/'review_plot_provenance.json').open('x') as stream:
        json.dump(dict(source_sha256={str(p):sha(p) for p in inputs},output_sha256=sha(path)),stream,indent=2)


if __name__=='__main__':main()
