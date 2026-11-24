#!/usr/bin/env python3
"""Render verified member coverage sensitivity tables as a standalone figure."""
import argparse
import json
import os
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze_s23b_member_outer_icl_relation import read_tsv, sha256

MODES = [('native','Each cut','#285fa6',-.022),
         ('common_clusters','Same clusters','#c36a1b',0.),
         ('common_candidates','Same candidates','#387b53',.022)]


def plot(analysis, output):
    if output.exists():
        raise FileExistsError(f'refusing to overwrite {output}')
    summary = json.loads((analysis/'summary.json').read_text())
    if summary.get('status') != 'complete_exploratory_fixed_support_sensitivity' or summary.get('science_measurement') is not False:
        raise ValueError('completed exploratory sensitivity products required')
    for name,digest in summary['output_sha256'].items():
        if Path(name).name != name or sha256(analysis/name) != digest:
            raise ValueError('sensitivity products changed')
    fits = read_tsv(analysis/'coverage_slope_layouts.tsv')
    primary = [r for r in fits if float(r['block_width_deg'])==5. and float(r['grid_phase_fraction'])==0.]
    if len(primary)!=18:
        raise ValueError('expected 18 primary threshold/mode/filter records')
    fig,axs=plt.subplots(2,2,figsize=(10.8,8.2),sharex=True,sharey='row')
    max_n=max(int(r['n_target_clusters']) for r in primary)
    complete=[r for r in primary if r['status']=='complete']
    extrema=[float(r[field]) for r in complete for field in ('slope','ci_low','ci_high')]
    low=min([0.]+extrema); high=max([0.]+extrema)
    span=max(high-low,.2)
    for col,band in enumerate(('i','i2')):
        axn,axslope=axs[:,col]
        for mode,label,color,offset in MODES:
            rows=sorted((r for r in primary if r['filter']==band and r['mode']==mode),key=lambda r:float(r['coverage_floor']))
            axn.plot([float(r['coverage_floor'])+offset for r in rows],
                     [int(r['n_target_clusters']) for r in rows],'o-',color=color,label=label,ms=4)
            for r in rows:
                if r['status']!='complete':
                    continue
                x=float(r['coverage_floor'])+offset
                axslope.vlines(x,float(r['ci_low']),float(r['ci_high']),color=color,lw=1.3)
                axslope.plot(x,float(r['slope']),'o',color=color,ms=4)
        axn.set_ylim(0,max(max_n,1)*1.65)
        axslope.set_ylim(low-.12*span,high+.5*span)
        axn.text(.04,.95,f'HSC-{band}: supported clusters',transform=axn.transAxes,va='top',fontsize=10)
        axslope.text(.04,.95,f'HSC-{band}: adjusted slope',transform=axslope.transAxes,va='top',fontsize=10)
        if not any(r['filter']==band for r in complete):
            axslope.text(.5,.5,'Insufficient sample / spatial support',ha='center',transform=axslope.transAxes,fontsize=10)
        axslope.axhline(0,color='.6',lw=.7,zorder=0)
        axslope.set_xlabel('Target and control coverage floor')
        for ax in (axn,axslope):
            ax.set_xlim(.32,.88)
            ax.set_xticks([.4,.6,.8])
            ax.tick_params(direction='in',top=True,right=True)
    axs[0,0].set_ylabel('Clusters with outer profiles')
    axs[1,0].set_ylabel('Outer-light / relative ICL slope')
    axs[0,1].legend(loc='upper right',fontsize=10,frameon=False)
    fig.subplots_adjust(wspace=0,hspace=.08,left=.1,right=.98,top=.98,bottom=.1)
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix=output.name+'.partial.',dir=output.parent))
    fig.savefig(tmp/'coverage_sensitivity.png',dpi=170)
    fig.savefig(tmp/'coverage_sensitivity.pdf')
    plt.close(fig)
    provenance=dict(status='complete_exploratory_visualization',science_measurement=False,
        source_summary_sha256=sha256(analysis/'summary.json'),script_sha256=sha256(Path(__file__)),
        output_sha256={p.name:sha256(p) for p in tmp.iterdir()},
        interpretation='5-degree phase-zero layout only; marginal spatial-block intervals, no detection or multiplicity claim.')
    (tmp/'summary.json').write_text(json.dumps(provenance,indent=2)+'\n')
    os.rename(tmp,output)
    return provenance


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(plot(a.analysis.resolve(),a.output.resolve()),indent=2))
