#!/usr/bin/env python3
"""Scientific plot of the stored planar diagnostic; no simulation rerun."""
import argparse,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path('/tmp')/'codex-ns-pass9-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def run(folder,out):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(11.8,4.7))
    fig.subplots_adjust(left=.08,right=.97,bottom=.16,top=.76,wspace=.3)
    fig.suptitle('Finite gains in the planar model',x=.08,ha='left',fontsize=18,fontweight='bold',y=.96)
    fig.text(.08,.875,'Nonlinear slice with prescribed strain; compact 3D evolution remains unverified.',fontsize=10,color='#4b5563')
    for name,color,label in [('slice-1536.json','#8b5cf6','A = 950'),('slice-fixed-A1020.json','#0369a1','A = 1020, selected datum')]:
        data=json.loads((folder/name).read_text())['records'];t=np.array([x['t'] for x in data])
        G=np.array([x['mean_maximum'] for x in data]);fixed=np.array([x['mean_fixed_initial_receiver'] for x in data]);K=np.array([x['fluctuation_energy_per_axial_length'] for x in data])
        axs[0].plot(t*1000,100*(G/G[0]-1),color=color,lw=2.4,label=label)
        axs[1].plot(t*1000,100*(K/K[0]-1),color=color,lw=2.4,label=label)
        if name=='slice-1536.json':
            axs[0].plot(t*1000,100*(fixed/fixed[0]-1),color=color,lw=1.3,ls=':',label='A = 950, fixed receiver')
    axs[0].set_title('Receiving mean angular momentum',loc='left',fontsize=11,pad=12)
    axs[1].set_title('Fluctuation energy within r < 0.8',loc='left',fontsize=11,pad=12)
    for ax in axs:
        ax.axhline(0,color='#6b7280',lw=.8)
        ax.axvline(1,color='#9ca3af',lw=.8,ls='--')
        ax.grid(axis='y',color='#e5e7eb',lw=.7)
        ax.set_xlabel('Time × 1000 (model units)')
        ax.set_ylabel('Change from initial value (%)')
        ax.legend(frameon=False,fontsize=8,loc='best')
    fig.text(.08,.035,'Mean curve follows the sampled radial critical branch. Energy is per unit axial length; exterior tails are excluded.',fontsize=8,color='#4b5563')
    for ext in ['png','svg']:fig.savefig(out/f'slice-pulse.{ext}',dpi=160,facecolor='white')
    plt.close(fig)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--folder',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();run(args.folder,args.output)
