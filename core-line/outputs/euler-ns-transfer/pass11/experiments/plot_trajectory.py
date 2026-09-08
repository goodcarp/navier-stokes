#!/usr/bin/env python3
"""Plot recorded finite-cylinder diagnostics, never inferred continuum data."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/forced-route-mpl')
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(paths,output):
    runs=[json.loads(Path(p).read_text()) for p in paths]
    if not all(r.get('complete') is True for r in runs):raise ValueError('Only completed trajectories may be plotted')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.8));fig.subplots_adjust(top=.72,bottom=.23,left=.075,right=.98,wspace=.36)
    fig.text(.075,.92,'A finite pulse, with a moving receiver',fontsize=21,weight='bold',color='#172a3a')
    fig.text(.075,.845,'Unchanged full 3D datum · viscosity 0.001 · finite no-slip cylinder · periodic axial direction',fontsize=10.5,color='#536471')
    colors=['#bd6b3b','#176980'];styles=['--','-']
    for index,run in enumerate(runs):
        h=run['history'];t=np.array([x['time'] for x in h])*1000
        c=colors[min(index,1)];style=styles[min(index,1)];label='dt = '+format(run['parameters']['dt'],'.0e')
        K=np.array([x['K'] for x in h]);G=np.array([x['mean_G_branch'] for x in h]);Gf=np.array([x['mean_G_fixed'] for x in h])
        cone=np.array([x['core_b']-x['core_Omega'] for x in h])
        axes[0].plot(t,100*(K/K[0]-1),style,color=c,lw=2,label=label)
        axes[1].plot(t,G-h[0]['mean_G_branch'],style,color=c,lw=2,label='Selected receiver, '+label)
        if index==len(runs)-1:axes[1].plot(t,Gf-h[0]['mean_G_fixed'],':',color='#8b5360',lw=2,label='Original receiver')
        axes[2].plot(t,cone*1e5,style,color=c,lw=2,label=label)
    titles=['Fluctuation energy','Mean circulation','Core strain minus rotation']
    labels=['Change in K (%)','Change in G','(b − Ω) × 10⁵']
    for ax,title,label in zip(axes,titles,labels):
        ax.set_title(title,loc='left',fontweight='bold',pad=13);ax.set_xlabel('Time × 10³');ax.set_ylabel(label)
        ax.axhline(0,color='#adb6bf',lw=.8);ax.grid(axis='y',alpha=.18);ax.set_xlim(0,1)
    axes[0].legend(frameon=False,fontsize=8.5)
    axes[1].legend(frameon=False,fontsize=7.5,loc='lower left')
    fig.text(.075,.065,'Numerical diagnostics only. Spatial, domain, reconstruction and residual errors remain unbounded; no NS stage or blowup is certified.',fontsize=9,color='#536471')
    fig.savefig(output,dpi=180,facecolor='white');plt.close(fig)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('runs',nargs='+');ap.add_argument('--output',required=True)
    args=ap.parse_args();plot(args.runs,args.output)
