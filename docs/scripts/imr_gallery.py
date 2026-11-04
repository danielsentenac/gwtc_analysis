"""Waveform gallery of the home page of the docs: docs/img/imr_gallery.png.

Eight representative binaries from 20 Hz (BNS, NSBH, equal and unequal BBH from 20 to 238 Msun, aligned and
precessing spins): left, the whole signal in band in seconds; right, the last orbits, merger and ringdown in units of
the total mass M. Durations and cycles come from the dominant mode of aligned-spin models (IMRPhenomXAS,
IMRPhenomNSBH, IMRPhenomD_NRTidalv2); the precessing case is drawn with IMRPhenomXPHM; ringdown frequency and damping
time are the fundamental quasi-normal mode from fits to the remnant mass and spin. Needs pycbc (not a gwtc_analysis
dependency). Run from anywhere: python docs/scripts/imr_gallery.py
(if pycbc fails at import on getconf, set LEVEL2_CACHE_SIZE=1048576 in the environment).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pycbc.waveform import get_td_waveform
from pycbc import conversions as cv
T=4.925491e-6  # G Msun / c^3 in s
cases=[
 ("BNS, GW170817-like","1.46 + 1.27",1.46,1.27,0,0,0,dict(approximant="IMRPhenomPv2_NRTidalv2",lambda1=400,lambda2=400),dict(approximant="IMRPhenomD_NRTidalv2",lambda1=400,lambda2=400)),
 ("NSBH","8 + 1.4",8.0,1.4,0,0,0,dict(approximant="IMRPhenomNSBH",lambda2=400),dict(approximant="IMRPhenomNSBH",lambda2=400)),
 ("Light BBH, equal","10 + 10",10,10,0,0,0,{},{}),
 ("BBH, GW150914-like","36 + 29",36,29,0,0,0,{},{}),
 ("Heavy BBH, GW231123-like","137 + 101",137,101,0,0,0,{},{}),
 ("Unequal BBH, GW190814-like","23 + 2.6",23,2.6,0,0,0,{},{}),
 ("Equal, aligned spins 0.9","30 + 30",30,30,0.9,0.9,0,{},{}),
 ("Unequal, precessing (χ₁ = 0.9 in plane)","30 + 6",30,6,0,0,0.9,{},{}),
]
res=[]
for name,lab,m1,m2,s1z,s2z,s1x,kwplot,kwmet in cases:
    M=m1+m2; eta=m1*m2/M**2; Mc=(m1*m2)**0.6/M**0.2
    dt=1/16384 if M<5 else (1/8192 if M<15 else 1/4096)
    # metrics: dominant mode, aligned
    km=dict(approximant="IMRPhenomXAS",spin1z=s1z,spin2z=s2z); km.update(kwmet)
    hp,hc=get_td_waveform(mass1=m1,mass2=m2,delta_t=dt,f_lower=20.0,distance=400,inclination=0,**km)
    t=hp.sample_times.numpy(); h=hp.numpy()-1j*hc.numpy(); a=np.abs(h); ip=int(np.argmax(a))
    ph=np.unwrap(np.angle(h)); f=np.abs(np.gradient(ph,dt))/(2*np.pi)
    good=a>1e-3*a.max(); k=max(1,int(0.01/dt)); fs=np.convolve(f,np.ones(k)/k,'same')
    idx=np.where((fs>=20)&good)[0]; i20=idx[0] if len(idx) else 0
    if i20>=ip: i20=0
    dur=t[ip]-t[i20]; cyc=abs(ph[ip]-ph[i20])/(2*np.pi)
    fpk=fs[ip]
    fisco_=1/(6**1.5*np.pi*M*T); ii=np.where((fs>=fisco_)&good)[0]; tisco=(t[ip]-t[ii[0]]) if len(ii) and ii[0]<ip else 0.0
    fisco=1/(6**1.5*np.pi*M*T)
    if not name.startswith("BNS"):
        Mf=cv.final_mass_from_initial(m1,m2,spin1z=s1z,spin2z=s2z,spin1x=s1x); chif=cv.final_spin_from_initial(m1,m2,spin1z=s1z,spin2z=s2z,spin1x=s1x)
        Mf=float(Mf); chif=float(chif)
        fqnm=(1.5251-1.1568*(1-chif)**0.1292)/(2*np.pi*Mf*T); Q=0.7+1.4187*(1-chif)**-0.4990; tau=Q/(np.pi*fqnm)
    else:
        Mf=chif=fqnm=tau=float('nan')
    # plot waveform: full model, face-on
    kp=dict(approximant=("IMRPhenomXPHM" if s1x else "IMRPhenomXAS"),spin1z=s1z,spin2z=s2z,spin1x=s1x); kp.update(kwplot)
    hp2,hc2=get_td_waveform(mass1=m1,mass2=m2,delta_t=dt,f_lower=20.0,distance=400,inclination=0,**kp)
    t2=hp2.sample_times.numpy(); y=hp2.numpy(); a2=np.abs(y-1j*hc2.numpy()); ip2=int(np.argmax(a2))
    t2=t2-t2[ip2]
    sel=t2>=-(dur+0.5)
    t2=t2[sel]; y=y[sel]
    late=t2>-1000*M*T
    step=max(1,int((~late).sum())//150000)
    tt=np.concatenate([t2[~late][::step],t2[late]]); yy=np.concatenate([y[~late][::step],y[late]])
    res.append(dict(name=name,lab=lab,M=M,q=m2/m1,Mc=Mc,dur=dur,cyc=cyc,fisco=fisco,tisco=tisco,fpk=fpk,Mf=Mf,chif=chif,fqnm=fqnm,tau=tau,
                    t=tt.tolist(),h=yy.tolist()))
    print(f"{name:42s} M={M:6.1f} q={m2/m1:.2f} Mc={Mc:5.2f} dur={dur:8.3f}s cyc={cyc:6.0f} fISCO={fisco:5.0f} fpk={fpk:5.0f} Mf={Mf:6.1f} chif={chif:.2f} fRD={fqnm:5.0f} tau={tau*1e3:5.1f}ms",flush=True)

T=4.925491e-6
r=res
NAVY,ORANGE,GREEN,GREY="#1f3a5f","#d9822b","#2e7d5b","#666"
fig,axes=plt.subplots(len(r),2,figsize=(13,15),gridspec_kw=dict(width_ratios=[1.25,1]))
for i,x in enumerate(r):
    t=np.array(x['t']); h=np.array(x['h']); hn=h/np.abs(h).max()
    M=x['M']; ti=min(x['tisco'],x['dur'])
    ax=axes[i,0]
    ax.plot(t,hn,color=NAVY,lw=0.35)
    ax.set_xlim(-x['dur']*1.02, max(0.02*x['dur'],60*M*T))
    ax.axvspan(-ti,0,color=ORANGE,alpha=0.25,lw=0); ax.axvspan(0,ax.get_xlim()[1],color=GREEN,alpha=0.18,lw=0)
    ax.set_yticks([]); ax.set_ylim(-1.15,1.15)
    ax.text(0.01,0.93,f"{x['name']}  ({x['lab']} M$_\\odot$)",transform=ax.transAxes,fontsize=9.5,fontweight='bold',color=NAVY,va='top')
    ax.text(0.01,0.08,f"{x['dur']:.3g} s and {x['cyc']:.0f} cycles from 20 Hz",transform=ax.transAxes,fontsize=8.5,color=GREY)
    ax.tick_params(labelsize=8)
    if i==len(r)-1: ax.set_xlabel("time to merger (s), from 20 Hz",fontsize=9)
    ax=axes[i,1]
    tm=t/(M*T); sel=(tm>-400)&(tm<150)
    ax.plot(tm[sel],hn[sel],color=NAVY,lw=0.9)
    ax.axvspan(-ti/(M*T),0,color=ORANGE,alpha=0.25,lw=0); ax.axvspan(0,150,color=GREEN,alpha=0.18,lw=0)
    ax.set_xlim(-400,150); ax.set_yticks([]); ax.set_ylim(-1.15,1.15); ax.tick_params(labelsize=8)
    rd = "post-merger (not modelled)" if x['name'].startswith('BNS') else f"ringdown {x['fqnm']:.0f} Hz, τ = {x['tau']*1e3:.2g} ms"
    ax.text(0.99,0.93,f"peak {x['fpk']:.0f} Hz · {rd}",transform=ax.transAxes,fontsize=8.5,color=GREY,ha='right',va='top')
    if i==len(r)-1: ax.set_xlabel("time from peak (units of total mass M, GM/c³)",fontsize=9)
axes[0,0].set_title("Whole signal in band (h₊, face-on, normalised)",fontsize=10.5,color=NAVY,fontweight='bold')
axes[0,1].set_title("Last orbits, merger, ringdown (same scale in M)",fontsize=10.5,color=NAVY,fontweight='bold')
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=ORANGE,alpha=0.4,label="from ISCO frequency to peak (plunge, merger)"),Patch(color=GREEN,alpha=0.3,label="after peak (ringdown)")],
           loc='lower center',ncol=2,fontsize=9,frameon=False,bbox_to_anchor=(0.5,0.0))
fig.text(0.5,0.018,"IMRPhenomXAS (BBH, dominant mode), IMRPhenomXPHM (precessing), IMRPhenomNSBH, IMRPhenomPv2_NRTidalv2 (BNS, Λ = 400); pycbc, f_low = 20 Hz; ringdown f and τ: fundamental QNM fits",ha='center',fontsize=7.5,color=GREY,style='italic')
fig.tight_layout(rect=(0,0.035,1,1))
fig.savefig(Path(__file__).resolve().parents[1] / "img" / "imr_gallery.png", dpi=130)
