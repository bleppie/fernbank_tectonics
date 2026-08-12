#!/usr/bin/env python3
"""Detailed island terrain heightfield: macro shape from the contour guideline,
procedural ridged/fBm detail, and a hydraulic-erosion pass to carve valleys.
Units: meters. +Y up, sea level 0, peak +1500, ocean-floor rim -1500."""
import json, numpy as np
from matplotlib.path import Path
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates

PEAK_Y, FLOOR_Y, SEA_Y = 1500.0, -1500.0, 0.0
LAND_LEVELS = [0.0, 500.0, 1000.0]          # pink / red / brightred contour edges
TARGET_MAX_EXTENT_M = 20000.0
N = 620                                       # raster resolution (larger axis)
SEED = 7

# ---------- contours -> model space ----------
def load_contours():
    C = json.load(open("contours.json"))
    dom = np.array(C["domain"])
    cx, cy = dom[:,0].mean(), dom[:,1].mean()
    x0,y0,x1,y1 = dom[:,0].min(),dom[:,1].min(),dom[:,0].max(),dom[:,1].max()
    scale = TARGET_MAX_EXTENT_M / max(x1-x0, y1-y0)
    f = lambda pl: np.column_stack([(np.asarray(pl)[:,0]-cx)*scale,
                                    (cy-np.asarray(pl)[:,1])*scale])
    return {k:f(v) for k,v in C.items()}, scale

# ---------- perlin noise (vectorised, float-coord sampling) ----------
def _perm(seed):
    rng = np.random.default_rng(seed)
    p = np.arange(256); rng.shuffle(p); return np.concatenate([p,p]).astype(int)
_G = np.array([[1,1],[-1,1],[1,-1],[-1,-1],[1,0],[-1,0],[0,1],[0,-1]],float)
def perlin(x, y, perm):
    xi = np.floor(x).astype(int) & 255; yi = np.floor(y).astype(int) & 255
    xf = x-np.floor(x); yf = y-np.floor(y)
    fade = lambda t: t*t*t*(t*(t*6-15)+10)
    u,v = fade(xf), fade(yf)
    def grad(h, X, Y):
        g = _G[perm[h] & 7]; return g[...,0]*X + g[...,1]*Y
    aa=perm[perm[xi]+yi]; ab=perm[perm[xi]+yi+1]
    ba=perm[perm[xi+1]+yi]; bb=perm[perm[xi+1]+yi+1]
    x1 = grad(aa,xf,yf)*(1-u)+grad(ba,xf-1,yf)*u
    x2 = grad(ab,xf,yf-1)*(1-u)+grad(bb,xf-1,yf-1)*u
    return x1*(1-v)+x2*v                       # ~[-1,1]

def fbm(x,y,perm,oct=7,lac=2.0,gain=0.5,ridged=False):
    amp=1.0; freq=1.0; s=np.zeros_like(x); norm=0.0
    for _ in range(oct):
        n=perlin(x*freq,y*freq,perm)
        if ridged: n=(1.0-np.abs(n)); n=n*n
        s+=amp*n; norm+=amp; amp*=gain; freq*=lac
    s/=norm
    return s

# ---------- hydraulic erosion (vectorised droplets) ----------
def erode(h, n_drops, seed, steps=36, inertia=0.05, cap=3.6, dep=0.09,
          ero=0.30, grav=4.0, evap=0.02, radius=1.6, min_slope=0.0006, land=None):
    H,W = h.shape
    h = h.astype(np.float64)
    lo, hi = h.min(), h.max(); span = max(hi-lo, 1e-6)
    h = (h.copy()-lo)/span                     # work in normalised [0,1]
    rng = np.random.default_rng(seed)
    # brush for spreading erosion
    r=int(np.ceil(radius)); yy,xx=np.mgrid[-r:r+1,-r:r+1]
    wk=np.maximum(0,radius-np.hypot(xx,yy)); wk/=wk.sum(); bh,bw=wk.shape
    def sample(hh, px, py):
        c=np.stack([py,px]); return map_coordinates(hh,c,order=1,mode="nearest")
    for b in range(0, n_drops, 40000):
        m=min(40000, n_drops-b)
        if land is not None:                    # spawn on land for valley cutting
            ys,xs=np.where(land)
            idx=rng.integers(0,len(xs),m); px=xs[idx]+rng.random(m); py=ys[idx]+rng.random(m)
        else:
            px=rng.random(m)*(W-1); py=rng.random(m)*(H-1)
        dx=np.zeros(m); dy=np.zeros(m); vel=np.ones(m); water=np.ones(m); sed=np.zeros(m)
        for _ in range(steps):
            gx0=sample(h,px+1,py)-sample(h,px-1,py)
            gy0=sample(h,px,py+1)-sample(h,px,py-1)
            gx=gx0*0.5; gy=gy0*0.5
            dx=dx*inertia-gx*(1-inertia); dy=dy*inertia-gy*(1-inertia)
            L=np.hypot(dx,dy)+1e-9; dx/=L; dy/=L
            npx=px+dx; npy=py+dy
            inb=(npx>=1)&(npx<W-2)&(npy>=1)&(npy<H-2)
            hold=sample(h,px,py); hnew=sample(h,npx,npy)
            dh=hnew-hold
            c=np.maximum(-dh,min_slope)*vel*water*cap
            over=sed-c
            depos=np.where(dh>0, np.minimum(sed,dh), np.where(over>0, over*dep, 0.0))
            erod=np.where((dh<=0)&(over<0), np.minimum(-over*ero,-dh*0.9), 0.0)
            # deposit at current cell (bilinear)
            ix=np.floor(px).astype(int); iy=np.floor(py).astype(int)
            fx=px-ix; fy=py-iy
            for (ox,oy,wgt) in ((0,0,(1-fx)*(1-fy)),(1,0,fx*(1-fy)),(0,1,(1-fx)*fy),(1,1,fx*fy)):
                np.add.at(h,(np.clip(iy+oy,0,H-1),np.clip(ix+ox,0,W-1)),depos*wgt)
            # erode spread over brush
            eff=erod.copy()
            for j in range(bh):
                for i in range(bw):
                    yy2=np.clip(iy+ j-r,0,H-1); xx2=np.clip(ix+ i-r,0,W-1)
                    np.add.at(h,(yy2,xx2),-eff*wk[j,i])
            sed=sed-depos+erod
            vel=np.sqrt(np.maximum(0.0, vel*vel + dh*(-grav)))
            water*=(1-evap)
            px,py=npx,npy
            px=np.clip(px,1,W-2); py=np.clip(py,1,H-2)
    return h*span+lo                            # back to meters

# ---------- assemble field ----------
def build_field(verbose=True):
    C,scale = load_contours()
    dom=C["domain"]; pink=C["pink"]; red=C["red"]; peak=C["peak"]
    Xmin,Xmax=dom[:,0].min(),dom[:,0].max(); Zmin,Zmax=dom[:,1].min(),dom[:,1].max()
    asp=(Zmax-Zmin)/(Xmax-Xmin)
    if asp>=1: nz,nx=N,max(2,int(round(N/asp)))
    else: nx,nz=N,max(2,int(round(N*asp)))
    xs=np.linspace(Xmin,Xmax,nx); zs=np.linspace(Zmin,Zmax,nz)
    XX,ZZ=np.meshgrid(xs,zs)
    def mask(poly):
        p=Path(poly); pts=np.column_stack([XX.ravel(),ZZ.ravel()])
        return p.contains_points(pts).reshape(nz,nx)
    m_dom=mask(dom); m_pink=mask(pink)&m_dom; m_red=mask(red)&m_dom; m_peak=mask(peak)&m_dom
    Din=lambda m: distance_transform_edt(m); Dout=lambda m: distance_transform_edt(~m)
    E=np.full((nz,nx),FLOOR_Y,float)
    ocean=m_dom&~m_pink
    dc=Dout(m_pink); de=Din(m_dom); t=dc/(dc+de+1e-9)
    E[ocean]=SEA_Y+(t[ocean]**1.3)*(FLOOR_Y-SEA_Y)        # gentle then steeper -> rim
    r1=m_pink&~m_red; a=Din(m_pink); bb=Dout(m_red); t=a/(a+bb+1e-9)
    E[r1]=LAND_LEVELS[0]+t[r1]*(LAND_LEVELS[1]-LAND_LEVELS[0])
    r2=m_red&~m_peak; a=Din(m_red); bb=Dout(m_peak); t=a/(a+bb+1e-9)
    E[r2]=LAND_LEVELS[1]+t[r2]*(LAND_LEVELS[2]-LAND_LEVELS[1])
    dp=Din(m_peak); E[m_peak]=LAND_LEVELS[2]+(dp[m_peak]/max(dp.max(),1))*(PEAK_Y-LAND_LEVELS[2])
    macro=gaussian_filter(E,2.0)

    # ---- procedural detail ----
    perm=_perm(SEED)
    freq=nx/6.0                                  # base spatial frequency
    u=XX/(Xmax-Xmin); v=ZZ/(Zmax-Zmin)
    wx=fbm(u*3+11,v*3+4,perm,oct=4); wy=fbm(u*3+7,v*3+19,perm,oct=4)  # domain warp
    U=(u+0.35*wx)*freq; V=(v+0.35*wy)*freq
    ridg=fbm(U,V,perm,oct=8,ridged=True)         # 0..1 ridge lines
    roll=fbm(U*0.6+3,V*0.6+3,perm,oct=6)         # -1..1 rolling hills
    detail=(ridg-0.5)*2.0*0.62 + roll*0.38       # ~[-1,1]

    land=m_pink
    # relief amplitude: strong on land (esp. upper slopes), gentle undersea
    landdist=Din(m_pink)/ (Din(m_pink).max()+1e-9)
    amp_land=430.0*(0.35+0.65*np.sqrt(landdist))
    amp=np.where(land, amp_land, 70.0)           # gentle submarine texture
    field=macro+detail*amp

    # keep peak high, coastline near 0, rim low
    field=np.where(m_dom, field, FLOOR_Y)
    if verbose: print(f"grid {nx}x{nz}; pre-erosion range {field.min():.0f}..{field.max():.0f}")

    # ---- erosion (land-spawned droplets carve dendritic valleys) ----
    ndrop=int(nx*nz*0.9)
    field=erode(field, ndrop, SEED+1, land=land)
    field=gaussian_filter(field,0.6)             # de-speckle

    # ---- normalise exact bounds, keep sea level at 0 ----
    field=np.where(m_dom, field, FLOOR_Y)
    pos=(field>0)&m_dom; neg=(field<0)&m_dom
    if pos.any(): field[pos]*=PEAK_Y/field[pos].max()
    if neg.any(): field[neg]*=FLOOR_Y/field[neg].min()
    field=np.clip(field,FLOOR_Y,PEAK_Y)
    # force domain rim to the lowest elevation
    rim=m_dom&(~gaussian_filter(m_dom.astype(float),1.5).astype(bool)==False)
    if verbose: print(f"post range {field.min():.0f}..{field.max():.0f}")
    return dict(E=field,xs=xs,zs=zs,nx=nx,nz=nz,m_dom=m_dom,m_pink=m_pink,
                Xmin=Xmin,Xmax=Xmax,Zmin=Zmin,Zmax=Zmax,dom=dom,scale=scale)

if __name__=="__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    d=build_field(); E=d["E"]
    fig=plt.figure(figsize=(18,6))
    ax=fig.add_subplot(1,3,1); im=ax.imshow(E,origin="lower",cmap="terrain",vmin=-1500,vmax=1500)
    ax.set_title("elevation"); plt.colorbar(im,ax=ax,fraction=.046)
    ax=fig.add_subplot(1,3,2); ls=LightSource(315,45)
    ax.imshow(ls.shade(E,cmap=plt.cm.gist_earth,vert_exag=2.2,blend_mode="soft",vmin=-1500,vmax=1500),origin="lower")
    ax.set_title("hillshade (detail check)"); ax.axis("off")
    ax=fig.add_subplot(1,3,3); ax.plot(E[E.shape[0]//2,:]); ax.axhline(0,color="b",ls="--",lw=.7)
    ax.axhline(1500,color="r",ls=":",lw=.6); ax.axhline(-1500,color="brown",ls=":",lw=.6)
    ax.set_title("W-E cross-section"); ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig("verify_field.png",dpi=105); print("saved verify_field.png")
