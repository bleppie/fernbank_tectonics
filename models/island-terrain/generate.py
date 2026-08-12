#!/usr/bin/env python3
"""
Detailed island-terrain 3D model -> OBJ.

Macro shape comes from the simplified contour guideline (peak location, coast,
ocean rim); fine relief (rugged ridges + dendritic erosion valleys, at roughly
the level of detail in the reference image) is generated procedurally and
carved with a hydraulic-erosion pass. See terrain.py for the field builder.

Units: METERS. Altitude on +Y. Sea level Y = 0.
Highest mountain = +1500 m. Ocean-floor rim (outer contour = model limit) =
-1500 m, with the sea bed gently sloping from the island out to that rim.

The mesh is clipped to the outer contour (the model limit) and closed into a
watertight solid: terrain surface on top, a vertical skirt around the rim, and
a flat base just below the ocean floor.
"""
import numpy as np, trimesh
from scipy.spatial import Delaunay
from matplotlib.path import Path
import terrain as T

BASE_BELOW = 80.0          # base plane below the ocean floor
EDGE_MARGIN = 0.45         # drop interior nodes this close (in cells) to the rim

def resample_loop(poly, step):
    poly = np.asarray(poly)
    if np.allclose(poly[0], poly[-1]): poly = poly[:-1]
    out = []
    n = len(poly)
    for i in range(n):
        a = poly[i]; b = poly[(i+1) % n]
        L = np.hypot(*(b-a)); k = max(1, int(round(L/step)))
        for t in np.linspace(0, 1, k, endpoint=False):
            out.append(a + (b-a)*t)
    return np.array(out)

def build_mesh(use_cache=False):
    import os, pickle
    cache = "field_cache.pkl"
    if use_cache and os.path.exists(cache):
        d = pickle.load(open(cache, "rb"))
    else:
        d = T.build_field()
        pickle.dump(d, open(cache, "wb"))
    E, xs, zs = d["E"], d["xs"], d["zs"]
    nx, nz = d["nx"], d["nz"]
    dom = d["dom"]
    dx = xs[1]-xs[0]

    XX, ZZ = np.meshgrid(xs, zs)
    inside = d["m_dom"]
    # shrink mask slightly so interior nodes stay clear of the rim loop
    from scipy.ndimage import binary_erosion
    core = binary_erosion(inside, iterations=1)
    ii = np.where(core.ravel())[0]
    ipts = np.column_stack([XX.ravel()[ii], ZZ.ravel()[ii]])
    iy = E.ravel()[ii]

    loop = resample_loop(dom, dx)                 # ordered rim, lowest elevation
    ly = np.full(len(loop), T.FLOOR_Y)

    pts2d = np.vstack([ipts, loop])
    yvals = np.concatenate([iy, ly])
    n_interior = len(ipts)
    loop_ids = np.arange(n_interior, n_interior+len(loop))

    # convex domain: Delaunay of (interior + rim loop) tiles the octagon exactly
    tri = Delaunay(pts2d)
    faces = tri.simplices.tolist()
    top = np.column_stack([pts2d[:,0], yvals, pts2d[:,1]])

    # actual boundary of the top surface = edges used by exactly one triangle
    from collections import defaultdict
    ec = defaultdict(int); en = defaultdict(list)
    for a,b,c in faces:
        for u,v in ((a,b),(b,c),(c,a)):
            k=(min(u,v),max(u,v)); ec[k]+=1
    bedges=[k for k,n in ec.items() if n==1]
    adj=defaultdict(list)
    for u,v in bedges: adj[u].append(v); adj[v].append(u)
    start=bedges[0][0]; loop_order=[start]; prev=None; cur=start
    while True:
        nxts=[w for w in adj[cur] if w!=prev]
        if not nxts: break
        nxt=nxts[0]
        if nxt==start: break
        loop_order.append(nxt); prev,cur=cur,nxt
    # base vertices under each boundary vertex
    base_off=len(top)
    base=top[loop_order].copy(); base[:,1]=T.FLOOR_Y-BASE_BELOW
    bmap={vid:base_off+i for i,vid in enumerate(loop_order)}
    F=[tuple(f) for f in faces]
    P=len(loop_order)
    for i in range(P):                             # skirt
        a=loop_order[i]; b=loop_order[(i+1)%P]
        F.append((a,b,bmap[b])); F.append((a,bmap[b],bmap[a]))
    for k in range(1,P-1):                          # base fan
        F.append((base_off,base_off+k+1,base_off+k))

    verts=np.vstack([top,base]); faces=np.array(F)
    mesh=trimesh.Trimesh(vertices=verts,faces=faces,process=True)
    mesh.fix_normals()
    trimesh.repair.fill_holes(mesh)
    return mesh, d

if __name__ == "__main__":
    mesh, d = build_mesh(use_cache=True)
    y = mesh.vertices[:,1]
    print(f"verts {len(mesh.vertices)} | faces {len(mesh.faces)} | watertight {mesh.is_watertight}")
    print(f"elevation min/max: {y.min():.0f} / {y.max():.0f} m (sea level 0)")
    print(f"extent: {mesh.bounds[1,0]-mesh.bounds[0,0]:.0f} x {mesh.bounds[1,2]-mesh.bounds[0,2]:.0f} m")
    mesh.export("island_terrain.obj")
    mesh.export("island_terrain.glb")
    print("wrote island_terrain.obj / .glb")
