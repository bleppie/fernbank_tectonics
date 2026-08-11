#!/usr/bin/env python3
"""
Plate Tectonics Counter -- 8-sided museum exhibit table.
Generates a simple watertight 3D solid from dimensions read off drawing
FM-EX-315 (REV 1, 13 JULY 26). All units are INCHES.

Geometry was extracted from the vector plan view (scale: 4.5 pt = 1 in,
verified against all eight labeled outer edges). In plan there are three
octagon outlines:
    A = outer outline        (table top surface at its lowest, z = 0)
    B = first inner line      (top of the rising border, z = RISE)
    C = second inner line     (flat lip ends here; edge of the 2" drop)
The center of the table (inside C) is the recessed flat area.

Cross-section (per the plan notes), moving inward from the outer edge:
    outer outline  -> first inner line : surface RISES 7.75" (gentle ramp)
    first inner    -> second inner     : FLAT
    second inner   -> central area     : straight vertical DROP of 2"
So heights are: A=0, B=7.75, C=7.75, central floor = 7.75-2 = 5.75.

Four 24"-class touchscreens (~20.62 x 10.22") are embedded (cut through)
the rising border, centered on alternating edges 1,3,5,7.

BASE_THICKNESS is the only value NOT on the drawing: it sets how thick the
slab is below the outer rim so the result is a solid (not a zero-thickness
sheet). Change it freely; it does not affect the top surface geometry.
"""
import json, math, numpy as np, trimesh

# ---- parameters from the drawing (inches) ----
RISE = 7.75            # outer outline -> first inner line
DROP = 2.0             # second inner line -> central flat area
CENTRAL = RISE - DROP  # = 5.75
BASE_THICKNESS = 4.0   # ASSUMED slab thickness below the outer rim
SCREEN_RECESS_DEPTH = 2.0  # ASSUMED pocket depth (perpendicular to the sloped
                           # border) for the embedded touchscreens

# ---- extracted plan geometry (inches, PDF coords; y flipped so z-up is right-handed) ----
A = [[84.0,123.111],[119.556,117.333],[165.333,117.333],[199.111,149.778],
     [196.889,205.778],[136.889,269.333],[88.444,251.111],[55.556,204.889]]
B = [[94.667,136.0],[120.444,131.556],[159.556,131.556],[184.889,155.556],
     [182.667,200.0],[132.889,252.889],[97.333,239.556],[71.111,202.667]]
C = [[96.0,137.778],[120.444,133.778],[159.111,133.778],[182.667,156.0],
     [180.889,199.111],[132.444,250.667],[98.667,237.778],[73.333,202.222]]
# Screen bezel openings (the black frame in the plan), centered on edges 1,3,5,7.
# `edge` = index of the outer edge (between corner[i] and corner[i+1]) the
# screen sits on; the recess is cut perpendicular to that border facet.
SCREENS = [
    {"edge":1,"center":[142.44,124.67],"length":20.89,"width":10.22},
    {"edge":3,"center":[190.89,177.34],"length":20.96,"width":10.29},
    {"edge":5,"center":[108.83,251.49],"length":20.85,"width":10.56},
    {"edge":7,"center":[76.28,166.47], "length":20.93,"width":10.26},
]

def flipy(p):  # PDF y grows downward; negate for a natural z-up model
    return [p[0], -p[1]]

A = [flipy(p) for p in A]; B = [flipy(p) for p in B]; C = [flipy(p) for p in C]

def centroid(poly):
    return [sum(p[0] for p in poly)/len(poly), sum(p[1] for p in poly)/len(poly)]

def build():
    V = []; F = []
    def add(x, y, z):
        V.append((x, y, z)); return len(V) - 1
    def quad(a, b, c, d):           # two tris, CCW
        F.append((a, b, c)); F.append((a, c, d))
    def fan(idx, cen_idx, flip=False):
        n = len(idx)
        for i in range(n):
            a, b = idx[i], idx[(i+1) % n]
            F.append((cen_idx, b, a) if flip else (cen_idx, a, b))

    Ai = [add(x, y, 0.0)     for x, y in A]                 # outer rim (top)
    Bi = [add(x, y, RISE)    for x, y in B]                 # top of rise
    Ct = [add(x, y, RISE)    for x, y in C]                 # lip inner edge (top of drop)
    Cb = [add(x, y, CENTRAL) for x, y in C]                 # bottom of drop (central rim)
    Ab = [add(x, y, -BASE_THICKNESS) for x, y in A]         # outer rim (bottom)
    n = 8
    for i in range(n):
        j = (i+1) % n
        quad(Ai[i], Ai[j], Bi[j], Bi[i])                   # rising border
        quad(Bi[i], Bi[j], Ct[j], Ct[i])                   # flat lip
        quad(Ct[i], Ct[j], Cb[j], Cb[i])                   # vertical drop wall
        quad(Ab[i], Ab[j], Ai[j], Ai[i])                   # outer vertical face
    cc = centroid(C)
    cen_top = add(cc[0], cc[1], CENTRAL)                    # central floor
    fan(Cb, cen_top)
    cen_bot = add(cc[0], cc[1], -BASE_THICKNESS)           # underside
    fan(Ab, cen_bot, flip=True)

    mesh = trimesh.Trimesh(vertices=np.array(V), faces=np.array(F), process=True)
    mesh.fix_normals()

    # Recess the four screens into the sloped border. Each recess is a blind
    # pocket cut PERPENDICULAR to its border facet (so a flat panel mounts
    # flush with the ramp), of uniform depth SCREEN_RECESS_DEPTH.
    Apts = [np.array([x, y, 0.0])  for x, y in A]      # outer corners, z=0
    Bpts = [np.array([x, y, RISE]) for x, y in B]      # inner corners, z=RISE
    over = 1.0                                         # top overshoot above surface
    for s in SCREENS:
        i = s["edge"]; j = (i + 1) % 8
        p0, p1, p2 = Apts[i], Apts[j], Bpts[i]         # three points on the facet
        u = p1 - p0; u /= np.linalg.norm(u)            # tangent (screen length axis)
        n = np.cross(u, p2 - p0); n /= np.linalg.norm(n)
        if n[2] < 0: n = -n                            # surface normal, pointing up
        w = np.cross(n, u)                             # in-plane up-slope (width axis)
        cx, cy = s["center"][0], -s["center"][1]
        # surface point above (cx,cy): intersect vertical line with facet plane
        t = np.dot(p0 - np.array([cx, cy, 0.0]), n) / n[2]
        surf = np.array([cx, cy, t])
        box = trimesh.creation.box(extents=(s["length"], s["width"],
                                            SCREEN_RECESS_DEPTH + over))
        R = np.eye(4); R[:3, 0] = u; R[:3, 1] = w; R[:3, 2] = n
        R[:3, 3] = surf + n * (over - SCREEN_RECESS_DEPTH) / 2.0
        box.apply_transform(R)
        mesh = mesh.difference(box)

    mesh.merge_vertices()
    return mesh

if __name__ == "__main__":
    m = build()
    print("watertight:", m.is_watertight, "| volume(in^3):", round(m.volume, 1),
          "| tris:", len(m.faces))
    ext = m.bounds[1] - m.bounds[0]
    print("bounding box (in):", [round(v, 2) for v in ext],
          "= %.1f' x %.1f'" % (ext[0]/12, ext[1]/12))
    for ext_name in ("obj", "stl", "glb"):
        m.export("plate_tectonics_counter." + ext_name)
    print("wrote plate_tectonics_counter.obj / .stl / .glb")
