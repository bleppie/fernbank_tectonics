#!/usr/bin/env python3
"""
Flat 8-sided ocean plane matching the island model.

A single flat polygon (the octagon outline that surrounds the island), with:
  * NO thickness (single-sided fan, 8 verts / 6 tris)
  * normals pointing straight up (+Y), baked as explicit vertex normals
  * y = 0 (sea level)  -- same datum as the island (sea level = 0)
  * the EXACT same horizontal transform as the island model, including the
    X-mirror (MIRROR_X) applied there, so the outline lines up perfectly.

Outputs ocean_plane.obj (meters). Change PLANE_Y if you want it at a different
height than sea level.
"""
import numpy as np, trimesh, terrain as T

PLANE_Y = 0.0   # sea level (same datum as the island)

C, scale = T.load_contours()          # center + (largest extent -> 20 km) + y-flip
dom = C["domain"].copy()              # (8,2) model coords: X, Z
dom[:, 0] *= -1.0                     # MIRROR_X, matching the island model

verts = np.column_stack([dom[:, 0], np.full(len(dom), PLANE_Y), dom[:, 1]])
faces = np.array([[0, i, i + 1] for i in range(1, len(dom) - 1)])   # convex fan
m = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
if m.face_normals[:, 1].mean() < 0:                                 # ensure +Y
    m = trimesh.Trimesh(vertices=verts, faces=m.faces[:, ::-1], process=False)

obj = trimesh.exchange.obj.export_obj(m, include_normals=True, include_texture=False)
open("ocean_plane.obj", "w").write(obj)
print(f"verts {len(m.vertices)} | faces {len(m.faces)} | Y={PLANE_Y}")
print("all face normals +Y:", bool(np.allclose(m.face_normals, [0, 1, 0])))
print("extent X %.0f  Z %.0f m (matches island footprint)"
      % (m.bounds[1,0]-m.bounds[0,0], m.bounds[1,2]-m.bounds[0,2]))
print("wrote ocean_plane.obj")
