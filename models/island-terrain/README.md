# Island Terrain — Detailed 3D Model

A detailed 3D island terrain (rugged mountains, dendritic erosion valleys, and
a gently sloping ocean floor), delivered as **OBJ**.

**Units: meters. Altitude on +Y. Sea level Y = 0.**
- Highest mountain: **+1500 m**
- Ocean-floor rim (outer contour = model limit): **−1500 m**
- Sea bed gently slopes from the island out to that rim.

## How it's built
The contour guideline PDF is **simplified** — it only marks the broad elevation
zones, so it is used as a *macro* guide (where the massif/peak sits, where the
coast is, the model boundary), **not** as exact contours. The fine relief is
generated to match the **level of detail** in the reference image
(`reference_detail_level.jpg`), which is a detail reference, **not** a texture.

Pipeline (`terrain.py` builds the field, `generate.py` meshes it):
1. **Macro shape** from the contours via distance-transform interpolation:
   ocean rim −1500 m → coast 0 m → interior rising toward the peak.
2. **Procedural detail** — domain-warped **ridged** multifractal + fractal
   (fBm) noise, amplitude strong on land (esp. upper slopes) and gentle
   undersea, for rugged ridges and rolling relief.
3. **Hydraulic erosion** — ~325k water droplets spawned on land carve realistic
   **dendritic valley networks** down the flanks.
4. Normalised so the peak is exactly **+1500 m**, the ocean floor exactly
   **−1500 m**, and sea level stays at **0**.

The mesh is **clipped to the outer contour** (the model limit) and closed into
a **watertight solid**: terrain surface on top, a vertical skirt around the rim,
and a flat base just below the ocean floor.

Horizontal size ≈ **18.8 × 20 km** (`TARGET_MAX_EXTENT_M` in `terrain.py`;
only the vertical scale is fixed by the spec). Grid resolution `N = 620`
(≈ 32 m/cell) → ~261k vertices. Increase `N` for finer detail (bigger file).

## Run
```
python3 -m pip install pymupdf trimesh numpy scipy matplotlib pillow networkx
python3 generate.py          # writes island_terrain.obj + .glb
python3 terrain.py           # optional: elevation/hillshade/section check
```

## Files
- `island_terrain.obj` — **the deliverable** (geometry only, ~23 MB)
- `island_terrain.glb` — compact single-file version for quick viewing
- `terrain.py` — heightfield generator (noise + erosion; tune params at top)
- `generate.py` — meshing, octagon clip, watertight solid, export
- `contours.json` — macro contour polygons from the guideline PDF
- `reference_detail_level.jpg` — the detail-level reference (not a texture)
- `verify_3d.png` / `verify_field.png` — self-check renders
