# Island Terrain — 3D Model

A textured 3D island terrain (mountains + ocean floor) generated from an
elevation-contour guideline PDF and a satellite texture. **Units: meters.
Altitude on +Y. Sea level at Y = 0.**

## Elevation scheme
The guideline PDF is a plan-view contour map with four nested regions. Each
region boundary is treated as an iso-elevation contour:

| contour (PDF region) | elevation |
|----------------------|-----------|
| outer polygon (blue / ocean) — model domain edge | **−1500 m** (lowest ocean floor) |
| coastline (pink edge) | **0 m** (sea level) |
| red edge | +500 m |
| bright-red edge | +1000 m |
| mountain core interior (highest point) | **+1500 m** (highest mountain) |

Between contours the surface is interpolated smoothly with Euclidean
distance transforms (distance-ratio to the two bounding contours), so the
terrain is continuous rather than stepped. A light Gaussian smoothing removes
raster stair-stepping; the highest point is pinned to exactly +1500 m and the
ocean floor bottoms out at exactly −1500 m.

## Geometry
Watertight terrain tile: the elevation surface on top, a vertical skirt around
the perimeter, and a flat base just below the ocean floor. The satellite image
(`island_terrain_texture.jpg`) is planar-projected (top-down) as the surface
texture. Default horizontal size ≈ **18.8 × 20 km** (`TARGET_MAX_EXTENT_M`,
edit freely — only the vertical scale is fixed by the spec).

## Run
```
python3 -m pip install pymupdf trimesh numpy scipy shapely pillow networkx
python3 generate.py
```
Outputs:
- `island_terrain.obj` + `material.mtl` + `material_0.png` — textured OBJ
- `island_terrain.glb` — compact single-file version (good for quick viewing)

## Files
- `generate.py` — parametric generator (elevations, extent, resolution at top)
- `contours.json` — contour polygons extracted from the guideline PDF
- `island_terrain_texture.jpg` — source satellite texture
- `render.html` — three.js viewer used for the textured verification render
- `verify_field.png` — elevation map / hillshade / cross-section check
- `verify_textured_3d.png` — textured 3D render
