# Plate Tectonics Counter — 3D Model

Simple watertight 3D solid of the 8-sided museum exhibit table from drawing
**FM-EX-315 (REV 1, 13 JULY 26)**. Units: **inches**.

## How it was built
The plan view is a scaled CAD drawing. The vector geometry was extracted
directly (rather than eyeballing fractions): scale came out to exactly
**4.5 pt = 1 in**, verified against all eight labeled outer edges
(7'-2⅝, 2'-11¾, 3'-10, 3'-10⁹⁄₁₆, 4'-8³⁄₁₆, 7'-3½, 4'-3⅝, 4'-8¹³⁄₁₆).

In plan there are **three** octagon outlines (the "second inner line" and the
edge of the central area are the same line in plan, because the 2" drop is
vertical):

| ring | meaning | height z |
|------|---------|----------|
| A — outer outline | lowest point of the top surface | 0 |
| B — first inner line | top of the 7.75" rising border | 7.75 |
| C — second inner line | flat lip ends; edge of the vertical drop | 7.75 |
| central area (inside C) | recessed flat floor after the 2" drop | 5.75 |

Four **24"-class touchscreens** (bezel opening ~20.9 × 10.3 in) are set into
the rising border, centered on alternating edges (1, 3, 5, 7), as **blind
recesses** cut perpendicular to the ramp so a flat panel sits flush with the
slope.

The model is kept deliberately simple: 8-vertex outer polygon, no legs, just
the rising border, the flat recessed center, and the screen recesses.

Two values are **not** on the drawing (parameters at the top of `generate.py`):
- `BASE_THICKNESS` (default **4 in**) — slab thickness below the outer rim, so
  the result is a closed solid. Does not affect the top surface.
- `SCREEN_RECESS_DEPTH` (default **2 in**) — pocket depth for the screens,
  measured perpendicular to the sloped border.

## Run
```
python3 -m pip install trimesh manifold3d numpy scipy
python3 generate.py
```
Outputs `plate_tectonics_counter.obj`, `.stl`, and `.glb`. Overall footprint
≈ 12'-0" × 12'-8", height 11.75" (7.75 rise + 4 base).

## Files
- `generate.py` — parametric generator (edit dimensions/heights/base/recess here)
- `plate_tectonics_counter.obj` / `.stl` / `.glb` — the model
- `verification_render.png` — three views used to self-check the geometry
