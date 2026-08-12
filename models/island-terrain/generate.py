#!/usr/bin/env python3
"""
Island terrain 3D model.

Built from two inputs:
  * TectonicDashboard_ElevationGuidelines_ForAI.pdf -- a plan-view contour map.
    Four nested regions define the elevation zones:
        blue  (outer polygon) ... ocean / model domain edge
        pink  ................... coastline (sea level)
        red   ................... mid slope
        bright-red core ......... the high mountains
  * island_terrain_texture.jpg -- satellite texture, planar-projected on top.

Elevation scheme (UNITS = METERS, altitude on +Y, sea level Y = 0):
    domain boundary (outer polygon) : -1500   (lowest ocean floor)
    coastline       (pink edge)     :     0   (sea level)
    red edge                        :  +500
    bright-red edge                 : +1000
    mountain peak (core interior)   : +1500   (highest mountain)

Between contours the surface is interpolated smoothly using Euclidean
distance transforms (ratio of distance to the bounding contours), giving a
continuous terrain rather than stepped terraces.

The result is a watertight terrain tile: the elevation surface on top, a skirt
around the perimeter, and a flat base just below the ocean floor.

Only the vertical scale is fixed by the spec. TARGET_MAX_EXTENT_M sets the
horizontal size (island ~= a real volcanic island; edit freely).
"""
import json, numpy as np, trimesh
from PIL import Image
from scipy.ndimage import distance_transform_edt

# ---------- parameters ----------
PEAK_Y   =  1500.0   # highest mountain
FLOOR_Y  = -1500.0   # lowest ocean floor
SEA_Y    =     0.0   # sea level
LAND_LEVELS = [0.0, 500.0, 1000.0]   # elevations at pink / red / brightred edges
TARGET_MAX_EXTENT_M = 20000.0        # horizontal size (largest domain dimension)
GRID = 420                           # grid resolution along the larger axis
BASE_BELOW = 100.0                   # base plane sits this far below FLOOR_Y
TEXTURE = "island_terrain_texture.jpg"

# ---------- load + place contours in model space (meters, Y up) ----------
C = json.load(open("contours.json"))
dom = np.array(C["domain"])
cx, cy = dom[:,0].mean(), dom[:,1].mean()
x0,y0,x1,y1 = dom[:,0].min(),dom[:,1].min(),dom[:,0].max(),dom[:,1].max()
scale = TARGET_MAX_EXTENT_M / max(x1-x0, y1-y0)      # meters per pdf-point

def to_model(pl):
    pl = np.asarray(pl)
    X = (pl[:,0]-cx)*scale
    Z = (cy-pl[:,1])*scale                            # flip y -> north up
    return np.column_stack([X, Z])

domain = to_model(C["domain"]); pink = to_model(C["pink"])
red    = to_model(C["red"]);    peak = to_model(C["peak"])

Xmin,Xmax = domain[:,0].min(),domain[:,0].max()
Zmin,Zmax = domain[:,1].min(),domain[:,1].max()
aspect = (Zmax-Zmin)/(Xmax-Xmin)
if aspect >= 1: nz, nx = GRID, max(2,int(round(GRID/aspect)))
else:           nx, nz = GRID, max(2,int(round(GRID*aspect)))

xs = np.linspace(Xmin,Xmax,nx)
zs = np.linspace(Zmin,Zmax,nz)

def mask(poly):
    img = Image.new("1",(nx,nz),0)
    from PIL import ImageDraw
    d = ImageDraw.Draw(img)
    cols = (poly[:,0]-Xmin)/(Xmax-Xmin)*(nx-1)
    rows = (poly[:,1]-Zmin)/(Zmax-Zmin)*(nz-1)
    d.polygon(list(zip(cols.tolist(),rows.tolist())), fill=1)
    return np.array(img,dtype=bool)

m_dom  = mask(domain)
m_pink = mask(pink) & m_dom
m_red  = mask(red)  & m_dom
m_peak = mask(peak) & m_dom

# ---------- elevation field via distance-transform interpolation ----------
def Din(m):  return distance_transform_edt(m)          # inside -> dist to boundary
def Dout(m): return distance_transform_edt(~m)         # outside -> dist to region

E = np.full((nz,nx), FLOOR_Y, float)                   # default: ocean floor

# ocean ring: inside domain, outside land (pink)
ocean = m_dom & ~m_pink
d_coast = Dout(m_pink); d_edge = Din(m_dom)
t = d_coast/(d_coast+d_edge+1e-9)
E[ocean] = SEA_Y + t[ocean]*(FLOOR_Y-SEA_Y)

# ring pink->red : LAND_LEVELS[0] .. LAND_LEVELS[1]
r1 = m_pink & ~m_red
a = Din(m_pink); b = Dout(m_red); t = a/(a+b+1e-9)
E[r1] = LAND_LEVELS[0] + t[r1]*(LAND_LEVELS[1]-LAND_LEVELS[0])

# ring red->brightred : LAND_LEVELS[1] .. LAND_LEVELS[2]
r2 = m_red & ~m_peak
a = Din(m_red); b = Dout(m_peak); t = a/(a+b+1e-9)
E[r2] = LAND_LEVELS[1] + t[r2]*(LAND_LEVELS[2]-LAND_LEVELS[1])

# peak dome : LAND_LEVELS[2] .. PEAK_Y
dp = Din(m_peak); dpmax = dp.max() if dp.max()>0 else 1.0
E[m_peak] = LAND_LEVELS[2] + (dp[m_peak]/dpmax)*(PEAK_Y-LAND_LEVELS[2])

# gentle smoothing to remove rasterization stairstep (keeps extremes)
from scipy.ndimage import gaussian_filter
E = gaussian_filter(E, sigma=1.2)
E[m_peak] = np.maximum(E[m_peak], LAND_LEVELS[2])      # don't sink the peak band
E = np.clip(E, FLOOR_Y, PEAK_Y)
pos = E > 0                                            # pin highest point to +1500 exactly
if pos.any() and E[pos].max() > 0: E[pos] *= PEAK_Y/E[pos].max()

# ---------- build watertight mesh ----------
XX,ZZ = np.meshgrid(xs,zs)
top = np.column_stack([XX.ravel(), E.ravel(), ZZ.ravel()])
def vid(i,j): return i*nx+j                            # i=row(z), j=col(x)
faces=[]
for i in range(nz-1):
    for j in range(nx-1):
        a,b,c,d = vid(i,j),vid(i,j+1),vid(i+1,j+1),vid(i+1,j)
        faces.append((a,b,c)); faces.append((a,c,d))
V=[top]; nfix=len(top)
BASE_Y=FLOOR_Y-BASE_BELOW
# perimeter (clockwise) top-vertex ids
perim=[vid(0,j) for j in range(nx)]+[vid(i,nx-1) for i in range(1,nz)]+\
      [vid(nz-1,j) for j in range(nx-2,-1,-1)]+[vid(i,0) for i in range(nz-2,0,-1)]
base=top[perim].copy(); base[:,1]=BASE_Y
base_off=nfix
V.append(base)
P=len(perim)
for k in range(P):
    a=perim[k]; b=perim[(k+1)%P]; c=base_off+((k+1)%P); d=base_off+k
    faces.append((a,b,c)); faces.append((a,c,d))       # skirt
for k in range(1,P-1):                                 # base fan
    faces.append((base_off,base_off+k+1,base_off+k))
verts=np.vstack(V); faces=np.array(faces)

# ---------- UVs (planar top-down; texture spans domain bbox) ----------
uv=np.zeros((len(verts),2))
uv[:nfix,0]=(top[:,0]-Xmin)/(Xmax-Xmin)
uv[:nfix,1]=(top[:,2]-Zmin)/(Zmax-Zmin)
uv[base_off:,0]=np.clip((base[:,0]-Xmin)/(Xmax-Xmin),0,1)
uv[base_off:,1]=np.clip((base[:,2]-Zmin)/(Zmax-Zmin),0,1)

img=Image.open(TEXTURE).convert("RGB")
mat=trimesh.visual.texture.SimpleMaterial(image=img)
mesh=trimesh.Trimesh(vertices=verts,faces=faces,
      visual=trimesh.visual.TextureVisuals(uv=uv,image=img,material=mat),
      process=False)
mesh.fix_normals()

if __name__=="__main__":
    print(f"grid {nx} x {nz} | horiz extent {Xmax-Xmin:.0f} x {Zmax-Zmin:.0f} m")
    print(f"elevation min/max: {E.min():.0f} / {E.max():.0f} m (sea level = 0)")
    print(f"verts {len(verts)} | faces {len(faces)} | watertight {mesh.is_watertight}")
    mesh.export("island_terrain.glb")
    mesh.export("island_terrain.obj")     # writes .obj + .mtl + texture
    np.save("elevation_grid.npy", E)
    print("wrote island_terrain.obj / .mtl / .glb")
