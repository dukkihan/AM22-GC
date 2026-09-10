#!/usr/bin/env python3
"""
Site location map for sediment coring stations on the eastern Amundsen Sea
continental shelf (Getz and Dotson ice-shelf sectors), West Antarctica.

Bathymetry : IBCSO v2, 500 m, EPSG:9354 (Dorschel et al., 2022)
             https://doi.org/10.1594/PANGAEA.937574          CC-BY 4.0
Coastline,
ice shelves: Natural Earth 1:10 m physical vectors (public domain)
             https://www.naturalearthdata.com
Projection : polar stereographic on WGS 84, true scale at 71 S (as EPSG:3031),
             rotated so that the 113.7 W meridian is vertical.

Required input files in the working directory
    IBCSO_v2_bed.tif
    ne/ne_10m_coastline.{shp,shx,dbf}
    ne/ne_10m_antarctic_ice_shelves_polys.{shp,shx,dbf}
    ne/ne_10m_glaciated_areas.{shp,shx,dbf}

Dependencies : numpy, matplotlib, pyproj, rasterio, pyshp
Usage        : python3 make_site_map.py
Output       : amundsen_site_map.png (600 dpi), amundsen_site_map.pdf
"""
import os
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
import shapefile                                   # pyshp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle, Patch
from matplotlib.collections import PathCollection
from matplotlib.colors import LightSource, Normalize, LinearSegmentedColormap
from pyproj import CRS, Transformer

# ----------------------------------------------------------------- settings
IBCSO     = "IBCSO_v2_bed.tif"
NE        = "ne"
ICE_SHELF = "ne_10m_antarctic_ice_shelves_polys"   # swap for SCAR ADD / MEaSUREs
GROUNDED  = "ne_10m_glaciated_areas"
COASTLINE = "ne_10m_coastline"

LON0    = -113.7                 # meridian pointing up
CENTER  = (-113.7, -73.8)        # map centre (lon, lat)
HALF    = (142_000, 116_000)     # half-width, half-height (m)
RES     = 500                    # output pixel size (m)

DEPTHS      = [-1200, -900, -600, -300]   # isobaths to draw
LABELLED    = [-900, -600]                # isobaths that get a number
VERT_EXAG   = 12                          # hillshade vertical exaggeration
DEPTH_RANGE = (-1500, -150)               # colour scale (vmin, vmax) in m
GRAT_LON    = np.arange(-122, -105, 2)
GRAT_LAT    = np.arange(-76, -71.9, 1)

SITES = [("Getz1",   -73.982820, -115.769063),
         ("Getz2",   -73.504947, -114.969610),
         ("Dotson1", -74.130567, -112.550132),
         ("Dotson2", -73.652503, -113.631410)]

# label offsets in metres: (dx, dy, horizontal alignment)
SITE_LABEL_OFFSET = {"Getz1":   (6500, 5000, "left"),
                     "Getz2":   (6500, 5000, "left"),
                     "Dotson1": (-6500, 5000, "right"),
                     "Dotson2": (-6500, 5000, "right")}

PLACE_NAMES = [("Getz\nIce Shelf",     -74.42, -115.7),
               ("Dotson\nIce Shelf",  -74.48, -112.7),
               ("Amundsen Sea",       -73.05, -116.3)]

DRAW_LABELS = True               # False -> symbols only, text added by hand later

# panel placement in figure coordinates [left, bottom, width, height]
CBAR_BOX   = [.088, .700, .018, .205]     # vertical colour bar, upper left
LEGEND_LOC = "lower right"                # ice shelf / grounded ice key

BBOX = (-128, -79, -100, -70)    # lon/lat pre-filter for the vector layers

MAPCRS = CRS.from_proj4(f"+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0={LON0} "
                        f"+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
T = Transformer.from_crs("EPSG:4326", MAPCRS, always_xy=True)

cx, cy = T.transform(*CENTER)
EXT = (cx - HALF[0], cx + HALF[0], cy - HALF[1], cy + HALF[1])


# ----------------------------------------------------------- vector helpers
def read_shapes(name, bbox=BBOX):
    """Read a Natural Earth shapefile, keeping only shapes near the map area."""
    sf = shapefile.Reader(os.path.join(NE, name))
    out = []
    for s in sf.shapes():
        if not s.points:
            continue
        x0, y0, x1, y1 = s.bbox
        if x1 < bbox[0] or x0 > bbox[2] or y1 < bbox[1] or y0 > bbox[3]:
            continue
        out.append((np.asarray(s.points, float), list(s.parts)))
    return out


def to_path(pts, parts):
    """Project lon/lat vertices and build a (possibly multi-part) Path."""
    xy = np.column_stack(T.transform(pts[:, 0], pts[:, 1]))
    codes = np.full(len(xy), Path.LINETO, np.uint8)
    for p in parts:
        codes[p] = Path.MOVETO
    return Path(xy, codes)


def add_polys(ax, name, bbox=BBOX, **kw):
    paths = [to_path(p, q) for p, q in read_shapes(name, bbox)]
    if paths:
        ax.add_collection(PathCollection(paths, **kw))


def add_lines(ax, name, bbox=BBOX, **kw):
    for p, q in read_shapes(name, bbox):
        ax.add_patch(PathPatch(to_path(p, q), fill=False, **kw))


# ------------------------------------------------------- raster: clip + warp
def load_grid():
    """Reproject the IBCSO window straight from EPSG:9354 into the map CRS."""
    w = int((EXT[1] - EXT[0]) / RES)
    h = int((EXT[3] - EXT[2]) / RES)
    dst = np.full((h, w), np.nan, "float32")
    dst_tr = from_origin(EXT[0], EXT[3], RES, RES)
    with rasterio.open(IBCSO) as src:
        reproject(source=rasterio.band(src, 1), destination=dst,
                  src_transform=src.transform, src_crs=src.crs,
                  src_nodata=src.nodata,
                  dst_transform=dst_tr, dst_crs=MAPCRS, dst_nodata=np.nan,
                  resampling=Resampling.cubic)
    return dst


def ice_mask(shape):
    """True where a grid cell lies under an ice shelf or grounded ice."""
    gx = np.linspace(EXT[0] + RES / 2, EXT[1] - RES / 2, shape[1])
    gy = np.linspace(EXT[3] - RES / 2, EXT[2] + RES / 2, shape[0])
    X, Y = np.meshgrid(gx, gy)
    pts = np.column_stack([X.ravel(), Y.ravel()])
    m = np.zeros(len(pts), bool)
    for nm in (ICE_SHELF, GROUNDED):
        for p, q in read_shapes(nm):
            path = to_path(p, q)
            bb = path.get_extents()
            if bb.x1 < EXT[0] or bb.x0 > EXT[1] or bb.y1 < EXT[2] or bb.y0 > EXT[3]:
                continue
            todo = ~m
            if todo.any():
                m[todo] = path.contains_points(pts[todo])
    return m.reshape(shape)


# --------------------------------------------------------------- decoration
def graticule(ax, lons, lats):
    for la in lats:
        lo = np.linspace(lons[0] - 6, lons[-1] + 6, 400)
        ax.plot(*T.transform(lo, np.full_like(lo, la)), lw=.5,
                color="#6f8490", alpha=.8, zorder=4)
    for lo in lons:
        la = np.linspace(lats[0] - 3, lats[-1] + 3, 300)
        ax.plot(*T.transform(np.full_like(la, lo), la), lw=.5,
                color="#6f8490", alpha=.8, zorder=4)


def edge_labels(ax, lons, lats):
    """Label the graticule where each line crosses the map frame."""
    x0, x1, y0, y1 = EXT
    st = dict(fontsize=7.5, color="#2d3b44", zorder=8)

    def hits(xs, ys):
        out = []
        for i in range(len(xs) - 1):
            for a, b, f, side in ((x0, x1, y0, "b"), (x0, x1, y1, "t"),
                                  (y0, y1, x0, "l"), (y0, y1, x1, "r")):
                if side in "bt" and (ys[i] - f) * (ys[i + 1] - f) < 0:
                    t = (f - ys[i]) / (ys[i + 1] - ys[i])
                    xx = xs[i] + t * (xs[i + 1] - xs[i])
                    if a <= xx <= b:
                        out.append((xx, f, side))
                if side in "lr" and (xs[i] - f) * (xs[i + 1] - f) < 0:
                    t = (f - xs[i]) / (xs[i + 1] - xs[i])
                    yy = ys[i] + t * (ys[i + 1] - ys[i])
                    if a <= yy <= b:
                        out.append((f, yy, side))
        return out

    for la in lats:
        lo = np.linspace(-130, -100, 600)
        for xx, yy, s in hits(*T.transform(lo, np.full_like(lo, la))):
            if s == "l":
                ax.text(xx - 4500, yy, f"{abs(la):g}\u00b0S", ha="right", va="center", **st)
            elif s == "r":
                ax.text(xx + 4500, yy, f"{abs(la):g}\u00b0S", ha="left", va="center", **st)
    for lo in lons:
        la = np.linspace(-79, -70, 600)
        for xx, yy, s in hits(*T.transform(np.full_like(la, lo), la)):
            if s == "b":
                ax.text(xx, yy - 4500, f"{abs(lo):g}\u00b0W", ha="center", va="top", **st)
            elif s == "t":
                ax.text(xx, yy + 4500, f"{abs(lo):g}\u00b0W", ha="center", va="bottom", **st)


def scalebar(ax, km=50, lat_ref=-74.0):
    """Scale bar corrected for the polar-stereographic scale factor at lat_ref."""
    x0, x1, y0, y1 = EXT
    k = (1 + np.sin(np.radians(71))) / (1 + np.sin(np.radians(-lat_ref)))
    L = km * 1000 / k
    bx, by = x0 + .055 * (x1 - x0), y0 + .07 * (y1 - y0)
    hgt = .011 * (y1 - y0)
    for i in range(2):
        ax.add_patch(Rectangle((bx + i * L / 2, by), L / 2, hgt,
                               facecolor="#12212a" if i == 0 else "white",
                               edgecolor="#12212a", lw=.7, zorder=9))
    for f, lab in ((0, "0"), (.5, str(km // 2)), (1, str(km))):
        ax.text(bx + f * L, by + hgt * 1.7, lab, ha="center", va="bottom",
                fontsize=7, color="#12212a", zorder=9)
    ax.text(bx + L * 1.06, by + hgt * .5, "km", ha="left", va="center",
            fontsize=7, color="#12212a", zorder=9)


# --------------------------------------------------------------------- main
def main():
    z = load_grid()
    z = np.where(z > 0, np.nan, z)          # drop anything above sea level
    z[ice_mask(z.shape)] = np.nan           # drop sea floor beneath ice

    cmap = LinearSegmentedColormap.from_list(
        "ibcso", ["#0b3d6b", "#1b5f96", "#3785bd", "#6fb0d8", "#a8d2ea", "#d6ecf7"])
    norm = Normalize(vmin=DEPTH_RANGE[0], vmax=DEPTH_RANGE[1])

    fig = plt.figure(figsize=(7.2, 5.9), dpi=600)
    ax = fig.add_axes([.055, .055, .90, .90])
    ax.set_facecolor("#eef4f8")

    ls = LightSource(azdeg=315, altdeg=45)
    rgb = ls.shade(np.ma.masked_invalid(z), cmap=cmap, norm=norm,
                   blend_mode="soft", vert_exag=VERT_EXAG, dx=RES, dy=RES)
    ax.imshow(rgb, extent=(EXT[0], EXT[1], EXT[2], EXT[3]), origin="upper",
              zorder=1, interpolation="bilinear")

    gx = np.linspace(EXT[0], EXT[1], z.shape[1])
    gy = np.linspace(EXT[3], EXT[2], z.shape[0])
    cs = ax.contour(gx, gy, z, levels=DEPTHS, colors="#2f4b5c",
                    linewidths=.3, alpha=.6, zorder=2)
    ax.clabel(cs, levels=LABELLED, fmt="%d", fontsize=5, inline_spacing=3)

    add_polys(ax, GROUNDED, facecolor="#f5f5f3", edgecolor="#b4babd",
              linewidths=.4, zorder=3)
    add_polys(ax, ICE_SHELF, facecolor="#e3ebef", edgecolor="#8f9fa8",
              linewidths=.5, zorder=3.2)
    add_lines(ax, COASTLINE, edgecolor="#55636a", lw=.6, zorder=3.4)

    graticule(ax, GRAT_LON, GRAT_LAT)

    if DRAW_LABELS:
        for lab, la, lo in PLACE_NAMES:
            ax.text(*T.transform(lo, la), lab, fontsize=8, style="italic",
                    color="#3f5460", ha="center", va="center", zorder=6)

    for name, la, lo in SITES:
        x, y = T.transform(lo, la)
        ax.plot(x, y, "o", ms=6.5, mfc="#c8102e", mec="white", mew=1.2, zorder=7)
        if DRAW_LABELS:
            dx, dy, ha = SITE_LABEL_OFFSET[name]
            ax.text(x + dx, y + dy, name, fontsize=8.5, fontweight="bold",
                    color="#12212a", ha=ha, va="bottom", zorder=8)

    if DRAW_LABELS:
        leg = ax.legend(
            handles=[Patch(fc="#e3ebef", ec="#8f9fa8", lw=.5, label="Ice shelf"),
                     Patch(fc="#f5f5f3", ec="#b4babd", lw=.5, label="Grounded ice")],
            loc=LEGEND_LOC, fontsize=7, frameon=True, framealpha=.94,
            edgecolor="#2d3b44", borderpad=.6, handlelength=1.4)
        leg.get_frame().set_linewidth(.6)
        leg.set_zorder(9)

    scalebar(ax)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_linewidth(.8); s.set_color("#2d3b44")
    edge_labels(ax, GRAT_LON, GRAT_LAT)

    cax = fig.add_axes(CBAR_BOX)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                      cax=cax, orientation="vertical")
    cb.set_ticks([DEPTH_RANGE[0], int(sum(DEPTH_RANGE) / 2), DEPTH_RANGE[1]])
    cb.ax.tick_params(labelsize=6.5, length=2, pad=2)
    cb.ax.yaxis.set_ticks_position("right")
    cb.ax.yaxis.set_label_position("left")     # axis title on the opposite side
    cb.outline.set_linewidth(.5)
    if DRAW_LABELS:
        cb.set_label("Depth (m)", fontsize=6.5, labelpad=4)

    ins = fig.add_axes([.735, .735, .215, .215])
    add_polys(ins, GROUNDED, (-180, -90, 180, -58), facecolor="#e9edef",
              edgecolor="#98a3a9", linewidths=.2)
    ins.add_patch(Rectangle((EXT[0], EXT[2]), EXT[1] - EXT[0], EXT[3] - EXT[2],
                            fill=False, ec="#c8102e", lw=1.0))
    ins.plot(*T.transform(*CENTER), "s", ms=3.2, color="#c8102e")
    R = 3.1e6
    ins.set_xlim(-R, R); ins.set_ylim(-R, R)
    ins.set_xticks([]); ins.set_yticks([]); ins.set_aspect("equal")
    ins.set_facecolor("white")
    for s in ins.spines.values():
        s.set_linewidth(.6); s.set_color("#2d3b44")

    for e in ("png", "pdf"):
        fig.savefig(f"amundsen_site_map.{e}", dpi=600, bbox_inches="tight",
                    facecolor="white")
    print("done: amundsen_site_map.png / .pdf")


if __name__ == "__main__":
    main()
