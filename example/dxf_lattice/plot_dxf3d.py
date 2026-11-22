# -*- coding: utf-8 -*-
"""3D の DXF (write_dxf(..., mode='3d') / DrawMagnetsDXF[..., Mode->"3D"] の出力) を読み戻して matplotlib の mplot3d で描き、PNG に保存する
(目視の確認用)。3DFACE をレイヤーの色 (ACI) の面にして描き、基準軌道 (z = 0) を線で重ねる。軸は x・y・z とも同じ縮尺。
使い方: py -3.12 example/dxf_lattice/plot_dxf3d.py in.dxf out.png [title]
関数: plot3d(dxf, png, title, panels=[(仰角, 方位角, 範囲 (xmin, xmax, ymin, ymax) または None), ...], size, dpi)
      各 panel が 1 つの 3D の図 (横に並べる)。範囲を与えると、重心がその中にある面だけを描き、軸もその範囲にする"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pysad.dxf_lattice import read_dxf

ACI = {1: '#d62728', 2: '#bcbd22', 3: '#2ca02c', 4: '#17becf', 5: '#1f4fd8', 6: '#c71585', 7: '#222222', 8: '#888888', 9: '#b0b0b0',
       30: '#ff7f0e'}


def plot3d(dxf, png, title='', panels=((25, -60, None),), size=(14, 6), dpi=130):
    d = read_dxf(dxf)
    faces = [e for e in d['entities'] if e['type'] == '3DFACE']
    trajs = [e for e in d['entities'] if e['type'] == 'POLYLINE' and e['layer'] == 'TRAJECTORY']
    fig = plt.figure(figsize=size)
    for k, (elev, azim, lim) in enumerate(panels):
        ax = fig.add_subplot(1, len(panels), k + 1, projection='3d')
        sel = faces
        if lim:      # 6 枚ずつの立体の重心が範囲の中にあるものだけ (立体の途中で切らない)
            sel = []
            for k0 in range(0, len(faces), 6):
                chunk = faces[k0:k0 + 6]
                cx = sum(p[0] for f in chunk for p in f['pts']) / (4 * len(chunk))
                cy = sum(p[1] for f in chunk for p in f['pts']) / (4 * len(chunk))
                if lim[0] <= cx <= lim[1] and lim[2] <= cy <= lim[3]:
                    sel += chunk
        bylayer = {}
        for f in sel:
            bylayer.setdefault(f['layer'], []).append([tuple(p) for p in f['pts']])
        for lay, polys in bylayer.items():
            col = ACI.get(d['layers'].get(lay, 7), 'k')
            ax.add_collection3d(Poly3DCollection(polys, facecolors=col, edgecolors='#222222', linewidths=0.25, alpha=0.92))
        pts = [p for f in sel for p in f['pts']]
        if not pts:
            continue
        for t in trajs:
            xs, ys = [p[0] for p in t['pts']], [p[1] for p in t['pts']]
            if lim:
                keep = [(x, y) for x, y in zip(xs, ys) if lim[0] <= x <= lim[1] and lim[2] <= y <= lim[3]]
                xs, ys = [q[0] for q in keep], [q[1] for q in keep]
            ax.plot(xs, ys, [0.0] * len(xs), color='k', lw=0.8)
        x0, x1 = (lim[0], lim[1]) if lim else (min(p[0] for p in pts), max(p[0] for p in pts))
        y0, y1 = (lim[2], lim[3]) if lim else (min(p[1] for p in pts), max(p[1] for p in pts))
        z0, z1 = min(p[2] for p in pts), max(p[2] for p in pts)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_zlim(z0, z1)
        ax.set_box_aspect((x1 - x0, y1 - y0, z1 - z0))
        ax.view_init(elev=elev, azim=azim)
        u = {6: 'm', 4: 'mm', 5: 'cm'}.get(d['insunits'], '?')
        ax.set_xlabel('X [%s]' % u)
        ax.set_ylabel('Y [%s]' % u)
        ax.set_zlabel('Z [%s]' % u)
        ax.tick_params(labelsize=6)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=dpi)
    plt.close(fig)
    return d


if __name__ == '__main__':
    a = sys.argv
    plot3d(a[1], a[2], a[3] if len(a) > 3 else '')
