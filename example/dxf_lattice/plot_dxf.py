# -*- coding: utf-8 -*-
"""DXF (pysad.dxf_lattice.write_dxf / DrawMagnetsDXF の出力) を読み戻して matplotlib で描き、PNG に保存する (目視の確認用)。
使い方: py -3.12 example/dxf_lattice/plot_dxf.py in.dxf out.png [title] [xmin xmax ymin ymax]"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysad.dxf_lattice import read_dxf

ACI = {1: '#d62728', 2: '#bcbd22', 3: '#2ca02c', 4: '#17becf', 5: '#1f4fd8', 6: '#c71585', 7: '#222222', 8: '#888888', 30: '#ff7f0e'}


def plot(dxf, png, title='', lim=None, labels=True, dpi=150, size=(9, 9)):
    d = read_dxf(dxf)
    fig, ax = plt.subplots(figsize=size)
    for e in d['entities']:
        col = ACI.get(d['layers'].get(e['layer'], 7), 'k')
        if e['type'] == 'POLYLINE':
            xs = [p[0] for p in e['pts']]; ys = [p[1] for p in e['pts']]
            if e['closed']:
                ax.fill(xs, ys, facecolor=col, alpha=0.35, edgecolor=col, lw=0.8)
            else:
                ax.plot(xs, ys, color=col, lw=0.6)
        elif e['type'] == 'LINE':
            ax.plot([e['pts'][0][0], e['pts'][1][0]], [e['pts'][0][1], e['pts'][1][1]], color=col, lw=0.8)
        elif e['type'] == 'TEXT' and labels:
            ax.text(e['at'][0], e['at'][1], e['text'], fontsize=5, rotation=e.get('rot', 0), ha='center', va='center', color=col)
    ax.set_aspect('equal')
    if lim:
        ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3])
    ax.set_xlabel('X [%s]' % {6: 'm', 4: 'mm', 5: 'cm'}.get(d['insunits'], '?')); ax.set_ylabel('Y')
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(png, dpi=dpi); plt.close(fig)
    return d


if __name__ == '__main__':
    a = sys.argv
    lim = [float(v) for v in a[4:8]] if len(a) >= 8 else None
    plot(a[1], a[2], a[3] if len(a) > 3 else '', lim)
