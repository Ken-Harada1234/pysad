"""pysad - Python translation of SAD (Strategic Accelerator Design, by K. Oide, KEK), K64 version.
pysad - SAD (K64 版) の Python への翻訳。

SAD lattices and command scripts run as they are. The default `v2` mode fixes K64 bugs; `PYSAD_COMPAT=k64` reproduces K64.
SAD のラティスとコマンドは、そのまま実行できる。既定の `v2` は K64 のバグを直した動作で、`PYSAD_COMPAT=k64` で K64 と同じ動作になる。

Entry points / 入口
-------------------
pysad.run(text)                  run SAD lattice and command text given as a string. 文字列の SAD のテキストを実行する
pysad.sad_help('Twiss')          help for SAD commands, functions, variables, flags, elements. SAD の名前の説明 (`search='tune'` でキーワード検索)
%load_ext pysad.notebook         Jupyter: `%%sad` cells, `%sad_help`, `sad.value(...)`. Jupyter でラティスとコマンドをセルに貼って実行する
pysad-run, pysad-help, pysad-share, pysad-warmup   command line (also `python -m pysad <command>`). コマンドライン
Python API (Lattice, fodo_cell, track_particles, ...): see docs/12_pysad_api.md. 下位の Python API は docs/12_pysad_api.md

Quick start / 最初の例
----------------------
>>> import pysad
>>> pysad.run(r'''MOMENTUM = 3GEV;
... DRIFT D1 = (L=1.0) D2 = (L=0.5);
... QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
... BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
... LINE CELL = (QF D2 B D2 QD D1);
... LINE RING = (24*CELL);
... FFS;
... USE RING; CELL; CALC;
... Print["NX=", Twiss["NX","$$$"]];
... ''')                          # prints NX=34.01412233018365 (after the optics table)

Documents, examples and tests are not in the pip package: https://github.com/Ken-Harada1234/pysad and
https://ken-harada1234.github.io/pysad/ (`pysad-share` prints the places). 文書・例・試験は pip には入っていない (上の URL)。
Not translated (see docs/31 for the list): e.g. EPICS device control. 未翻訳のものは docs/31 の一覧を見ること。
"""

from .constants import (
    SpeedOfLight, ElectronMass, ElectronRadius,
    ProtonMass, Cgamma, Cq,
)

from .elements import (
    Element, Drift, Quadrupole, Bend, Sextupole,
    Cavity, Monitor, Marker, Multipole,
    make_element,
)

from .optics import (
    symplectic_j, symplectic_inverse,
    normal_coordinates, transfer_matrix,
    beam_matrix, periodic_r, periodic_optics,
    dispersion_from_6d, geo_base,
)

from .lattice import (
    Lattice,
    fodo_cell,
)

from .parser import (
    SADScript,
    load_sad_file,
    load_sad_string,
    load_sad_files,
    SADEnvironment,
)

from .tracking import (
    pxy2dpz,
    track_drift_exact,
    track_thin_quad,
    track_thin_sext,
    track_multipole_kick,
    track_thick_quad,
    track_thick_bend,
    track_cavity,
    radiation_loss_per_turn,
    damping_times,
    track_element,
    track_particles,
    TrackResult,
    dynamic_aperture,
    tune_by_fft,
)

from .matching import (
    Vary,
    Constraint,
    Match,
    MatchResult,
    DynamicApertureSurvey,
)

from .ffs import (
    FFSState,
    SADList,
)

__version__ = "2.6.0"
__author__ = "SAD Python translation (KEK-Acc6th)"
__all__ = [
    # constants
    "SpeedOfLight", "ElectronMass", "ElectronRadius",
    "ProtonMass", "Cgamma", "Cq",
    # elements
    "Element", "Drift", "Quadrupole", "Bend", "Sextupole",
    "Cavity", "Monitor", "Marker", "Multipole", "make_element",
    # optics
    "symplectic_j", "symplectic_inverse",
    "normal_coordinates", "transfer_matrix",
    "beam_matrix", "periodic_r", "periodic_optics",
    "dispersion_from_6d", "geo_base",
    # lattice
    "Lattice", "fodo_cell",
    # tracking
    "pxy2dpz", "track_drift_exact",
    "track_thin_quad", "track_thin_sext", "track_multipole_kick",
    "track_thick_quad", "track_thick_bend",
    "track_element", "track_particles", "TrackResult",
    "dynamic_aperture", "tune_by_fft",
    # parser
    "SADScript", "load_sad_file", "load_sad_string", "load_sad_files", "SADEnvironment",
    # matching
    "Vary", "Constraint", "Match", "MatchResult",
    "DynamicApertureSurvey",
    # ffs
    "FFSState", "SADList",
    # run SAD text from Python (notebook), SAD help
    "run", "sad_help",
]


def run(src, **kwargs):
    """Run SAD lattice/command text given as a Python string (see pysad.notebook.run; use r-strings for backslashes).
    文字列で渡した SAD のラティスとコマンドを実行する (詳しくは pysad.notebook.run。バックスラッシュは r 文字列で)。"""
    from .notebook import run as _run
    return _run(src, **kwargs)


def sad_help(name=None, search=None, full=False, kind=None, out=None):
    """Show help for a SAD/FFS/SADScript name, e.g. sad_help('Twiss') or sad_help(search='tune') (see pysad.sadhelp).
    SAD の FFS コマンド・関数・特別変数・フラグ・要素の説明を出す。名前で引く (`sad_help('Twiss')`)、キーワードで探す (`search='tune'`)、引数なしで概要。"""
    from .sadhelp import sad_help as _sad_help
    return _sad_help(name, search, full, kind, out)
