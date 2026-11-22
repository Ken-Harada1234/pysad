"""
SAD の Packages/*.n を SADScript のまま使うための自動読み込み。
Packages ディレクトリは環境変数 PYSAD_SAD_PACKAGES、なければ pysad の隣の SAD-master/Packages (見つからなければ何もしない)。

- correction.n・bump.n・corutility.n (CorrectOrbit、MakeBump、Monitor、Steer ...): そのファイルが定義する関数名
  (`Name[...]:=` と `Name=` の最上位) を登録する (init.n の AutoLoad の表より広い。correction.n は bump.n と corutility.n も読む)。
- そのほか init.n の AutoLoad[名前..., PackagesDirectory//"file.n"] の表のうち、下の許可リストのファイル (D.n、rational.n、tableform.n、LDSolve.n、
  Couple.n、optimize.n、emit.n、Optics.n) の名前を登録する。GUI・EPICS・Canvas 系は読まない。
- 最初に使われたときだけ読む。pysad が組み込みで持つ関数は、パッケージの定義で上書きしない (読んだ後で組み込みと同名の定義を捨てる)。
"""
from __future__ import annotations

import os
import re

from .sadscript import Sym

_FILES = ('correction.n', 'bump.n', 'corutility.n')
_TARGET = 'correction.n'          # この 1 本が bump.n と corutility.n を Get する
_DEF = re.compile(r'^([A-Z][A-Za-z0-9$]*)\s*(?:\[[^\n]*\]\s*:=|=)', re.M)
_UPDEF = re.compile(r'^\(\s*[a-z][A-Za-z0-9$]*_\s*=\s*([A-Z][A-Za-z0-9$]*)\[', re.M)      # (a_=Bump[…])^:= … (UpValue で定義する関数名)
_SKIP = {'Protect', 'Chop'}
_ALLOW = ('D.n', 'rational.n', 'tableform.n', 'LDSolve.n', 'Couple.n', 'optimize.n', 'emit.n', 'Optics.n',
          'ListPlot.n', 'Plots.n', 'OpticsPlot.n', 'GeometryPlot.n', 'CSR.n')
_INIT = re.compile(r'AutoLoad(?:Set)?\[(.*?)PackagesDirectory//"([^"]+)"\s*\]', re.S)


def packages_dir():
    p = os.environ.get('PYSAD_SAD_PACKAGES')
    if p and os.path.isdir(p):
        return p.rstrip('/' + chr(92)) + '/'
    here = os.path.dirname(os.path.abspath(__file__))
    for rel in ('..', '../..'):
        for base in (os.path.join(here, rel, 'SAD-master', 'Packages'), os.path.join(here, rel, '..', 'SAD-master', 'Packages')):
            if os.path.isfile(os.path.join(base, _TARGET)):
                return os.path.abspath(base).replace(chr(92), '/') + '/'
    return None


def register(I):
    pk = packages_dir()
    if pk is None or not os.path.isfile(pk + _TARGET):
        return False
    if 'PackagesDirectory' not in I.globals:
        I.globals['PackagesDirectory'] = pk
    n = 0
    path = pk + _TARGET
    for f in _FILES:
        try:
            src = open(pk + f, encoding='utf-8', errors='replace').read()
        except OSError:
            continue
        for m in _DEF.finditer(src):
            name = m.group(1)
            if name in _SKIP or name in I.builtins or name in I.downvalues:
                continue
            I.autoload.setdefault(name, path)
            n += 1
        for m in _UPDEF.finditer(src):
            name = m.group(1)
            if name not in _SKIP and name not in I.builtins and name not in I.downvalues:
                I.autoload.setdefault(name, path)
                n += 1
    try:
        init = open(pk + 'init.n', encoding='utf-8', errors='replace').read()
    except OSError:
        init = ''
    for m in _INIT.finditer(init):
        f = m.group(2)
        if f not in _ALLOW or not os.path.isfile(pk + f):
            continue
        for nm in m.group(1).split(','):
            nm = nm.strip().strip('`')
            if re.fullmatch(r'[A-Za-z][A-Za-z0-9$]*', nm) and nm not in I.builtins and nm not in I.downvalues:
                I.autoload.setdefault(nm, pk + f)
                n += 1
    for nm in ('EVCompile', 'EVFun', 'EVMatrix', 'ClearEV', 'ElementValueExceptions', 'IgnoreCouplingMatrix'):      # Couple.n の補助 (init.n の表にあるのは SetCoupledElements・CouplingMatrix・CoupledKeys だけ)
        if os.path.isfile(pk + 'Couple.n') and nm not in I.builtins and nm not in I.downvalues:
            I.autoload.setdefault(nm, pk + 'Couple.n')
            n += 1
    # init.n の `Get[PackagesDirectory//"PkgConfig.n"]` (OpticsPlot.n などが PkgConfig$DefaultList / PkgConfig を使う): 最初に使われたときに読む
    if os.path.isfile(pk + 'PkgConfig.n'):
        for nm in ('PkgConfig', 'PkgConfigList', 'PkgConfigMatch', 'PkgConfig$DefaultList', 'PkgConfig$TypeList',
                   'PkgConfig$Default', 'PkgConfig$Type'):
            if nm not in I.builtins and nm not in I.downvalues:
                I.autoload.setdefault(nm, pk + 'PkgConfig.n')
                n += 1
    # $HaveLaTeX は SAD 本体が latex・dvips を探して設定する (見つからなければ "")。pysad は LaTeX の描画を持たないので "" (ユーザーが設定すれば尊重する)
    if '$HaveLaTeX' not in I.globals:
        I.globals['$HaveLaTeX'] = ''
    # init.n: $DisplayFunction=If[Environment["DISPLAY"]<=>"", CanvasDrawer, Identity]。pysad は Windows などの画面があるときも CanvasDrawer
    # (tkinter に描く。tclbridge が Canvas.n を読む)。画面なし (PYSAD_NO_DISPLAY=1、Linux で DISPLAY なし) は Identity (Graphics を返すだけ)
    # v2.4.1: 組み込み関数 MatplotlibDrawer (mpl_display.py。Graphics を matplotlib に描く) を足す。環境変数 PYSAD_DISPLAY=mpl で
    # $DisplayFunction をそれにする (tk = CanvasDrawer、none = Identity)。未設定は従来どおり
    from . import mpl_display
    mpl_display.install(I)
    if '$DisplayFunction' not in I.globals:
        mode = mpl_display.display_mode_from_env()
        if mode == 'mpl':
            I.globals['$DisplayFunction'] = Sym('MatplotlibDrawer')
        elif mode in ('tk', 'none'):
            I.globals['$DisplayFunction'] = Sym('CanvasDrawer' if mode == 'tk' else 'Identity')
        else:
            gui = not os.environ.get('PYSAD_NO_DISPLAY') and (os.name == 'nt' or os.environ.get('DISPLAY'))
            I.globals['$DisplayFunction'] = Sym('CanvasDrawer' if gui else 'Identity')
    _protect_builtins(I)
    return n > 0


def _protect_builtins(I):
    """自動読み込みしたパッケージが、pysad の組み込みと同名の関数 (LinearSolve、Binomial ...) の定義を持ち込んでも捨てる"""
    if getattr(I, '_pkg_protect', False):
        return
    I._pkg_protect = True
    orig = I._load_autoload

    I._pkg_keep = set()

    def load(name, orig=orig):
        path_ = I.autoload.get(name)
        before_b = set(I.builtins)
        before_d = {k: list(v) for k, v in I.downvalues.items() if k in before_b}
        if path_ is not None and os.path.basename(path_) not in (_FILES + _ALLOW):
            r = orig(name)             # Tkinter.n・Canvas.n・KBFrame.n など GUI は、組み込みと同名の定義 (Window・Canvas・Options …) が本物なので捨てない
            for k in before_b:         # ListPlot.n の読み込みの途中で GUI が読まれても、外側の復元でその定義を消さないよう覚えておく
                if k in I.downvalues and I.downvalues[k] != before_d.get(k):
                    I._pkg_keep.add(k)
            return r
        saved_keep = I._pkg_keep
        I._pkg_keep = set()
        try:
            orig(name)
            keep = I._pkg_keep
        finally:
            I._pkg_keep = saved_keep
        for k in before_b:
            if k in keep:
                continue
            if k in ('Set', 'SetDelayed', 'Plus', 'Times', 'Power', 'Minus', 'Divide', 'Subtract') and k in I.downvalues:
                old_ = before_d.get(k, [])          # `(a_ = Bump[…]) ^:= …` (Set)、`Plus[Rational[…], …] ^:= …` (rational.n) などの UpValue は組み込みの定義として持ち込まれる。捨てずに加える
                I.downvalues[k] = [d for d in I.downvalues[k] if not any(d is o for o in old_)] + old_
                continue
            if k in I.downvalues:
                if k in before_d:
                    I.downvalues[k] = before_d[k]
                else:
                    del I.downvalues[k]
    I._load_autoload = load
