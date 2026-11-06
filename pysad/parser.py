"""
parser.py - SADScript (.n file) interpreter for pysad

SADScript is a Mathematica-like language used by the SAD accelerator design code.
This module implements a practical interpreter sufficient to load SAD lattice files
and run standard beam optics calculations.

Supported subset:
  - Variable assignments:  x = value;
  - Element definitions:   SetElement[name, type, K1->val, L->val, ...]
  - Element parameter set: Element["K1", "QF"] = 1.5;
  - BeamLine definitions:  LINE = BeamLine[QF, D1, QD, D1];
  - Arithmetic expressions: standard operators, Sin, Cos, Sqrt, Abs, ...
  - Comments:  (* ... *)  and  ! line comment
  - Twiss[beamline]
  - Emittance[beamline, energy]

Not supported (no GUI, no device control):
  - Pattern matching rules (^:=)
  - Module[], With[], If[] with complex patterns
  - Tcl/Tk, X display
  - EPICS
"""

import re
import math
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from .elements import Element, Drift, Quadrupole, Bend, Sextupole, Cavity, Monitor, Marker, Multipole, make_element
from .lattice import Lattice, fodo_cell
from .constants import SpeedOfLight, ElectronMass


# ---------------------------------------------------------------------------
# Tokenizer
# ---------------------------------------------------------------------------

class Token:
    """A lexical token from a SADScript file."""
    __slots__ = ('type', 'value', 'line')

    def __init__(self, type_: str, value: str, line: int = 0):
        self.type = type_
        self.value = value
        self.line = line

    def __repr__(self):
        return f"Token({self.type!r}, {self.value!r}, line={self.line})"


_KEYWORDS = {
    'If', 'While', 'Do', 'For', 'Module', 'With', 'Return', 'Break',
    'BeamLine', 'SetElement', 'Element', 'Table',
    'Sin', 'Cos', 'Tan', 'ArcSin', 'ArcCos', 'ArcTan',
    'Sqrt', 'Abs', 'Exp', 'Log', 'Floor', 'Ceiling', 'Round',
    'Max', 'Min', 'Sign', 'Mod',
    'Pi', 'E', 'True', 'False', 'Null',
    'Print', 'Twiss', 'Emittance', 'RadiationIntegrals',
}

# SAD ネイティブ要素タイプキーワード (大文字・小文字両方)
_SAD_ELEMENT_TYPES = {
    'DRIFT', 'QUAD', 'QUADRUPOLE', 'BEND', 'DIPOLE', 'SBEND', 'RBEND',
    'SEXT', 'SEXTUPOLE', 'OCT', 'OCTUPOLE',
    'MONI', 'MONITOR', 'BPM', 'HMON', 'VMON',
    'MARK', 'MARKER',
    'MULT', 'KICKER', 'HKICK', 'VKICK',
    'CAVI', 'CAVITY', 'RFCA', 'LCAV',
    'CORS', 'CORR', 'SOLENOID', 'SOL',
    'DECA', 'DECAPOLE', 'DODECA', 'DODECAPOLE', 'TCAV', 'UND', 'WIG', 'APERT', 'ST', 'STEER',
    'COORD', 'MAP', 'BEAMBEAM', 'PROF', 'INS', 'SPCH', 'TEST',
}

# 単位変換テーブル (SAD 単位 → SI/eV 基本単位)
# (initb1.f:1618-1670 の icUNIT 表。大文字小文字は大文字化して引く)
_SAD_UNITS = {
    'RAD': 1.0, 'MRAD': 1e-3, 'DEG': math.pi / 180.0,
    'M': 1.0, 'CM': 1e-2, 'MM': 1e-3,
    'T': 1.0, 'GAUSS': 1e-4,
    'JOULE': 1.0 / 1.602176634e-19, 'EV': 1.0, 'KEV': 1e3, 'MEV': 1e6, 'GEV': 1e9,
    'V': 1.0, 'KV': 1e3, 'MV': 1e6, 'GV': 1e9,
    'HZ': 1.0, 'KHZ': 1e3, 'MHZ': 1e6, 'GHZ': 1e9,
}

_TOKEN_RE = re.compile(
    r'(?P<COMMENT_PAREN>\(\*.*?\*\))'                        # (* ... *)
    r'|(?P<COMMENT_BANG>![^\n]*)'                            # ! comment to end of line
    r'|(?P<STRING>"[^"]*")'                                  # "string"
    # 数値: 先頭に小数点がある場合も許容 (SAD では .05 のような記法が多い)
    r'|(?P<NUMBER>\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?)'
    r'|(?P<IDENT>[A-Za-z$][A-Za-z0-9$_`]*)'                # identifier / keyword
    r'|(?P<ARROW>->|:>)'                                     # rule arrows
    # ASSIGN: =. は次が数字でない場合のみ (=.05 を誤認識しないため)
    r'|(?P<ASSIGN>:=|=\.(?!\d)|\^:=|\^=)'
    r'|(?P<SAMEQ><=>)'                                       # SameQ operator
    r'|(?P<OP>[+\-*/^=<>!&|,;:@~#])'                       # operators
    r'|(?P<LBRACE>\[)'
    r'|(?P<RBRACE>\])'
    r'|(?P<LPAREN>\()'
    r'|(?P<RPAREN>\))'
    r'|(?P<LBRACKET>\{)'
    r'|(?P<RBRACKET>\})'
    r'|(?P<WS>\s+)',                                         # whitespace
    re.DOTALL
)


def tokenize(source: str) -> List[Token]:
    """Tokenize SADScript source, stripping comments and whitespace."""
    tokens = []
    line = 1
    for m in _TOKEN_RE.finditer(source):
        kind = m.lastgroup
        val = m.group()
        if kind in ('COMMENT_PAREN', 'COMMENT_BANG', 'WS'):
            line += val.count('\n')
            continue
        tok = Token(kind, val, line)
        line += val.count('\n')
        tokens.append(tok)
    return tokens


# ---------------------------------------------------------------------------
# Environment (symbol table)
# ---------------------------------------------------------------------------

class SADEnvironment:
    """Execution environment for SADScript evaluation."""

    # Physical constants available in SAD
    _BUILTINS = {
        'Pi': math.pi,
        'E': math.e,
        'True': True,
        'False': False,
        'Null': None,
        'Infinity': math.inf,
        'NaN': math.nan,
        # SAD physical constants
        'SpeedOfLight': 2.99792458e8,
        'ElectronMass': 0.51099895e6,   # eV/c^2
        'ElectronRadius': 2.8179403227e-15,  # m
        'ProtonMass': 938.27208816e6,    # eV/c^2
        'PlanckHbar': 1.054571817e-34,   # J·s
        # Degree conversion
        'Degree': math.pi / 180.0,
        # Math functions (mapped to Python)
        'Sin': math.sin,
        'Cos': math.cos,
        'Tan': math.tan,
        'ArcSin': math.asin,
        'ArcCos': math.acos,
        'ArcTan': math.atan,
        'Sqrt': math.sqrt,
        'Abs': abs,
        'Exp': math.exp,
        'Log': math.log,
        'Floor': math.floor,
        'Ceiling': math.ceil,
        'Round': round,
        'Sign': lambda x: (1 if x > 0 else (-1 if x < 0 else 0)),
        'Max': max,
        'Min': min,
        'Mod': math.fmod,
    }

    def __init__(self):
        self._vars: Dict[str, Any] = dict(self._BUILTINS)
        # Element parameters: {name -> {param -> value}}
        self._elements: Dict[str, Dict[str, Any]] = {}
        # BeamLine definitions: {name -> [element_names]}
        self._beamlines: Dict[str, List[str]] = {}
        # Built lattices
        self._lattices: Dict[str, 'Lattice'] = {}
        # FFS 状態 (Twiss["BX","Q*"] / LINE["K1","Q*"] のディスパッチに使用)
        self._ffs_state = None

    # ---- Variable access ----

    def get(self, name: str) -> Any:
        return self._vars.get(name)

    def set(self, name: str, value: Any) -> None:
        self._vars[name] = value

    def has(self, name: str) -> bool:
        return name in self._vars

    # ---- Element registry ----

    def define_element(self, name: str, etype: str, params: dict) -> None:
        """Register a new element with type and parameters."""
        # SAD doelem.f:57-60: 同じ型で再定義したときは書いたキーワードだけ上書き (rdkwdl)。別の型なら SAD はエラー、ここでは置き換える
        old = self._elements.get(name)
        if old is not None and old.get('TYPE') == etype.upper():
            old.update(params)
            return
        entry = dict(params)
        entry['TYPE'] = etype.upper()
        self._elements[name] = entry

    def set_element_param(self, name: str, param: str, value: Any) -> None:
        """Set a parameter on an existing element."""
        if name not in self._elements:
            raise KeyError(f"Element '{name}' not defined")
        self._elements[name][param.upper()] = value

    def get_element_param(self, name: str, param: str) -> Any:
        if name not in self._elements:
            return None
        return self._elements[name].get(param.upper())

    def get_element(self, name: str) -> Optional[Dict]:
        return self._elements.get(name)

    def element_names(self) -> List[str]:
        return list(self._elements.keys())

    # ---- BeamLine registry ----

    def define_beamline(self, name: str, element_refs: List) -> None:
        self._beamlines[name] = element_refs

    def get_beamline(self, name: str) -> Optional[List]:
        return self._beamlines.get(name)

    def resolve_beamline(self, name: str) -> Optional[str]:
        """LINE 名を登録名に直す: 書かれたとおり → 大文字 → 大文字小文字を無視 (SAD は書かれたとおりの名前で引くが、
        pysad は FFS 内で名前を大文字にして扱う所があるので、その両方を受ける)。無ければ None"""
        if name in self._beamlines:
            return name
        u = name.upper()
        if u in self._beamlines:
            return u
        for k in self._beamlines:
            if k.upper() == u:
                return k
        return None

    def beamline_names(self) -> List[str]:
        return list(self._beamlines.keys())

    # ---- Build pysad Lattice ----

    def build_lattice(self, name_or_refs) -> 'Lattice':
        """
        BeamLine (名前またはリスト) を pysad Lattice に変換する。

        HARM パラメータ付き RF 空洞は周長 C を計算してから
        FREQ = HARM × c / C を設定する。
        """
        if isinstance(name_or_refs, str):
            refs = self._beamlines.get(self.resolve_beamline(name_or_refs) or name_or_refs, [])
        else:
            refs = name_or_refs

        elements = []
        self._expand_refs(refs, elements)

        pysad_elements = []
        for el_ref in elements:
            if isinstance(el_ref, tuple) and el_ref[0] == 'REV_EL':
                # 逆向き通過: BEND では E1↔E2 を入れ替える
                pysad_el = self._build_element(el_ref[1], reversed_traversal=True)
            else:
                pysad_el = self._build_element(el_ref)
            if pysad_el is not None:
                pysad_elements.append(pysad_el)

        # 周長 C を計算し、HARM を持つ空洞の周波数を設定
        C = sum(el.L for el in pysad_elements if hasattr(el, 'L'))
        if C > 0:
            for el in pysad_elements:
                if isinstance(el, Cavity) and hasattr(el, '_harm') and el._harm:
                    el.FREQ = el._harm * SpeedOfLight / C

        return Lattice(pysad_elements)

    def _expand_refs(self, refs, out: list, depth: int = 0,
                     is_reversed: bool = False) -> None:
        """
        BeamLine 参照を再帰的に展開する。

        is_reversed=True のとき、要素は逆向き通過を意味する。
        逆向き通過の BEND では E1 と E2 を入れ替えるため、
        ('REV_EL', name) タプルとして記録する。
        """
        if depth > 50:
            raise RecursionError("BeamLine expansion too deep")
        for ref in refs:
            if isinstance(ref, str):
                if ref in self._beamlines:
                    bl_refs = self._beamlines[ref]
                    if is_reversed:
                        bl_refs = list(reversed(bl_refs))
                    self._expand_refs(bl_refs, out, depth + 1,
                                      is_reversed=is_reversed)
                elif ref in self._elements:
                    if is_reversed:
                        out.append(('REV_EL', ref))
                    else:
                        out.append(ref)
                # else: unknown ref (skip silently)
            elif isinstance(ref, list):
                self._expand_refs(list(reversed(ref)) if is_reversed else ref, out, depth + 1, is_reversed)
            elif isinstance(ref, tuple) and ref[0] == 'REPEAT':
                n, inner = ref[1], ref[2]
                for _ in range(n):
                    self._expand_refs([inner], out, depth + 1, is_reversed)
            elif isinstance(ref, tuple) and ref[0] == 'REV_EL':
                # 既に逆向きマーク済み → 更に逆にしたら正順に戻る
                inner = ref[1]
                if is_reversed:
                    # 二重逆 = 正順
                    if inner in self._elements:
                        out.append(inner)
                else:
                    out.append(ref)
            elif isinstance(ref, tuple) and ref[0] == 'REVERSE':
                # -NAME: NAME を逆方向に通過
                inner = ref[1]
                flip = not is_reversed
                if isinstance(inner, (list, tuple)):
                    self._expand_refs([inner], out, depth + 1, is_reversed=flip)
                elif inner in self._beamlines:
                    bl_refs = list(reversed(self._beamlines[inner]))
                    self._expand_refs(bl_refs, out, depth + 1,
                                      is_reversed=flip)
                elif inner in self._elements:
                    if flip:
                        out.append(('REV_EL', inner))
                    else:
                        out.append(inner)

    def _build_element(self, name: str,
                       reversed_traversal: bool = False) -> Optional[Element]:
        """
        Build a pysad Element from a registered element definition.

        reversed_traversal=True のとき、BEND 要素では入口/出口端部角 (E1/E2) を
        入れ替える。これは SAD の -BEAMLINE 構文 (逆向き通過) に対応する。
        """
        info = self._elements.get(name)
        if info is None:
            return None

        etype = info.get('TYPE', 'DRIFT').upper()
        params = {k: v for k, v in info.items() if k != 'TYPE'}
        el = self._build_element0(name, etype, params, reversed_traversal)
        if el is not None:
            el._raw = {k: v for k, v in params.items() if isinstance(v, (int, float))}   # FFS の TYPE/LINE 用に元の値
            el._rev = bool(reversed_traversal)       # 逆向き通過 (BEND の E1/E2 属性は入れ替わっている。FFS の LINE["E1"] は定義の側を返す)
            self._apply_extra_params(el, etype, params)
        return el

    # SAD キーワード → pysad 属性 (lattice/emit6d/track6 が参照する名前)。ROTATE は BEND では設計ロール (Bend.ROTATE)、
    # それ以外ではロール (ROTATION)。
    _EXTRA_COMMON = ('DX', 'DY', 'DZ', 'DROTATE', 'K0', 'F1', 'F2', 'FB1', 'FB2', 'FRINGE', 'DISFRIN', 'EPS', 'K0FR',
                     'F1K1F', 'F1K1B', 'F2K1F', 'F2K1B', 'AE1', 'AE2')

    @staticmethod
    def _numv(v):
        return isinstance(v, (int, float)) and not isinstance(v, bool)

    def _apply_extra_params(self, el, etype, params):
        tn = type(el).__name__
        if tn in ('Drift', 'Marker', 'Monitor') and getattr(el, '_sad_type', None) != 'BEND':
            for k in ('DX', 'DY'):
                if self._numv(params.get(k)):
                    setattr(el, k, float(params[k]))
            return
        for k in self._EXTRA_COMMON:
            v = params.get(k)
            if self._numv(v):
                setattr(el, k, float(v))
        rot = params.get('ROTATE')
        if self._numv(rot) and tn != 'Bend':
            el.ROTATION = float(rot)
        if tn == 'Cavity':
            for k in ('PHI', 'DPHI', 'DVOLT'):
                v = params.get(k)
                if self._numv(v):
                    setattr(el, k, float(v))
        elif tn == 'Multipole':
            knl = {}
            ksl = {}
            for k, v in params.items():
                if not self._numv(v):
                    continue
                m = re.fullmatch(r'K(\d+)', k)
                if m:
                    knl[int(m.group(1))] = float(v)
                m = re.fullmatch(r'SK(\d+)', k)
                if m:
                    ksl[int(m.group(1))] = float(v)
            if knl:
                el.KnL = knl
            if ksl:
                el.KsL = ksl
            for k in ('ANGLE', 'E1', 'E2', 'VOLT', 'FREQ', 'PHI', 'HARM', 'BZ'):
                v = params.get(k)
                if self._numv(v):
                    setattr(el, k, float(v))
        elif tn == 'Solenoid':
            for k in ('BZ', 'BOUND', 'GEO'):
                v = params.get(k)
                if self._numv(v):
                    setattr(el, k, float(v))

    def _build_element0(self, name, etype, params, reversed_traversal):

        # Normalize parameter names (SAD -> MAD-style)
        L = float(params.get('L', params.get('LENGTH', 0.0)) or 0.0)
        angle = float(params.get('ANGLE', 0.0) or 0.0)
        K1 = float(params.get('K1', 0.0) or 0.0)
        K2 = float(params.get('K2', 0.0) or 0.0)
        K3 = float(params.get('K3', 0.0) or 0.0)
        K1L = float(params.get('K1L', 0.0) or 0.0)
        K2L = float(params.get('K2L', 0.0) or 0.0)
        E1 = float(params.get('E1', 0.0) or 0.0)
        E2 = float(params.get('E2', 0.0) or 0.0)
        HGAP = float(params.get('HGAP', 0.0) or 0.0)
        FINT = float(params.get('FINT', 0.5) or 0.5)
        # SAD の BEND E1/E2 は常に ANGLE に対する比 (tbende.f: pe1 = psi1*phi0 + apsi1)。
        # AE1/AE2 は絶対角 [rad] を加える。pysad の Bend.E1/E2 は端部角 [rad]。
        AE1 = float(params.get('AE1', 0.0) or 0.0)
        AE2 = float(params.get('AE2', 0.0) or 0.0)
        E1 = E1 * angle + AE1
        E2 = E2 * angle + AE2
        voltage = float(params.get('VOLT', params.get('VOLTAGE', 0.0)) or 0.0)
        freq = float(params.get('FREQ', params.get('FREQUENCY', 0.0)) or 0.0)
        harm = int(float(params.get('HARM', 0) or 0))

        if etype in ('DRIFT', 'DR', 'O'):
            return Drift(name=name, L=L)

        elif etype in ('QUAD', 'Q', 'QF', 'QD', 'QUADRUPOLE'):
            # SAD ファイル中の K1 は積分勾配 [m⁻¹] = K1_local * L
            # pysad 要素の K1 は局所勾配 [m⁻²] → L で割る
            if K1L != 0.0 and L > 0:
                K1 = K1L / L   # K1L は既に局所勾配
            elif L > 0:
                K1 = K1 / L   # SAD積分 → 局所勾配
            return Quadrupole(name=name, L=L, K1=K1)

        elif etype in ('BEND', 'B', 'SBEND', 'RBEND', 'DIPOLE'):
            # BEND の ROTATE は設計ロール (tbende の theta0 = ROT)。ANGLE=0 でも ROTATE≠0 (垂直ステアリング VD/PV/VS 等) の
            # BEND は Drift 化せず Bend (ANGLE=0) として保持し、K0/K1 を後から与えても効くようにする (docs/history/work_rotate.md)。
            rotate = float(params.get('ROTATE', 0.0) or 0.0)
            # SAD ファイル中の K1 (BEND) は積分勾配 [m⁻¹] → 局所勾配 [m⁻²] へ変換
            if L > 0 and abs(K1) > 0:
                K1 = K1 / L
            # 逆向き通過では入口/出口端部角を入れ替える
            E1_use = E2 if reversed_traversal else E1
            E2_use = E1 if reversed_traversal else E2
            return Bend(name=name, L=L, ANGLE=angle, K1=K1,
                        E1=E1_use, E2=E2_use, HGAP=HGAP, FINT=FINT,
                        ROTATE=rotate)

        elif etype in ('SEXT', 'S', 'SEXTUPOLE'):
            # SAD ファイル中の K2 は積分値 [m⁻²]。pysad の Sextupole.K2 は長さあたり → L で割る (L=0 ならそのまま)
            if K2L != 0.0 and L > 0:
                K2 = K2L / L
            elif L > 0:
                K2 = K2 / L
            return Sextupole(name=name, L=L, K2=K2)

        elif etype in ('OCT', 'OCTUPOLE'):
            # SAD の OCT は tthin (nord=8)。K3 は積分値 → pysad の Octupole.K3 は長さあたり
            from .elements import Octupole
            return Octupole(name=name, L=L, K3=(K3 / L if L > 0 else K3))

        elif etype in ('DECA', 'DECAPOLE'):
            from .elements import Decapole
            K4 = float(params.get('K4', 0.0) or 0.0)
            return Decapole(name=name, L=L, K4=(K4 / L if L > 0 else K4))

        elif etype in ('DODECA', 'DODECAPOLE'):
            from .elements import Dodecapole
            K5 = float(params.get('K5', 0.0) or 0.0)
            return Dodecapole(name=name, L=L, K5=(K5 / L if L > 0 else K5))

        elif etype in ('SOL', 'SOLENOID'):
            from .elements import Solenoid
            return Solenoid(name=name, L=L)

        elif etype in ('CAVI', 'CAVITY', 'RFCA', 'LCAV'):
            # HARM: 高調波数 → FREQ = HARM * c / C (C = 周長)
            # C が未知の場合は HARM を属性として保存し後で設定
            cav = Cavity(name=name, L=L, VOLT=voltage, FREQ=freq)
            if harm and not freq:
                # 周長が判明次第 FREQ を更新するためにメモしておく
                cav._harm = harm  # type: ignore[attr-defined]
            return cav

        elif etype in ('MONI', 'MONITOR', 'BPM', 'HMON', 'VMON'):
            return Monitor(name=name, L=L)

        elif etype in ('MARK', 'MARKER'):
            return Marker(name=name)

        elif etype in ('MULT', 'KICKER', 'HKICK', 'VKICK', 'CORS', 'CORR'):
            return Multipole(name=name, L=L)

        else:
            # Unknown type – treat as drift
            return Drift(name=name, L=L)

    def summary(self) -> str:
        lines = [
            f"SADEnvironment: {len(self._elements)} elements, "
            f"{len(self._beamlines)} beamlines, "
            f"{len(self._vars) - len(self._BUILTINS)} user variables"
        ]
        if self._elements:
            lines.append("Elements:")
            for n, d in self._elements.items():
                lines.append(f"  {n}: {d.get('TYPE','?')}  {d}")
        if self._beamlines:
            lines.append("BeamLines:")
            for n, refs in self._beamlines.items():
                lines.append(f"  {n}: {refs[:8]}{'...' if len(refs)>8 else ''}")
        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Expression evaluator (recursive descent)
# ---------------------------------------------------------------------------

class _Parser:
    """
    Recursive-descent parser/evaluator for a practical subset of SADScript.

    Grammar (simplified):
      program   := statement*
      statement := assignment ';'
                 | expr ';'
      assignment := ident ('=' | ':=') expr
                  | 'Element' '[' string ',' ident ']' '=' expr
                  | 'SetElement' '[' ident ',' string ',' rules ']'
                  | ident '=' 'BeamLine' '[' refs ']'
      expr      := comparison
      comparison := add (('==' | '<' | '>' | '<=>' | '<=' | '>=') add)*
      add       := mul (('+' | '-') mul)*
      mul       := power (('*' | '/') power)*
      power     := unary ('^' unary)*
      unary     := ('-' | '+') primary | primary
      primary   := number | string | ident | ident '[' args ']'
                 | '{' list '}'  | '(' expr ')'
    """

    def __init__(self, tokens: List[Token], env: SADEnvironment):
        self.tokens = tokens
        self.pos = 0
        self.env = env

    # ---- Token navigation ----

    def peek(self, offset: int = 0) -> Optional[Token]:
        i = self.pos + offset
        return self.tokens[i] if i < len(self.tokens) else None

    def consume(self) -> Token:
        tok = self.tokens[self.pos]
        self.pos += 1
        return tok

    def expect(self, type_: str, value: str = None) -> Token:
        tok = self.peek()
        if tok is None:
            raise SyntaxError(f"Expected {type_!r} but reached end of input")
        if tok.type != type_:
            raise SyntaxError(f"Line {tok.line}: expected {type_!r}, got {tok.type!r} ({tok.value!r})")
        if value is not None and tok.value != value:
            raise SyntaxError(f"Line {tok.line}: expected {value!r}, got {tok.value!r}")
        return self.consume()

    def match(self, type_: str, value: str = None) -> bool:
        tok = self.peek()
        if tok is None:
            return False
        if tok.type != type_:
            return False
        if value is not None and tok.value != value:
            return False
        self.consume()
        return True

    def match_val(self, *values) -> Optional[str]:
        tok = self.peek()
        if tok and tok.value in values:
            self.consume()
            return tok.value
        return None

    # ---- Top-level ----

    def parse_program(self) -> List[Any]:
        results = []
        while self.peek() is not None:
            try:
                r = self.parse_statement()
                if r is not None:
                    results.append(r)
            except (SyntaxError, KeyError, ValueError, TypeError) as e:
                # Skip to next semicolon on error
                tok = self.peek()
                line = tok.line if tok else '?'
                print(f"[pysad parser] Warning line {line}: {e}")
                self._skip_to_semicolon()
        return results

    def _skip_to_semicolon(self):
        depth = 0
        while self.peek() is not None:
            tok = self.peek()
            if tok.type in ('LBRACE', 'LPAREN', 'LBRACKET'):
                depth += 1
            elif tok.type in ('RBRACE', 'RPAREN', 'RBRACKET'):
                depth -= 1
            elif tok.type == 'OP' and tok.value == ';' and depth == 0:
                self.consume()
                return
            self.consume()

    def parse_statement(self) -> Any:
        """Parse one statement (ends with ;)."""
        # Skip standalone semicolons
        while self.match('OP', ';'):
            pass

        tok = self.peek()
        if tok is None:
            return None

        # MAIN レベルのコマンド (TYPE Q*; PRINT ...; 等) は読み飛ばす
        if tok.type == 'IDENT' and tok.value.upper() in ('TYPE', 'PRINT', 'SHOW', 'DISP'):
            nxt = self.peek(1)
            if nxt is None or not (nxt.type == 'OP' and nxt.value in ('=', ':=')) and nxt.type not in ('LBRACE',):
                self._skip_to_semicolon()
                return None

        # Try to detect special forms before general expression
        result = self._try_special_form()
        if result is not None:
            self.match('OP', ';')
            return result

        # General expression / assignment
        result = self.parse_expr()

        # Check for assignment
        tok = self.peek()
        if tok and tok.type in ('OP', 'ASSIGN') and tok.value in ('=', ':=', '^=', '^:='):
            self.consume()
            rhs = self.parse_expr()
            if isinstance(result, str):
                self.env.set(result, rhs)
                result = rhs
            elif isinstance(result, tuple) and result[0] == 'ELEMENT_ACCESS':
                # Element["PARAM", name] = value
                _, param, name = result
                self.env.set_element_param(name, param, rhs)
                result = rhs
            elif isinstance(result, tuple) and result[0] == 'LINE_QUERY':
                # LINE["PARAM", elem] = value  →  FFSState.LINE(param, elem, value)
                _, param, elem = result
                ffs_state = getattr(self.env, '_ffs_state', None)
                if ffs_state is not None:
                    try:
                        ffs_state.LINE(param, elem, rhs)
                    except Exception:
                        pass
                result = rhs
        else:
            # 代入なし → LINE_QUERY を get として評価
            if isinstance(result, tuple) and result[0] == 'LINE_QUERY':
                _, param, elem = result
                ffs_state = getattr(self.env, '_ffs_state', None)
                if ffs_state is not None:
                    try:
                        result = ffs_state.LINE(param, elem)
                    except Exception:
                        result = None

        self.match('OP', ';')
        return result

    def _try_special_form(self) -> Optional[Any]:
        """Try to parse known SADScript special forms."""
        tok = self.peek()
        if tok is None:
            return None

        # SetElement[name, type, K1->v, L->v, ...]
        if tok.type == 'IDENT' and tok.value == 'SetElement':
            return self._parse_set_element()

        # SAD ネイティブ: ETYPE name1=(params) name2=(params) ...
        if tok.type == 'IDENT' and tok.value.upper() in _SAD_ELEMENT_TYPES:
            return self._parse_native_element_block()

        # SAD ネイティブ: LINE name = (ref1 ref2 ...)
        # ただし LINE["param", ...] は FFS 関数呼び出しなのでスキップ
        if tok.type == 'IDENT' and tok.value.upper() == 'LINE':
            next_tok = self.peek(1)
            if next_tok and next_tok.type == 'LBRACE':
                pass  # LINE["param", elem] → _parse_func_call で処理
            else:
                return self._parse_line_block()

        # ident = BeamLine[...]  (Mathematica style)
        if tok.type == 'IDENT':
            next2 = self.peek(1)
            if next2 and next2.type == 'OP' and next2.value == '=':
                next3 = self.peek(2)
                if next3 and next3.type == 'IDENT' and next3.value == 'BeamLine':
                    return self._parse_beamline_assign()

        # Element["PARAM", name] = value  (handled in parse_statement via ELEMENT_ACCESS)
        # Twiss[...] — 2 引数 FFS 形式は _parse_calc_call 内で処理
        if tok.type == 'IDENT' and tok.value in ('Twiss', 'Emittance', 'RadiationIntegrals'):
            return self._parse_calc_call()

        return None

    # ------------------------------------------------------------------
    # SAD ネイティブ構文: ETYPE name1=(params) name2=(params) ...;
    # ------------------------------------------------------------------

    def _parse_native_element_block(self):
        """
        SAD ネイティブ要素定義ブロック。

        例:
          DRIFT DRID1 =(L=2.664) DR05=(L=.05) DR40=(L=.4);
          QUAD  QF1=(L=.42 K1=0.667) QD2H=(L=.16 K1=-0.262);
          BEND  B1=(L=1.0 ANGLE=0.1 E1=0.05);
        """
        etype_tok = self.consume()        # DRIFT / QUAD / BEND etc.
        etype = etype_tok.value.upper()
        defined = []

        while True:
            tok = self.peek()
            # ブロック終了条件
            if tok is None:
                break
            if tok.type == 'OP' and tok.value == ';':
                break
            # 次の要素タイプキーワードが来たら終了 (セミコロンなしで連続定義)
            # ただし、IDENT の次が '=' ならそれは要素「名」なので終了しない
            # (例: "CAVI CAVI=(…)" のように要素名がタイプ名と同じ場合)
            if tok.type == 'IDENT' and tok.value.upper() in _SAD_ELEMENT_TYPES:
                next_tok = self.peek(1)
                if next_tok is None or next_tok.value != '=':
                    break  # 新しいタイプキーワード → このブロック終了
                # else: 次が '=' → 要素名として扱い続行
            # 'LINE' キーワードが来たら終了
            if tok.type == 'IDENT' and tok.value.upper() == 'LINE':
                break
            # 名前以外の OP が来たら終了
            if tok.type == 'OP' and tok.value not in (',',):
                break

            # カンマをスキップ
            if tok.type == 'OP' and tok.value == ',':
                self.consume()
                continue

            # IDENT が要素名
            if tok.type != 'IDENT':
                self.consume()
                continue

            name = self.consume().value

            # オプションの '='
            if not self.match('OP', '='):
                # '=' がない → 要素名だけの参照 (スキップ)
                continue

            # '(' ... ')' の中にパラメータ
            if self.peek() and self.peek().type == 'LPAREN':
                self.consume()  # '('
                params = self._parse_sad_params()
                if self.peek() and self.peek().type == 'RPAREN':
                    self.consume()  # ')'
            else:
                params = {}

            self.env.define_element(name, etype, params)
            defined.append(name)

        return ('ELEMENT_DEF_BLOCK', etype, defined)

    def _parse_sad_params(self) -> dict:
        """
        SAD ネイティブ要素パラメータ '(key=val key=val ...)' の中身をパース。

        キーと値は空白で区切られ、コンマは不要。
        単位サフィックス (DEG, GEV 等) を処理する。
        """
        params = {}
        while True:
            tok = self.peek()
            if tok is None or tok.type == 'RPAREN':
                break
            # セミコロンで終了
            if tok.type == 'OP' and tok.value == ';':
                break
            # カンマをスキップ
            if tok.type == 'OP' and tok.value == ',':
                self.consume()
                continue
            # キー = 値
            if tok.type == 'IDENT':
                key = self.consume().value.upper()
                if not self.match('OP', '='):
                    # 値なし → キーだけ (例: ROTATE)
                    params[key] = True
                    continue
                val = self.parse_expr()
                # 単位サフィックス処理
                unit_tok = self.peek()
                if unit_tok and unit_tok.type == 'IDENT':
                    unit = unit_tok.value.upper()
                    if unit in _SAD_UNITS:
                        val = float(val) * _SAD_UNITS[unit]
                        self.consume()
                params[key] = val
            else:
                # 予期しないトークンはスキップ
                self.consume()
        return params

    # ------------------------------------------------------------------
    # SAD ネイティブ構文: LINE name = (ref1 ref2 ref3 ...)
    # ------------------------------------------------------------------

    def _parse_line_block(self):
        """
        SAD ネイティブ LINE ブロック。複数の LINE を一つのブロックで定義できる。

        例:
          LINE ARCAH = (QF D1 QD D2)
               RING  = (10*ARCAH RFC 10*ARCAH)
          ;
        """
        self.consume()   # 'LINE' トークンを消費
        defined = []

        while True:
            tok = self.peek()
            if tok is None or (tok.type == 'OP' and tok.value == ';'):
                break
            # 次の要素タイプキーワードや 'LINE' が来たら終了 (ただし `TEST = (...)` のように型名と同じ名前の LINE 定義は続き)
            if tok.type == 'IDENT' and (tok.value.upper() in _SAD_ELEMENT_TYPES
                                         or tok.value.upper() == 'LINE'):
                nxt = self.peek(1)
                if not (nxt is not None and nxt.type == 'OP' and nxt.value == '='):
                    break
            if tok.type == 'OP' and tok.value == ',':
                self.consume()
                continue
            if tok.type != 'IDENT':
                break

            name = self.consume().value
            if not self.match('OP', '='):
                break
            if not (self.peek() and self.peek().type == 'LPAREN'):
                break
            self.consume()   # '('
            refs = self._parse_beamline_refs_sad()
            if self.peek() and self.peek().type == 'RPAREN':
                self.consume()   # ')'

            self.env.define_beamline(name, refs)
            self.env.set(name, ('BEAMLINE', name))
            defined.append(name)

        return ('LINE_DEF_BLOCK', defined)

    def _parse_set_element(self):
        """SetElement[name, type, Rule1, Rule2, ...]"""
        self.expect('IDENT', 'SetElement')
        self.expect('LBRACE')
        name_tok = self.consume()
        name = name_tok.value
        self.expect('OP', ',')
        type_tok = self.consume()
        etype = type_tok.value

        params = {}
        while self.match('OP', ','):
            # key -> value  or  key :> value
            key_tok = self.consume()
            key = key_tok.value
            if not self.match('ARROW'):
                # Try bare param=value fallback
                if self.match('OP', '='):
                    pass
                else:
                    continue
            val = self.parse_expr()
            params[key.upper()] = val

        self.expect('RBRACE')
        self.env.define_element(name, etype, params)
        return ('ELEMENT_DEF', name)

    def _parse_beamline_assign(self):
        """name = BeamLine[ref1, ref2, ...]"""
        name = self.consume().value  # ident
        self.consume()              # '='
        self.expect('IDENT', 'BeamLine')
        self.expect('LBRACE')
        refs = self._parse_beamline_refs()
        self.expect('RBRACE')
        self.env.define_beamline(name, refs)
        # Also set as a variable
        self.env.set(name, ('BEAMLINE', name))
        return ('BEAMLINE_DEF', name)

    def _parse_beamline_refs(self) -> List:
        """Parse comma-separated element/beamline references inside BeamLine[...]."""
        refs = []
        while True:
            tok = self.peek()
            if tok is None or (tok.type == 'RBRACE'):
                break

            # n * element  or  -n * element  (repetition)
            if tok.type == 'NUMBER':
                n = int(float(tok.value))
                self.consume()
                if self.match('OP', '*'):
                    inner_tok = self.peek()
                    if inner_tok and inner_tok.type == 'IDENT':
                        inner_name = self.consume().value
                        refs.append(('REPEAT', n, inner_name))
                    else:
                        refs.append(n)
                else:
                    refs.append(n)
            elif tok.type == 'OP' and tok.value == '-':
                # -name  means reversed element (ignore reversal for now)
                self.consume()
                inner_tok = self.consume()
                refs.append(inner_tok.value)
            elif tok.type == 'IDENT':
                refs.append(tok.value)
                self.consume()
            else:
                self.consume()

            if not self.match('OP', ','):
                break
        return refs

    def _parse_beamline_refs_sad(self) -> List:
        """
        SAD ネイティブ LINE 構文内の要素参照をパース。
        コンマなし空白区切り。 ')' で終了。

        例: QF D1 QD D2  または  10*ARCAH RFC 10*ARCAH
            -ARCAH  (逆順 → ('REVERSE', 'ARCAH') として記録)
        """
        refs = []
        while True:
            tok = self.peek()
            if tok is None or tok.type == 'RPAREN':
                break
            if tok.type == 'OP' and tok.value == ';':
                break

            # コンマをスキップ
            if tok.type == 'OP' and tok.value == ',':
                self.consume()
                continue

            # [-] n * [-] name / -name / (…) の繰り返しと反転 (expln.f: 個数 ilist(1,i) の符号が負なら逆順、
            # 要素の向き dir = sign(direct*count) で BEND の E1↔E2 を入れ替える)
            if (tok.type == 'OP' and tok.value in ('-', '+')) or tok.type == 'NUMBER' or tok.type == 'LPAREN':
                neg = False
                if tok.type == 'OP':
                    neg = tok.value == '-'
                    self.consume()
                    tok = self.peek()
                n = 1
                if tok is not None and tok.type == 'NUMBER':
                    n = int(float(tok.value))
                    self.consume()
                    if not self.match('OP', '*'):
                        continue   # 裸の数
                    tok = self.peek()
                    if tok is not None and tok.type == 'OP' and tok.value in ('-', '+'):
                        neg = neg != (tok.value == '-')
                        self.consume()
                        tok = self.peek()
                if tok is None:
                    break
                if tok.type == 'IDENT':
                    inner = self.consume().value
                elif tok.type == 'LPAREN':
                    self.consume()
                    inner = self._parse_beamline_refs_sad()
                    self.match('RPAREN')
                else:
                    self.consume()
                    continue
                if neg:
                    inner = ('REVERSE', inner)
                if n == 1:
                    refs.append(inner)
                else:
                    refs.append(('REPEAT', n, inner))
                continue

            # 通常の要素/ビームライン名
            if tok.type == 'IDENT':
                refs.append(tok.value)
                self.consume()
                continue

            # 予期しないトークンはスキップ
            self.consume()
        return refs

    def _parse_calc_call(self):
        """Handle Twiss[beamline], Emittance[beamline, energy], etc."""
        func_tok = self.consume()
        func = func_tok.value
        self.expect('LBRACE')
        args = []
        while self.peek() and self.peek().type != 'RBRACE':
            args.append(self.parse_expr())
            if not self.match('OP', ','):
                break
        self.expect('RBRACE')

        if func == 'Twiss' and args:
            # FFS スタイル: Twiss["BX", element]  (2 引数)
            if len(args) >= 2:
                param = str(args[0])
                # 2 番目引数: 数値はそのまま、文字列はそのまま保持
                elem_raw = args[1]
                if isinstance(elem_raw, (int, float)):
                    elem = elem_raw          # 数値インデックス
                else:
                    elem = str(elem_raw)     # 要素名 / パターン
                ffs_state = getattr(self.env, '_ffs_state', None)
                if ffs_state is not None:
                    try:
                        return ffs_state.Twiss(param, elem)
                    except Exception:
                        pass
                return ('TWISS_QUERY', param, elem)
            # 従来スタイル: Twiss[beamline]  (1 引数)
            bl_ref = args[0]
            if isinstance(bl_ref, tuple) and bl_ref[0] == 'BEAMLINE':
                lat = self.env.build_lattice(bl_ref[1])
                return lat.twiss()
            elif isinstance(bl_ref, str) and bl_ref in self.env._beamlines:
                lat = self.env.build_lattice(bl_ref)
                return lat.twiss()

        elif func == 'Emittance' and len(args) >= 2:
            bl_ref, energy = args[0], float(args[1])
            if isinstance(bl_ref, tuple) and bl_ref[0] == 'BEAMLINE':
                lat = self.env.build_lattice(bl_ref[1])
            elif isinstance(bl_ref, str):
                lat = self.env.build_lattice(bl_ref)
            else:
                return None
            return lat.emittance(energy)

        elif func == 'RadiationIntegrals' and args:
            bl_ref = args[0]
            if isinstance(bl_ref, tuple) and bl_ref[0] == 'BEAMLINE':
                lat = self.env.build_lattice(bl_ref[1])
            elif isinstance(bl_ref, str):
                lat = self.env.build_lattice(bl_ref)
            else:
                return None
            return lat.radiation_integrals()

        return None

    # ---- Expression parsing ----

    def parse_expr(self) -> Any:
        return self._parse_add()

    def _parse_add(self) -> Any:
        left = self._parse_mul()
        while True:
            op = self.match_val('+', '-')
            if op is None:
                break
            right = self._parse_mul()
            left = left + right if op == '+' else left - right
        return left

    def _parse_mul(self) -> Any:
        left = self._parse_power()
        while True:
            op = self.match_val('*', '/')
            if op is None:
                break
            right = self._parse_power()
            if op == '*':
                left = left * right
            else:
                left = left / right if right != 0 else float('nan')
        return left

    def _parse_power(self) -> Any:
        base = self._parse_unary()
        if self.match_val('^'):
            exp = self._parse_unary()
            return base ** exp
        return base

    def _parse_unary(self) -> Any:
        tok = self.peek()
        if tok and tok.type == 'OP' and tok.value == '-':
            self.consume()
            return -self._parse_primary()
        if tok and tok.type == 'OP' and tok.value == '+':
            self.consume()
        return self._parse_primary()

    def _parse_primary(self) -> Any:
        tok = self.peek()
        if tok is None:
            return 0.0

        # Number literal (SAD では数値の直後に単位サフィックスが来ることがある)
        if tok.type == 'NUMBER':
            self.consume()
            v = tok.value
            val = int(v) if '.' not in v and 'e' not in v.lower() else float(v)
            # 単位サフィックスのチェック (例: 3GEV, 500MEV, 90DEG)
            unit_tok = self.peek()
            if unit_tok and unit_tok.type == 'IDENT':
                unit = unit_tok.value.upper()
                if unit in _SAD_UNITS:
                    val = float(val) * _SAD_UNITS[unit]
                    self.consume()
            return val

        # String literal
        if tok.type == 'STRING':
            self.consume()
            return tok.value[1:-1]  # strip quotes

        # List literal {a, b, c}
        if tok.type == 'LBRACKET':
            return self._parse_list()

        # Parenthesized expression
        if tok.type == 'LPAREN':
            self.consume()
            val = self.parse_expr()
            self.match('RPAREN')
            return val

        # Identifier or function call
        if tok.type == 'IDENT':
            name = tok.value
            self.consume()

            # Function call: name[args]
            if self.peek() and self.peek().type == 'LBRACE':
                return self._parse_func_call(name)

            # Plain identifier -> look up in environment
            val = self.env.get(name)
            if val is not None:
                return val
            # Unknown symbol -> return as string reference
            return name

        # Skip unknown token
        self.consume()
        return 0.0

    def _parse_list(self) -> list:
        self.expect('LBRACKET')
        items = []
        while self.peek() and self.peek().type != 'RBRACKET':
            items.append(self.parse_expr())
            if not self.match('OP', ','):
                break
        self.match('RBRACKET')
        return items

    def _parse_func_call(self, name: str) -> Any:
        self.expect('LBRACE')
        args = []
        while self.peek() and self.peek().type != 'RBRACE':
            # Handle key -> value rule syntax inside func args
            arg = self.parse_expr()
            tok = self.peek()
            if tok and tok.type == 'ARROW':
                self.consume()
                val = self.parse_expr()
                args.append((arg, val))
            else:
                args.append(arg)
            if not self.match('OP', ','):
                break
        self.match('RBRACE')

        # Evaluate known math functions
        fn = self.env.get(name)
        if callable(fn):
            numeric_args = [a for a in args if not isinstance(a, tuple)]
            try:
                return fn(*numeric_args)
            except Exception:
                return None

        # Element["PARAM", element_name]  -> ELEMENT_ACCESS node
        if name == 'Element' and len(args) == 2:
            param = args[0] if isinstance(args[0], str) else str(args[0])
            el_name = args[1] if isinstance(args[1], str) else str(args[1])
            val = self.env.get_element_param(el_name, param)
            if val is not None:
                return val
            return ('ELEMENT_ACCESS', param, el_name)

        # Table[expr, {i, n}]  -> list
        if name == 'Table' and len(args) == 2:
            return None  # Not critical for physics

        # If[cond, t, f]
        if name == 'If' and len(args) >= 2:
            cond = args[0]
            return args[1] if cond else (args[2] if len(args) > 2 else None)

        # Print[...]
        if name == 'Print':
            print(*[a for a in args if not isinstance(a, tuple)])
            return None

        # ------------------------------------------------------------------
        # FFS コマンド文字列: FFS["cell calc"] / FFS["fit NX 9.5", 3]
        # ------------------------------------------------------------------
        if name == 'FFS' and len(args) >= 1:
            cmd = str(args[0]) if not isinstance(args[0], str) else args[0]
            nsteps = int(args[1]) if len(args) >= 2 else 6
            ffs_state = getattr(self.env, '_ffs_state', None)
            if ffs_state is not None:
                try:
                    return _dispatch_ffs(ffs_state, cmd, nsteps)
                except Exception as e:
                    if self.env.get('__verbose__'):
                        print(f"  [FFS] {e}")
            return None

        # ------------------------------------------------------------------
        # FFS スタイル: Twiss["PARAM", element] → FFSState.Twiss()
        # ------------------------------------------------------------------
        if name == 'Twiss' and len(args) >= 2:
            param = str(args[0]) if not isinstance(args[0], str) else args[0]
            elem  = str(args[1]) if not isinstance(args[1], str) else args[1]
            ffs_state = getattr(self.env, '_ffs_state', None)
            if ffs_state is not None:
                try:
                    return ffs_state.Twiss(param, elem)
                except Exception:
                    pass
            # FFSState が未設定なら (TWISS_QUERY, param, elem) ノードを返す
            return ('TWISS_QUERY', param, elem)

        # ------------------------------------------------------------------
        # FFS スタイル: LINE["PARAM", element] → ('LINE_QUERY', param, elem)
        # 代入検出のため常にタプルを返す。get/set は parse_statement で判断。
        # ------------------------------------------------------------------
        if name == 'LINE' and len(args) >= 1:
            param = str(args[0]) if not isinstance(args[0], str) else args[0]
            elem_raw = args[1] if len(args) >= 2 else None
            if isinstance(elem_raw, (int, float)):
                elem = elem_raw
            elif elem_raw is not None:
                elem = str(elem_raw)
            else:
                elem = None
            return ('LINE_QUERY', param, elem)

        return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class SADScript:
    """
    SADScript interpreter.

    Load and evaluate SAD .n lattice files, then build pysad Lattice objects.

    Usage
    -----
    >>> interp = SADScript()
    >>> interp.load_file("atf.n")
    >>> lat = interp.build_lattice("RING")
    >>> tw = lat.twiss()
    """

    def __init__(self, verbose: bool = False):
        self.env = SADEnvironment()
        self.verbose = verbose

    def load_file(self, path: str) -> None:
        """Parse and evaluate a SADScript .n file."""
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            source = f.read()
        self.load_string(source)

    def load_string(self, source: str) -> None:
        """Parse and evaluate SADScript source text."""
        tokens = tokenize(source)
        parser = _Parser(tokens, self.env)
        results = parser.parse_program()
        if self.verbose:
            print(f"[pysad] Parsed {len(results)} statements")

    def set_variable(self, name: str, value: Any) -> None:
        """Set a variable in the interpreter environment."""
        self.env.set(name, value)

    def get_variable(self, name: str) -> Any:
        """Get a variable value from the interpreter environment."""
        return self.env.get(name)

    def define_element(self, name: str, etype: str, **params) -> None:
        """Manually define an element (alternative to loading from file)."""
        self.env.define_element(name, etype, params)

    def define_beamline(self, name: str, refs: List[str]) -> None:
        """Manually define a beamline."""
        self.env.define_beamline(name, refs)

    def build_lattice(self, name: str) -> 'Lattice':
        """Build a pysad Lattice from a BeamLine defined in the environment."""
        return self.env.build_lattice(name)

    def eval(self, expr: str) -> Any:
        """Evaluate a single SADScript expression."""
        tokens = tokenize(expr + ';')
        parser = _Parser(tokens, self.env)
        return parser.parse_statement()

    def element_names(self) -> List[str]:
        return self.env.element_names()

    def beamline_names(self) -> List[str]:
        return self.env.beamline_names()

    def summary(self) -> str:
        return self.env.summary()

    # ------------------------------------------------------------------
    # FFS インターフェース
    # ------------------------------------------------------------------

    @property
    def ffs(self) -> 'FFSState':  # type: ignore[name-defined]
        """
        FFS 状態オブジェクトを返す。

        Returns
        -------
        FFSState
            このインタープリタの環境に紐付いた FFSState。

        Examples
        --------
        >>> interp = load_sad_file("pf_ring.n")
        >>> f = interp.ffs
        >>> f.use("RING").cell_calc()
        >>> print(f.Twiss("BX", "$$$"))
        """
        from .ffs import FFSState
        if not hasattr(self, '_ffs_state'):
            self._ffs_state = FFSState(self.env)
            # env にも紐付けて Twiss["BX","Q*"] 等のパーサーから参照できるようにする
            self.env._ffs_state = self._ffs_state
        return self._ffs_state

    def run_ffs(self, cmd: str, nsteps: int = 6) -> Any:
        """
        SAD FFS コマンド文字列を実行する。

        SAD スクリプト内の ``FFS["cell calc"]`` 等を Python から
        直接呼び出すためのインターフェース。

        Parameters
        ----------
        cmd : str
            FFS コマンド文字列 (例: ``"cell calc"``, ``"fit NX 9.5"``).
        nsteps : int
            ``go`` コマンドの反復回数。

        Returns
        -------
        Any
            コマンドに応じた結果。クエリ系は値を返す。

        Examples
        --------
        >>> interp.run_ffs("cell calc emit")
        >>> interp.run_ffs("fit NX 9.5975 3")
        >>> interp.run_ffs("free Q* K1")
        >>> interp.run_ffs("go", 6)
        """
        return _dispatch_ffs(self.ffs, cmd, nsteps)


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

def load_sad_file(path: str, verbose: bool = False) -> SADScript:
    """
    Load a SAD .n lattice file and return a SADScript interpreter ready to use.

    Parameters
    ----------
    path : str
        Path to the .n file
    verbose : bool
        Print parsing progress

    Returns
    -------
    interp : SADScript
        Interpreter with all definitions loaded

    Example
    -------
    >>> interp = load_sad_file("ring.n")
    >>> lat = interp.build_lattice("RING")
    >>> tw = lat.twiss()
    """
    interp = SADScript(verbose=verbose)
    interp.load_file(path)
    return interp


def load_sad_string(source: str, verbose: bool = False) -> SADScript:
    """
    Parse a SAD script given as a string and return a SADScript interpreter.

    This is the string-input equivalent of :func:`load_sad_file`.  Use it
    when you want to write (or generate) a SAD lattice directly in Python
    without creating a temporary file.

    SAD ネイティブ記法の書き方::

        source = '''
        MOMENTUM=2.5 GEV;

        QUAD QF=(K1= 5.0  L=0.3)
             QD=(K1=-4.8  L=0.3)

        BEND B=(L=1.942696 ANGLE=0.22439947 E1=0.5 E2=0.5)

        DRIFT D1=(L=0.5)
              D2=(L=1.2)

        LINE
          CELL = (QF D1 B D2 QD D2 B D1)
          RING = (14*CELL)
        ;
        '''

        interp = load_sad_string(source)
        lat = interp.build_lattice("RING")

    Parameters
    ----------
    source : str
        SAD script text (same format as a .n file).
    verbose : bool
        Print parsing progress if True.

    Returns
    -------
    interp : SADScript
        Interpreter with all definitions loaded, ready for
        ``build_lattice()``, ``element_names()``, etc.
    """
    interp = SADScript(verbose=verbose)
    interp.load_string(source)
    return interp


def load_sad_files(*paths: str, verbose: bool = False) -> SADScript:
    """
    複数の SAD ファイルを同一インタープリタに順番にロードする。

    SAD では通常、メインファイル（要素・ビームライン定義）と
    コマンドファイル（FFS コマンド: Twiss, Emit 等）が分かれている。
    この関数を使うとその 2 ファイルを分けて渡せる。

    Parameters
    ----------
    *paths : str
        SAD ファイルパス（1 個以上）。先頭から順番にパースされる。
        後のファイルは前のファイルで定義された変数・要素を参照できる。
    verbose : bool
        True にするとパース進捗を表示する。

    Returns
    -------
    interp : SADScript
        全ファイルをロード済みのインタープリタ。

    Example
    -------
    >>> interp = load_sad_files("pf_ring.n", "ffs_commands.n")
    >>> lat = interp.build_lattice("RING")

    コマンドラインから::

        python plot_optics.py pf_ring.n ffs_commands.n
    """
    if not paths:
        raise ValueError("load_sad_files: 少なくとも 1 つのファイルパスが必要です")
    interp = SADScript(verbose=verbose)
    for path in paths:
        interp.load_file(path)
        if verbose:
            print(f"[pysad] ロード完了: {path}")
    return interp


# ---------------------------------------------------------------------------
# FFS コマンド文字列ディスパッチャー
# ---------------------------------------------------------------------------

def _dispatch_ffs(ffs, cmd: str, nsteps: int = 6) -> Any:
    """
    SAD FFS コマンド文字列を FFSState メソッドにディスパッチする。

    対応コマンド
    ------------
    cell [calc] [emit]   — cell_calc() [+ emit_calc()]
    ins [calc]           — ins_calc()
    emit                 — emit_calc()
    disp [element]       — disp(element)
    free PATTERN [PARAM] — free(pattern, param)
    fit ELEMENT PARAM TARGET [NSTEPS]   — fit_add(...)
    fit ELEM1 ELEM2 PARAM TARGET        — fit_add (2要素範囲; 簡易)
    go [N]               — go(n)
    FIT; NX target N;    — 複数行 fit/free/go ブロック
    """
    from .ffs import FFSState  # noqa: F401

    # セミコロン・改行で分割して複数コマンドを処理
    sub_cmds = [c.strip().rstrip(';') for c in re.split(r'[;\n]', cmd)]
    result = None
    for sub in sub_cmds:
        if not sub:
            continue
        result = _dispatch_single_ffs(ffs, sub, nsteps)
    return result


def _dispatch_single_ffs(ffs, cmd: str, nsteps: int = 6) -> Any:
    """FFS 単一コマンドをディスパッチ。"""
    tokens = cmd.strip().split()
    if not tokens:
        return None

    kw = tokens[0].upper()

    # ---------- USE ----------
    if kw == 'USE':
        if len(tokens) >= 2:
            ffs.use(tokens[1].upper())
        return None

    # ---------- cell ----------
    if kw == 'CELL':
        # "cell calc", "cell calc emit"
        rest = [t.upper() for t in tokens[1:]]
        ffs.cell_calc()
        if 'EMIT' in rest:
            try:
                ffs.emit_calc()
            except Exception as e:
                print(f"  [emit_calc] {e}")
        return None

    # ---------- ins ----------
    if kw == 'INS':
        ffs.ins_calc()
        return None

    # ---------- emit ----------
    if kw == 'EMIT':
        try:
            return ffs.emit_calc()
        except Exception as e:
            print(f"  [emit_calc] {e}")
        return None

    # ---------- calc ----------
    if kw == 'CALC':
        ffs.cell_calc()
        return None

    # ---------- disp ----------
    if kw == 'DISP':
        elem = tokens[1] if len(tokens) >= 2 else None
        return ffs.disp(elem, verbose=True)

    # ---------- free ----------
    if kw == 'FREE':
        # free Q*  /  free Q* K1  /  free SD*,SF*
        if len(tokens) >= 2:
            # カンマ区切りパターン対応
            patterns = tokens[1].split(',')
            param = tokens[2].upper() if len(tokens) >= 3 else 'K1'
            # 複数パターン → _free_params を累積する
            ffs._free_params = []
            for pat in patterns:
                pat = pat.strip()
                if not pat:
                    continue
                indices = ffs._resolve_indices(pat)
                for i in indices:
                    ffs._free_params.append((i, ffs._names[i], param))
        return None

    # ---------- fit ----------
    if kw == 'FIT':
        # SAD の fit 構文バリエーション:
        # (A) fit ELEMENT PARAM TARGET [NSTEPS]
        #       fit MKC NX 9.5  3
        # (B) fit ELEM1 ELEM2 PARAM TARGET
        #       fit MKSFA1.1 MKSFA1.2 NX 1.5;
        # (C) FIT; — フィット開始マーカー (go で実行)
        # (D) fit MKS BY 3  — "BY" が PARAM, "3" が TARGET 扱い

        # 引数なし → 無視
        if len(tokens) < 2:
            return None

        # ---- FIT; (セミコロンだけ) → 無視 ----
        second = tokens[1].upper()

        # free / NX / NY などの制約コマンドが来るケース
        # fit NX target nsteps  (element 省略 → $$$)
        _TWISS_PARAMS = {'NX', 'NY', 'BX', 'BY', 'AX', 'AY', 'EX', 'EPX', 'EY', 'EPY'}
        if second in _TWISS_PARAMS:
            # fit NX 9.5975 [3]
            param = second
            target = float(tokens[2]) if len(tokens) >= 3 else 0.0
            ns = int(tokens[3]) if len(tokens) >= 4 else nsteps
            ffs.fit_add('$$$', param, target, ns)
            return None

        # 2 番目がパターン/要素名、3 番目がパラメータかどうか判定
        # (B) fit ELEM1 ELEM2 PARAM TARGET
        if len(tokens) >= 5:
            try:
                float(tokens[4])
                # tokens: fit ELEM1 ELEM2 PARAM TARGET [NSTEPS]
                elem = tokens[2]   # 終端要素を使う
                param = tokens[3].upper()
                target = float(tokens[4])
                ns = int(tokens[5]) if len(tokens) >= 6 else nsteps
                ffs.fit_add(elem, param, target, ns)
                return None
            except (ValueError, IndexError):
                pass

        # (A) fit ELEM PARAM TARGET [NSTEPS]
        if len(tokens) >= 4:
            try:
                float(tokens[3])
                elem  = tokens[1]
                param = tokens[2].upper()
                target = float(tokens[3])
                ns = int(tokens[4]) if len(tokens) >= 5 else nsteps
                ffs.fit_add(elem, param, target, ns)
                return None
            except (ValueError, IndexError):
                pass

        # (D) fit ELEM PARAM TARGET (TARGET は 3 番目)
        if len(tokens) == 3:
            try:
                float(tokens[2])
                elem  = '$$$'
                param = tokens[1].upper()
                target = float(tokens[2])
                ffs.fit_add(elem, param, target, nsteps)
                return None
            except (ValueError, IndexError):
                pass

        return None

    # ---------- go ----------
    if kw == 'GO':
        n = int(tokens[1]) if len(tokens) >= 2 else nsteps
        return ffs.go(n)

    # ---------- fitp ----------
    if kw == 'FITP':
        # FITP 1 → fit clear + 再実行 (SAD の FITP コマンド)
        ffs.fit_clear()
        return None

    # ---------- type ----------
    if kw == 'TYPE':
        # TYPE element — 簡易表示 (disp と同様)
        elem = tokens[1] if len(tokens) >= 2 else None
        return ffs.disp(elem, verbose=True)

    # 未対応コマンドは無視
    return None


# ---------------------------------------------------------------------------
# MAIN レベルのコマンド列 (FFS; の前の ON/OFF, EMIT, TRACK ... STOP)
# ---------------------------------------------------------------------------

# MAIN で `NAME=値;` と書いて大域変数を設定するもの (track.f / tracka.f / tinip が rgetgl1 で読む)。
# PHSPLOTS は描画ファイルの番号なので値だけ保持して使わない
_MAIN_SET_NAMES = {'PHSPLOTS', 'EMITX', 'EMITY', 'EMITZ', 'SIGE', 'SIGZ', 'NPARA', 'LOSSAMPL', 'LOSSDZ', 'FSHIFT', 'ESHIFT',
                   'TDXI', 'TDYI', 'TDZI', 'MINCOUP', 'DP', 'DP0', 'MAXORDER', 'ADDTERMS', 'CHARGE', 'MASS',
                   'EMITDIV', 'EMITDIVB', 'EMITDIVQ', 'EMITDIVS', 'FRINGDIV', 'PBUNCH', 'NBUNCH'}
_MAIN_STOP = {'STOP', 'EXIT', 'END', 'QUIT', 'ABORT'}


def _strip_comments(src: str) -> str:
    """SAD の注釈 `! ... 行末` と `(* ... *)` を落とす (文字列リテラルの中の ! は残す)"""
    src = re.sub(r'\(\*.*?\*\)', ' ', src, flags=re.S)
    out = []
    for ln in src.split('\n'):
        q = False
        for i, c in enumerate(ln):
            if c == '"':
                q = not q
            elif c == '!' and not q:
                ln = ln[:i]
                break
        out.append(ln)
    return '\n'.join(out)


def _split_statements(text: str) -> List[str]:
    """`;` で文に分ける ((), [], {}, "" の中の ; では切らない)"""
    out, buf, depth, q = [], [], 0, False
    for c in text:
        if q:
            buf.append(c)
            if c == '"':
                q = False
            continue
        if c == '"':
            q = True
        elif c in '([{':
            depth += 1
        elif c in ')]}':
            depth = max(0, depth - 1)
        elif c == ';' and depth == 0:
            out.append(''.join(buf).strip())
            buf = []
            continue
        buf.append(c)
    tail = ''.join(buf).strip()
    if tail:
        out.append(tail)
    return [s for s in out if s]


def _eval_main_num(tok: str) -> float:
    """MAIN の数値リテラル/簡単な式 (1, 1E-3, 2*3, 1.5GEV など)。評価できなければ ValueError"""
    t = tok.strip()
    try:
        return float(t)
    except ValueError:
        pass
    m = re.fullmatch(r'([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)([A-Za-z]+)', t)
    if m and m.group(2).upper() in _SAD_UNITS:
        return float(m.group(1)) * _SAD_UNITS[m.group(2).upper()]
    if re.fullmatch(r'[-+*/^().\deE ]+', t):
        try:
            return float(eval(t.replace('^', '**'), {'__builtins__': {}}, {}))
        except Exception:
            pass
    raise ValueError(tok)


def _parse_main_list(body: str) -> List[float]:
    """`(1,10,20)` / `(1 10 20)` / `{1,2}` の中身 → 数値リスト。MAIN のリストは区切りが空白でもコンマでもよい"""
    b = body.strip()
    if b and b[0] in '({[' and b[-1] in ')}]':
        b = b[1:-1]
    return [_eval_main_num(t) for t in re.split(r'[,\s]+', b.strip()) if t]


_TRACK_KV = re.compile(r'([A-Za-z_$][\w$]*)\s*=\s*(\([^)]*\)|\{[^}]*\}|[^\s=(){}]+)')


def parse_track_statement(stmt: str) -> dict:
    """`TRACK USE=RING NPART=21 TURNS=1000 NX=(1,10,20) NY=(0) NZ=(1)` → dict (track.f の引数表: initbl.f の TRACK 行)。

    キーワード (USE TURNS NPART XEMI YEMI SIGS SIGE ELEC POSI PROTON RADI SPAN CENTER AFTER AXIS SYNC NX NY NZ) のうち
    追跡に効くのは USE / NPART / TURNS / NX / NY / NZ と粒子種 (ELEC POSI PROTON)。XEMI YEMI SIGS SIGE は ActTra.f が読むだけで
    使わない (初期分布の大きさは大域変数 EMITX EMITY EMITZ SIGE から tracka.f が決める)。同じキーワードを 2 回書くと後のもの
    (ActTra の引数解析は上書き)。値の評価できないキーワードは無視する"""
    body = re.sub(r'^\s*TRACK\b', '', stmt, flags=re.I)
    spec = {'cmd': 'TRACK', 'use': None, 'npart': None, 'turns': None, 'nx': None, 'ny': None, 'nz': None, 'ptype': None}
    kvs = _TRACK_KV.findall(body)
    for k, v in kvs:
        K = k.upper()
        try:
            if K == 'USE':
                spec['use'] = v.upper()
            elif K in ('NPART', 'TURNS'):
                spec[K.lower()] = int(_eval_main_num(v))
            elif K in ('NX', 'NY', 'NZ'):
                spec[K.lower()] = _parse_main_list(v)
        except ValueError:
            pass
    rest = _TRACK_KV.sub(' ', body)
    for w in re.findall(r'[A-Za-z_$][\w$]*', rest):
        if w.upper() in ('ELEC', 'POSI', 'PROTON'):
            spec['ptype'] = w.upper()
    return spec


def parse_main_commands(src: str) -> List[dict]:
    """MAIN 部 (FFS; より前) の命令を順に拾う。拾うもの:
        ON/OFF フラグ…、EMIT [USE=名前]、TRACK …、`NAME=値` (_MAIN_SET_NAMES)、STOP 類 (これ以降は読まない)
    要素/ビームラインの定義など他の文は SADScript の読み込み (load_sad_string) に任せるのでここでは無視する。"""
    cmds: List[dict] = []
    for stmt in _split_statements(_strip_comments(src)):
        m = re.match(r'[A-Za-z_$][\w$]*', stmt)
        if not m:
            continue
        W = m.group(0).upper()
        rest = stmt[m.end():].strip()
        # `track[...]:=...` (SADScript の関数定義)、`emit = Emittance[]` (代入) などは MAIN の命令ではない
        if W in ('TRACK', 'EMIT') and rest[:1] in ('[', '=', ':', '(', '{'):
            continue
        if W == 'TRACK':
            cmds.append(parse_track_statement(stmt))
        elif W == 'EMIT':
            u = re.search(r'USE\s*=\s*([\w$.]+)', rest, flags=re.I)
            cmds.append({'cmd': 'EMIT', 'use': u.group(1).upper() if u else None})
        elif W in ('ON', 'OFF') and rest and re.fullmatch(r'[A-Za-z_$][\w$]*(?:\s+[A-Za-z_$][\w$]*)*', rest):
            cmds.append({'cmd': W, 'flags': rest.upper().split()})
        elif W in _MAIN_STOP and not rest:
            cmds.append({'cmd': 'STOP'})
            break
        elif W in _MAIN_SET_NAMES and rest.startswith('='):
            try:
                cmds.append({'cmd': 'SET', 'name': W, 'value': _eval_main_num(rest[1:])})
            except ValueError:
                pass
    return cmds
