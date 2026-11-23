# -*- coding: utf-8 -*-
"""SAD の台本 (FFS・SADScript) の簡易な Pygments レキサー。文法の完全な解析ではなく、読みやすさのための色分け。

    コメント  ! から行末まで、(* ... *)
    文字列    "..."
    数        1  1.5  .5  1E-3  1.2D3  (単位つきの数 1.7MV は 数 + 名前)
    キーワード  FFS の命令 (USE・CALC・EMIT・DRAW ...)、要素名の頭 (DRIFT・BEND・QUAD ...)、
               制御 (Module・Do・If・With ...)
    組み込み   Print・Fit・Show・Plot・TrackParticles など (大文字小文字は書いたとおりに一致)
    その他     [ ] { } ( ) , ;  と演算子 (-> := == // など)

どの文字列も色分けのときに落とさない・変えない (トークンをつなげると元の文字列に戻る)。
"""
import re

from pygments.lexer import RegexLexer, bygroups, words
from pygments.token import (Comment, Keyword, Name, Number, Operator, Punctuation, String, Text,
                            Whitespace)

# FFS の命令・要素の種類 (SADHelp にある、文書と例でよく使うもの。大文字で書かれる)
COMMANDS = (
    'FFS', 'USE', 'CALC', 'CALCULATE', 'CELL', 'INS', 'EMIT', 'DRAW', 'DISPLAY', 'SHOW', 'GO', 'FIT', 'FREE', 'FIX',
    'FITFUNCTION', 'FITPOINTS', 'MATCH', 'DISP', 'OPTICS', 'RING', 'EXPAND', 'RFSW', 'NORFSW', 'RAD', 'NORAD',
    'FLUC', 'NOFLUC', 'COD', 'NOCOD', 'INTRA', 'NOINTRA', 'RADCOD', 'RADTAPER', 'CODPLC', 'CONVERGENCE', 'TRPT',
    'TRACK', 'TRACKP', 'DAPERT', 'DAPWIDTH', 'SAVE', 'LOAD', 'READ', 'OUTPUT', 'APPEND', 'CLOSE', 'PRINT', 'ABORT',
    'QUIT', 'END', 'EXIT', 'STOP', 'GEO', 'GEOCAL', 'VISIT', 'ORG', 'SPLIT', 'MARK', 'SAVE', 'RESET', 'REJECT',
    'RENUM', 'STATUS', 'TYPE', 'SCALE', 'PROFILE', 'REAL', 'COUPLE', 'MOMENTUM', 'CHARGE', 'MASS', 'MINCOUP',
    'INS', 'WSPACSC', 'SPAC', 'MAXI', 'RADLIGHT', 'POL', 'POLARIZATION', 'SYNCHRON',
)
ELEMENTS = (
    'DRIFT', 'BEND', 'QUAD', 'SEXT', 'OCT', 'DECA', 'DODECA', 'MULT', 'SOL', 'CAVI', 'TCAVI', 'MARK', 'APERT',
    'COORD', 'LINE', 'UND', 'MAP', 'BEAMBEAM', 'MONI', 'MONITOR', 'MOVE', 'WAKE', 'LWAKE', 'TWAKE', 'PROF',
    'PROT', 'DISP', 'SPCH', 'STEER',
)
CONTROL = (
    'Module', 'Do', 'If', 'While', 'For', 'Switch', 'Which', 'With', 'Block', 'Return', 'Break', 'Continue',
    'Function', 'Scan', 'Table', 'Catch', 'Throw', 'Check', 'Clear', 'Null', 'True', 'False', 'Return',
)
BUILTINS = (
    'Print', 'Fit', 'Show', 'Plot', 'ListPlot', 'Emittance', 'TrackParticles', 'DynamicApertureSurvey', 'Twiss',
    'Element', 'LINE', 'Length', 'Range', 'Map', 'Apply', 'Select', 'Sqrt', 'Abs', 'Sin', 'Cos', 'Tan', 'Exp', 'Log',
    'Min', 'Max', 'Sum', 'Mod', 'Floor', 'Round', 'Flatten', 'Append', 'Join', 'Reverse', 'Sort', 'Part', 'First',
    'Last', 'Rest', 'ToString', 'ToExpression', 'StringJoin', 'StringLength', 'Open', 'Close', 'Read', 'Write',
    'WriteString', 'Get', 'FFS', 'Real', 'Pi', 'Sign', 'Re', 'Im', 'Conjugate', 'Random', 'Gamma', 'Erf', 'Erfc',
    'Bessel', 'Interpolate', 'Dot', 'Inverse', 'Transpose', 'Det', 'Eigensystem', 'SingularValues', 'Fold',
    'FitFunction', 'FFSKeywords', 'FrequencyMap', 'MomentumAcceptance', 'TouschekLifetime',
    'TouschekLifetimeLMA', 'ResonanceDrivingTerms', 'AmplitudeTuneShift', 'ChromaticityTerms', 'TransferMatrix',
    'BeamMatrix', 'OpticsPlot', 'Twiss', 'SetElement', 'LINE', 'GaussianCoulomb', 'PolynomialFit', 'Chop',
    'NaN', 'Infinity', 'Rad', 'Deg', 'Nest', 'Take', 'Drop', 'Position', 'Insert', 'Delete', 'Count', 'Cases',
    'MemberQ', 'Head', 'Dimensions', 'Array', 'Union', 'Intersection', 'Complement', 'Tally', 'Total', 'Mean',
    'Variance', 'StandardDeviation', 'Fix', 'Ceiling', 'ArcTan', 'ArcSin', 'ArcCos', 'Sinh', 'Cosh', 'Tanh',
)


class SADLexer(RegexLexer):
    name = 'SAD'
    aliases = ['sad', 'sadscript', 'ffs']
    filenames = ['*.sad']
    mimetypes = ['text/x-sad']
    flags = re.MULTILINE

    tokens = {
        'root': [
            (r'[ \t\r\f]+', Whitespace),
            (r'\n', Whitespace),
            (r'!(?!=)[^\n]*', Comment.Single),
            (r'\(\*', Comment.Multiline, 'comment'),
            (r'"(?:\\.|[^"\\\n])*"', String.Double),
            (r'(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?', Number),
            (words(ELEMENTS, prefix=r'(?<![\w$.%])', suffix=r'(?![\w$.%])'), Keyword.Type),
            (words(COMMANDS, prefix=r'(?<![\w$.%])', suffix=r'(?![\w$.%])'), Keyword),
            (words(CONTROL, prefix=r'(?<![\w$.%])', suffix=r'(?![\w$.%])'), Keyword),
            (words(sorted(set(BUILTINS)), prefix=r'(?<![\w$.%])', suffix=r'(?![\w$.%])'), Name.Builtin),
            (r'[A-Za-z_$%][\w$%.]*', Name),
            (r'->|:=|==|!=|<=|>=|//|&&|\|\||\+\+|--|[-+*/^<>=&|@#~?]', Operator),
            (r'[\[\]{}(),;:]', Punctuation),
            (r'.', Text),
        ],
        'comment': [
            (r'[^*)]+', Comment.Multiline),
            (r'\*\)', Comment.Multiline, '#pop'),
            (r'[*)]', Comment.Multiline),
        ],
    }
