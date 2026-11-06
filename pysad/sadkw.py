"""
sadkw.py — SAD の要素型ごとのキーワード一覧 (SetElement / Element["KEYWORDS"] の順序)。
SAD (K64 gs) で `SetElement["Zname",type]` を型ごとに実行して得た一覧 (bench の gen_kw.sad)。値の既定はすべて 0 (下の NONZERO 以外)。
"""

# 型名 → キーワード (SAD の順)
KEYWORDS = {
    'MARK': ['AX', 'BX', 'PSIX', 'AY', 'BY', 'PSIY', 'EX', 'EPX', 'EY', 'EPY', 'R1', 'R2', 'R3', 'R4', 'DETR', 'DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP', 'AZ', 'BZ', 'PSIZ', 'ZX', 'ZPX', 'ZY', 'ZPY', 'DP', 'OFFSET', 'SIGMAZ', 'SIGE', 'GEO', 'JDX', 'JDPX', 'JDY', 'JDPY', 'JDZ', 'JDPZ', 'EMITX', 'EMITY', 'EMITZ', 'COUPLE'],
    'DRIFT': ['L', 'DISKIN', 'COUPLE', 'RADIUS', 'AX', 'AY', 'DISRAD'],
    'BEND': ['L', 'ANGLE', 'E1', 'E2', 'ROTATE', 'DISRAD', 'K1', 'DX', 'DY', 'K0', 'DISFRIN', 'EPS', 'RANKICK', 'F1', 'FRINGE', 'COUPLE', 'DROTATE', 'AE1', 'AE2', 'FB1', 'FB2', 'INDEX', 'AX', 'AY', 'LDEV', 'LRAD', 'CHI2', 'DZ'],
    'QUAD': ['L', 'K1', 'ROTATE', 'DX', 'DY', 'DISRAD', 'ACHROMA', 'DISFRIN', 'F1', 'F2', 'FRINGE', 'EPS', 'DISKIN', 'COUPLE', 'INDEX', 'AX', 'AY', 'LDEV', 'LRAD', 'F1K1F', 'F2K1F', 'F1K1B', 'F2K1B'],
    'SEXT': ['L', 'K2', 'ROTATE', 'DX', 'DY', 'DISRAD', 'DISFRIN', 'COUPLE', 'INDEX', 'AX', 'AY', 'LDEV', 'LRAD'],
    'OCT': ['L', 'K3', 'ROTATE', 'DX', 'DY', 'DISRAD', 'DISFRIN', 'COUPLE', 'INDEX', 'AX', 'AY', 'LDEV', 'LRAD'],
    'DECA': ['L', 'K4', 'ROTATE', 'DX', 'DY', 'DISRAD', 'DISFRIN', 'COUPLE', 'INDEX', 'AX', 'AY', 'LDEV', 'LRAD'],
    'DODECA': ['L', 'K5', 'ROTATE', 'DX', 'DY', 'DISRAD', 'DISFRIN', 'COUPLE', 'INDEX', 'AX', 'AY', 'LDEV', 'LRAD'],
    'MULT': ['L', 'DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'ROTATE', 'EPS', 'DISRAD', 'DISFRIN', 'F1', 'F2', 'FRINGE', 'VOLT', 'HARM', 'PHI', 'FREQ', 'COUPLE', 'RADIUS', 'DPHI', 'W1', 'ANGLE', 'E1', 'E2', 'DROTATE', 'DISK0FR', 'FB1', 'FB2', 'K0', 'SK0', 'K1', 'SK1', 'K2', 'SK2', 'K3', 'SK3', 'K4', 'SK4', 'K5', 'SK5', 'K6', 'SK6', 'K7', 'SK7', 'K8', 'SK8', 'K9', 'SK9', 'K10', 'SK10', 'K11', 'SK11', 'K12', 'SK12', 'K13', 'SK13', 'K14', 'SK14', 'K15', 'SK15', 'K16', 'SK16', 'K17', 'SK17', 'K18', 'SK18', 'K19', 'SK19', 'K20', 'SK20', 'K21', 'SK21', 'AE1', 'AE2', 'AX', 'AY', 'LDEV', 'LRAD', 'AUTOPHI', 'F1K1F', 'F2K1F', 'F1K1B', 'F2K1B', 'DVOLT', 'PROFILE'],
    'SOL': ['L', 'BZ', 'DX', 'DY', 'DZ', 'DPX', 'DPY', 'BOUND', 'CHI1', 'CHI2', 'CHI3', 'GEO', 'F1', 'DBZ', 'COUPLE', 'DISFRIN', 'DISRAD', 'FLAT'],
    'CAVI': ['L', 'VOLT', 'HARM', 'PHI', 'FREQ', 'RANVOLT', 'RANPHASE', 'LWAKE', 'TWAKE', 'DX', 'DY', 'ROTATE', 'V1', 'V20', 'V11', 'V02', 'COUPLE', 'DPHI', 'DISFRIN', 'FRINGE', 'AX', 'AY', 'LDEV', 'AUTOPHI', 'DVOLT'],
    'TCAVI': ['L', 'K0', 'HARM', 'PHI', 'FREQ', 'DX', 'DY', 'ROTATE', 'RANKICK', 'RANPHASE', 'LWAKE', 'TWAKE', 'COUPLE', 'AX', 'AY', 'LDEV'],
    'MONI': ['DX', 'DY', 'OFFSET', 'COUPLE', 'ROTATE'],
    'APERT': ['DX1', 'DX2', 'DY1', 'DY2', 'JDPX', 'JDPY', 'DP', 'COUPLE', 'ROTATE', 'AX', 'AY', 'DX', 'DY'],
    'COORD': ['DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3', 'DIR', 'COUPLE'],
}

# 0 でない既定値
NONZERO = {}

# 型名 → 既定キーワード (Element["DEFAULT",name])。MARK / MONI は空
DEFAULT_KEY = {'MARK': '', 'DRIFT': 'L', 'BEND': 'ANGLE', 'QUAD': 'K1', 'SEXT': 'K2', 'OCT': 'K3', 'DECA': 'K4', 'DODECA': 'K5',
               'MULT': 'K1', 'CAVI': 'VOLT', 'TCAVI': 'K0', 'MONI': ''}

# 型名 → 内部コード (Element["TYPE",name])
TYPECODE = {'DRIFT': 1.0, 'BEND': 2.0, 'QUAD': 4.0, 'SEXT': 6.0, 'OCT': 8.0, 'DECA': 10.0, 'DODECA': 12.0, 'SOL': 20.0, 'MULT': 22.0,
            'CAVI': 31.0, 'TCAVI': 32.0, 'MARK': 41.0, 'MONI': 42.0}
