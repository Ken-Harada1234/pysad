# var_PhysicalConstants (物理定数)

種別: 変数 (SADHelp の physical-constants)・フォルダ 10_math_data・ラティスなし

## 何を見る例か
FFS から参照できる物理定数 (`SpeedOfLight` `PlanckConstant` `PlanckHbar` `ElectronCharge` `FineStructureConstant` `ElectronMass` `ElectronRadius` `ProtonMass` `ProtonRadius` `SIMu0` `SIEpsilon0` `ElectronGminus2over2` `BoltzmannConstant`) を表示し、SADHelp が書く定義値・定義式と照合する。

## スクリプト
```
Print["R:SpeedOfLight=",SpeedOfLight];  ...  (13 個)
Print["R:re_formula=",FineStructureConstant*PlanckHbar*SpeedOfLight/(ElectronCharge*ElectronMass)];
```
`$FORM="S17.12"` (12 桁)。質量は eV 単位 (`ElectronMass` = 510998.95 eV)。

## 結果の読み方
`R:名前=値` の行。単位は SADHelp の表のとおり: `SpeedOfLight` m/s、`PlanckConstant` Js、`ElectronCharge` C、`ElectronMass`・`ProtonMass` eV、`ElectronRadius`・`ProtonRadius` m、`BoltzmannConstant` J/K。

## 物理・解析式
SADHelp の定義: $c=299792458$ m/s、$h=6.62607015\times10^{-34}$ Js、$\hbar=h/2\pi$、$e=1.602176634\times10^{-19}$ C、$\alpha=1/137.035999084$、$m_ec^2=0.51099895$ MeV、$m_pc^2=938.27208816$ MeV、
$$r_e=\frac{\alpha\hbar c}{e\,m_e}\ [\mathrm m],\qquad \mu_0=\frac{2\alpha h}{c\,e^2},\qquad \epsilon_0=\frac1{\mu_0c^2}$$
($m_e$ は eV 単位なので $e\,m_e$ が J になる)。`ProtonRadius` は $r_p=\alpha\hbar c/(e\,m_p)=r_e\,m_e/m_p$ と値が合う (定義の式は SADHelp に無いので、値から推定した)。

| 項目 | 式・SADHelp の定義値 (Python で計算) | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| SpeedOfLight $c$ [m/s] | 299792458 | 299792458 | 299792458 | 0.0e+00 | 0.0e+00 |
| PlanckConstant $h$ [Js] | 6.62607015e-34 | 6.62607015E-34 | 6.62607015E-34 | 0.0e+00 | 0.0e+00 |
| PlanckHbar $\hbar=h/2\pi$ | 1.054571818e-34 | 1.05457181765E-34 | 1.05457181765E-34 | 3.6e-12 | 3.6e-12 |
| ElectronCharge $e$ [C] | 1.602176634e-19 | 1.602176634E-19 | 1.602176634E-19 | 0.0e+00 | 0.0e+00 |
| FineStructureConstant $\alpha=1/137.035999084$ | 0.00729735256928 | .007297352569 | .007297352569 | 3.9e-11 | 3.9e-11 |
| ElectronMass $m_e$ [eV] | 510998.95 | 510998.95 | 510998.95 | 0.0e+00 | 0.0e+00 |
| ElectronRadius $r_e=\alpha\hbar c/(e m_e)$ [m] | 2.817940326e-15 | 2.81794032619E-15 | 2.81794032619E-15 | 1.4e-12 | 1.4e-12 |
| ProtonMass $m_p$ [eV] | 938272088.2 | 938272088.1599999 | 938272088.1599999 | 1.3e-16 | 1.3e-16 |
| ProtonRadius $=\alpha\hbar c/(e m_p)$ [m] | 1.534698267e-18 | 1.53469826718E-18 | 1.53469826718E-18 | 8.3e-13 | 8.3e-13 |
| SIMu0 $\mu_0=2\alpha h/(c e^2)$ | 1.256637062e-06 | 1.256637062121E-6 | 1.256637062121E-6 | 3.8e-14 | 3.8e-14 |
| SIEpsilon0 $=1/(\mu_0c^2)$ | 8.854187813e-12 | 8.85418781279E-12 | 8.85418781279E-12 | 3.4e-13 | 3.4e-13 |
| ElectronGminus2over2 $(g-2)/2$ | 0.00115965218128 | .001159652181 | .001159652181 | 2.4e-10 | 2.4e-10 |
| BoltzmannConstant [J/K] | 1.380649e-23 | 1.380649E-23 | 1.380649E-23 | 0.0e+00 | 0.0e+00 |
| $r_e$ を定義式で再計算 (SAD の変数から) | 2.817940326e-15 | 2.81794032619E-15 | 2.81794032619E-15 | 1.4e-12 | 1.4e-12 |
| $\mu_0$ を定義式で再計算 | 1.256637062e-06 | 1.256637062121E-6 | 1.256637062121E-6 | 3.8e-14 | 3.8e-14 |
| $\epsilon_0$ を定義式で再計算 | 8.854187813e-12 | 8.85418781279E-12 | 8.85418781279E-12 | 3.4e-13 | 3.4e-13 |

## gs と pysad の比較
全 16 行が 12 桁で一致 (○)。`ProtonMass` は gs・pysad とも `938272088.1599999` (浮動小数点の表現による。値は 938272088.16)。

## 未確認・限界
SADHelp の定義値は CODATA 2018 相当。`ProtonRadius` の定義式は推定。旧版の `ElectronMass` (510999.06) からの変更 (docs/32) は、この例の範囲では gs の値が 510998.95 であることまで確認。
