# SAD/FFS & SADScript のコマンド・関数・変数一覧 (SADHelp 1.4.4.2k64 Updated: 10/19/2026 から機械的に抽出)

出典: ローカル `<SAD>\Documents\SADHelp.html` (1.4.4.2k64 Updated: 10/19/2026。CERN の https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html は 1.4.0.0k64 2025-07-14 の同系統)。各行は見出しと、その説明の最初の行 (Usage があればその形)。括弧内は最短の省略形。pysad での対応状況は docs/03_ffs_commands.md §7 と docs/04_sadscript.md を参照。

## 1. FFS コマンド (56)

| コマンド | 説明 (最初の行) |
|---|---|
| `ABORT` | Terminates SAD immediately. |
| `APPEND(APP)` | APP {filename \| file-number} switches the output stream to the specified file or the file number. The output is appended to th |
| `ATTRIBUTE(ATTR)` | Usage: ATTR element-pattern |
| `BEAMSIZE(BEAM)` | Calculate the beam size with the current Twiss parameters. The calculation is in 5D and not correct if the synchrotron motion |
| `BYE` | Exits from the current beam line and returns to the original beam line where VISIT command was issued. All information specifi |
| `CALCULATE(CAL)` | Usage: (1) CALCULATE(CAL) [[NO]EXPAND]] |
| `CHROMATICITY(CHRO)` | CHRO prints out the chromaticity of QUAD and SEXT in the entire beam line using the simplest formula: |
| `CLOSE(CLO)` | CLOSE [INPUT(IN)] closes the current input stream and switches it to the previous input stream. |
| `COUPLE(COUP)` | Usage: COUP slave-element master-element coefficient |
| `DISPLAY(DISP)` | Usage: DISP_LAY [keywords] [pattern-string] [region] |
| `DRAW` | Usage: DRAW [begin end] fun1 [fun2..] [& fun11 [fun12..]] [element-pattern] |
| `DUMP` | Usage: DUMP component-pattern [component-pattern1..] |
| `ELSE` | Usage: IF expr1 body1 [ELSEIF expr2 body2 [ELSEIF ..]] [ELSE body3] ENDIF |
| `ELSEIF` | Usage: IF expr1 body1 [ELSEIF expr2 body2 [ELSEIF ..]] [ELSE body3] ENDIF |
| `EMITTANCE(EMIT)` | Usage: (1) EMIT |
| `END` | Closes the current output-stream and set the output stream to the standard output(6). It also suspends all the input streams a |
| `ENDIF` | Usage: IF expr1 body1 [ELSEIF expr2 body2 [ELSEIF ..]] [ELSE body3] ENDIF |
| `EXECUTE(EXEC)` | Usage: EXEC character-string-expression |
| `EXPAND` | EXPAND distributes all variable keys to each corresponding components. Individual deviations (machine errors) are cleared. |
| `FIT` | Usage: (1) FIT [component] |
| `FITPOINTS(FITP)` | Usage: FITP n |
| `FIX` | Usage: (1) FIX element-pattern [keyword] [element-pattern1 [keyword1]..] |
| `FREE` | Usage: FREE element-pattern [keyword] [element-pattern1 [keyword1]..] |
| `GO` | Usage: GO [[NO]EXPAND] |
| `IF` | Usage: IF expr1 body1 [ELSEIF expr2 body2 [ELSEIF ..]] [ELSE body3] ENDIF |
| `INPUT(IN)` | IN {filename \| file-number} switches the input stream to the specified file or the file-number. The original stream is kept an |
| `MATRIX(MAT)` | Usage: MATRIX [{SYMPLECTIC(S) \| PHYSICAL(P)}] [from to] |
| `MEASURE(MEA)` | Usage: MEA [end-component] [OUT file plot-spaces] |
| `ORG` | Usage: ORG location, dgx, dgy, dgz, dchi1, dchi2, dch3 |
| `OUTPUT(OUT)` | OUT {filename \| file-number} switches the output stream to the specified file or the file-number. The file is written from the |
| `PRINT(PRI)` | PRI expression evaluates expression and prints out the result. |
| `QUIT` | Exits FFS and return to SAD/MAIN level, without saving the values of the elements. |
| `RADINT` | RADINT prints out the radiation integrals involving the x-coupling for all components of the beam line. |
| `READ` | READ {filename \| file-number} switches the input stream to the specified file or the file-number. The original stream is kept |
| `RECOVER(REC)` | REC exchanges the values of FREEd elements with those when the last GO command was issued. FIXed elements are not affected. |
| `REFERENCE(REF)` | REFERENCE(REF) sets the current optics as the reference optics. |
| `REJECT(REJ)` | Usage: (1) REJ matching-function-pattern [matching-function-pattern1..] |
| `RENUMBER(RENUM)` | RENUM comp renumbers the component number starting from a component comp. |
| `REPEAT(REP)` | Usage: REP [n] body UNTIL [expr1] |
| `RESET` | Usage: RESET [ALL] [element-pattern] |
| `RESUME(RES)` | Resumes reading from the previous input stream suspended by SUSPEND(SUSP) or END. |
| `REVERSE(REV)` | REV changes the sign of AX, AY, EPX, EPY, R2, R3, DPX, DPY at the |
| `SAVE` | Usage: SAVE [element-pattern] |
| `SEED` | The SEED command is obsolete. Use SeedRandom[] function instead of SEED. |
| `SHOW` | SHOW prints out the current matching conditions. |
| `SPLIT` | Usage: SPLIT component length |
| `STATUS(STAT)` | STAT shows the current settings of flags, fit points, |
| `STOP` | Exits FFS and returns to SAD/MAIN level, with saving the values of the elements. |
| `SUSPEND(SUSP)` | Suspend reading from the current input stream, and wait for input from the console. Resumes by RESUME(RES) |
| `TERMINATE(TERM)` | TERM [INPUT(IN)] suspends the current input stream and switches it to the previous input stream. |
| `TYPE(T)` | Usage: TYPE [element-pattern [element-pattern1..]] |
| `UNTIL` | Usage: REP [n] body UNTIL [expr1] |
| `USE` | Usage: USE [[NO]EXPAND] beam-line |
| `VARIABLES(VAR)` | VARIABLES displays a list of current matching-variables and their present, previous, saved, minimum, and maximum values togeth |
| `VARY` | Usage: VARY keyword element-pattern |
| `VISIT` | Usage: VISIT [[NO]EXPAND] beam-line |

## functions (108)

FFS functions:

### Data-Manipulation (6)

| 名前 | 説明 (最初の行) |
|---|---|
| `Fit` | Usage: Fit[data, expr, var, |
| `FitEmit` | Obtain the emittance, beta, alpha, etc. from data of particles x and px using FitGaussian. |
| `FitGaussian` | Performs Gaussian fit of 1D list data: |
| `NIntegrate` | NIntegrate returns numerical integration of a Real or Complex |
| `PolynomialFit` | Usage: PolynomialFit[data, n] |
| `Spline` | Spline returns data for cubic-spline interpolation. |

### DownhillSimplex

| 名前 | 説明 (最初の行) |
|---|---|
| `DownhillSimplex` | Usage: DownhillSimplex[initial, f, options] |

### functional-operations (15)

| 名前 | 説明 (最初の行) |
|---|---|
| `Apply` | Apply[f, x [, level]] or f@@x [f@@[x, level]] |
| `Cases` | Cases[l, pat [, level [,n]]] |
| `DeleteCases` | DeleteCases[l, pat [, level [,n]]] |
| `Difference` | Difference[list] returns Rest[list] - Drop[list, -1] . |
| `FixedPoint` | FixedPoint[f, e] starts from f[e] then applies f repeatedly until the result no longer changes. FixedPoint[f, e, n] specifies |
| `FixedPointList` | FixedPointList[f, e] returns the intermediate values of FixedPoint[f, e] as a list. FixedPointList[f, e, n] and options for Fi |
| `level-spec` | Functions Map, Apply, Scan, Cases, Position, Count, Level, DeleteCases takes an optional argument level to specify which level |
| `Map` | Map[f, x [, level]] or f/@x [f/@[x, level]] |
| `MapThread` | MapThread[f, l] ===> f@@[Thread[l], {1}] |
| `Position` | Position[l, pat [, level [,n]]] |
| `Scan` | Scan[f, x [, level]] |
| `ScanThread` | ScanThread[f, l] ===> (f@@[Thread[l], {1}];Null) , |
| `SelectCases` | Usage: SelectCases[list, {test1,..}] |
| `SwitchCases` | Usage: SwitchCases[list, {case1,..}] |
| `Thread` | Thread[l] returns a threaded list of l: |

### FFS-dedicated-functions (31)

Functions dedicated to the optics calculations and simulations in FFS.

| 名前 | 説明 (最初の行) |
|---|---|
| `AccelerateParticles` | AccelerateParticles does TrackParticles with acceleration in a ring for a given number of turns. The adiabatic damping is auto |
| `BeamMatrix` | BeamMatrix[i] returns a 6 by 6 beam-matrix (i.e., <x_k x_l>) at location i. The index i can have a fraction to specify interme |
| `DynamicApertureSurvey` | Usage: DynamicApertureSurvey[range,nturn,options] |
| `Element` | Element[key-string, {element-pattern-string \| element-position}] returns values for key-string of elements which match element |
|   `key-strings:Element` | The key-string is not case-sensitive. Available key-strings are: |
| `Emittance` | Emittance[option] returns a set of rules as |
| `ExternalMap` | With MAP elements, ExternalMap defines a user-defined map of particles. It also allows a user to do anything (doing statistics |
|   `CompiledMap` | ExternapMa["TRACK", n, nt_, x_]:=CompiledMap[nt, x, src, prm, opt___]; |
| `FFS` | FFS[command-string] executes command-string as FFS commands. Any commands can be used. Some commands CALCULATE(CAL), GO, VARIA |
| `FFS$SHOW` | FFS["SHOW"] or FFS$SHOW[] returns the current matching conditions as a list. Each element has a form of |
| `FitValue` | (1) FitValue[component, function, {id_,dp_}, goal_, now_] := body |
| `FitWeight` | A defined function to modify the weight of matching of particular function at particular component with particular momentum o |
| `GaussianCoulomb` | EquAlign["Text[{\\tt GaussianCoulomb[x, y, sigx, sigy]} returns the electromagnetic force $f_x+if_y$ at the coordinates $(x,y\ |
| `GeoBase` | GeoBase[{chi1, chi2, chi3}] converts the rotation angles of the coordinate base vectors to a transformation matrix {{x_gx,x_gy |
| `LINE` | LINE[key-string, {component-pattern-string \| component-position}] |
|   `key-strings:LINE` | The key-string is not case-sensitive. Available key-strings are: |
| `OptimizeOptics` | Usage: OptimizeOptics[options] |
| `OrbitGeo` | OrbitGeo[location] returns the geometry {GX, GY, GZ} of the current (not design) orbit. |
| `RadiationField` | To calculate the field of the synchrotron radiation from particles, first record trajectories of particles. This is done by th |
| `RadiationSpectrum` | To calculate the field of the synchrotron radiation from particles, first record trajectories of particles. This is done by th |
| `SetElement` | Create/set/read a MAIN-level element. |
| `SurvivedParticles` | SurvivedParticles[x] |
| `SymplecticJ` | SymplecticJ[n] returns an n by n symplectic matrix: |
| `SynchroBetaEmittance` | SynchroBetaEmittance calculates equilibrium emittance under influence of synchrotron motion and chromaticity. |
| `TouschekLifetime` | TouschekLifetime interpolates the data calculated by EMIT or Emittance[]. There are three ways of usage: |
| `TrackParticles` | TrackParticles[beam, destination-component, nbegin, nend] |
| `TransferMatrix` | TransferMatrix[from_, to_, options___] returns the 4x4 or 6x6 transfer matrix (TM) from component 'from' to component 'to', de |
| `Twiss` | Twiss[optical-function, component] returns the value of the optical-function at the entrance of component. The values are tho |
| `VariableRange` | Usage: VariableRange[element, keyword ,v_] := expression |
| `VariableWeight` | Usage: VariableWeight[element, keyword ,v_] := expression |
| `WakeFunction` | WakeFunction[Longitudinal, comp]={{z1, wl1}, ..., {zn, wln}}; |

### Graphics (13)

Graphics represents an object for 2D graphics with the form

| 名前 | 説明 (最初の行) |
|---|---|
| `BeamPlot` | Usage: BeamPlot[loc, axes, options] |
| `ColumnPlot` | Usage: ColumnPlot[data, options, ...] |
| `FitPlot` | Usage: FitPlot[data, fun, var, {par1,ini1}, .. , {parn, inin}, opt]. |
| `GeometryPlot` | Usage: GeometryPlot[options] |
| `HistoPlot` | Usage: HistoPlot[data, options, ...] |
| `LaTeX` | Some Graphics options and the Text primitive accept a LaTeX input, if the string |
| `ListContourPlot` | Usage: ListContourPlot[list, options, ...] |
| `ListDensityPlot` | Usage: ListDensityPlot[list, options, ...] |
| `ListPlot` | Usage: ListPlot[{{x1,y1},..,{xn,yn}}, options] |
| `OpticsPlot` | Usage: OpticsPlot[fun_list, options] |
| `ParametricPlot` | Usage: ParametricPlot[{funx,funy}, range, options] or Plot[{{funx1, funy1} .. }, range, options], |
| `Plot` | Usage: Plot[fun, range, options] or Plot[{fun1, .. }, range, options], |
| `PointSymbol` | Notations for symols to be plotted: |

### Input/Output (15)

| 名前 | 説明 (最初の行) |
|---|---|
| `$FORM` | $FORM is a character-string to specify the format of the output of a real number. |
| `$Input` | $Input holds the file number for console input stream. The default is -1. |
| `$Output` | $Output holds the file number for console output. The default is -1. |
| `Close` | Close[f] closes file number [f]. It is necessary to complete the output to an external file. |
| `OpenAppend` | f = OpenAppend[file] |
| `OpenRead` | f = OpenRead[file] |
| `OpenWrite` | f = OpenWrite[file] |
| `PageWidth` | PageWidth is the number of columns of the output. The default is set from GetEnv["WIDTH"]. |
| `Print` | Print[expr1 [,expr2 ...]] |
| `Read` | Read[f, item [, item1...] [, opts..]] |
| `ReadString` | ReadString[f] |
| `StandardForm` | StandardForm[expr] |
| `StringToStream` | f = StringToStream[string] |
| `Write` | Write[f, expr1 [,expr2 ...]] |
| `WriteString` | WriteString[f, expr1 [,expr2 ...]] |

### Multiprocessing (4)

| 名前 | 説明 (最初の行) |
|---|---|
| `Fork` | Forks the process into a parent and a child processes. |
| `OpenShared` | Allocated shared memory of n bytes. |
| `Shared` | Read/Write to the shared memory. |
| `SharedSize` | Returns the size of an object for OpenShared. |

### Object-oriented-programing (1)

Environment for an object-oriented-programming is supplied by:

| 名前 | 説明 (最初の行) |
|---|---|
| `Class` | Class sets up a class of objects. |

### Random-number-functions (5)

The random number functions use common seed given by SEED command or SeedRandom function. It has an initial value 17 at the be

| 名前 | 説明 (最初の行) |
|---|---|
| `Random` | Random[] gives a uniform random number between 0 and 1. |
| `GaussRandom` | GaussRandom[] gives a Gaussian random number with average 0, |
| `ParabolaRandom` | ParabolaRandom[] gives a parabola random number between -1 and 1. |
| `SeedRandom` | SeedRandom[plugin_String] selects new pseudo random-number generator |
| `ListRandom` | ListRandom[] returns List of available pseudo random-number generator plugins. |

### System-interface (2)

| 名前 | 説明 (最初の行) |
|---|---|
| `System` | System[string] executes string as a shell command. It returns the system return code. |
| `TemporaryName` | TemporaryName[] returns a unique file name for a temporary use. |

### Utilities (5)

| 名前 | 説明 (最初の行) |
|---|---|
| `DateString` | DateString[] returns the current date and time as string "mm/dd/CCYY HH:MM:SS". |
| `MemoryCheck` | MemoryCheck[n] |
| `TimeUsed` | TimeUsed[] returns the cputime since the start of SAD in seconds. |
| `Timing` | Timing[fun] |
| `TracePrint` | TracePrint[fun] |


## special-variables (48)

There are several variables which have special rolls in FFS. Some of them are also accessible in the MAIN level.

| 名前 | 説明 (最初の行) |
|---|---|
| `$Line` | %Line counts the number of the results in FFS Shown with Out[]:= . Out remembers all outputs up to Out[$Line]. $Line = 0 rese |
| `CASE` | CASE is a character-string to be attached with TITLE to issue the CASE command of TopDrawer in DRAW or GEO commands. |
| `CHARGE` | CHARGE contains the charge of the particle. The default is +1. |
| `CONVERGENCE` | CONVERGENCE is the goal of the convergence(==MatchingResidual) in the matching. If MatchingResidual becomes smaller than CONVE |
| `DAPWIDTH` | DAPWIDTH is the successive width of the x-amplitudes to terminate the tracking assuming that the aperture is enough. The defau |
| `DP` | DP represents the relative momentum spread of the beam. It is automatically set by the keyword DP of the MARK element at the b |
| `DP0` | DP0 represents the central value of the relative momentum offset in the optics calculation, or the center of the momentum dis |
| `DTSYNCH` | DTSYNCH is the equilibrium arrival time advance (= -c*delta t, in meter) where the synchrotron radiation loss balances with t |
| `EFFRFFREQ` | The effective angular frequency weff (= 2 Pi EFFRFFREQ) is obtained using the second derivative |
| `EFFVC` | Effective peak rf voltage EFFVC is given by |
| `EFFVCRATIO` | Ratio of (effective voltage)/(Sum[VOLT_k,{k}]) (default: 1). Set by EMITTANCE(EMIT) or Emittance[]. Effective with RING only. |
| `ElementValues` | ElementValues is a symbol to assign rules to determine values of keywords of elements or components. This is used to give a de |
| `EMITX` | EMITX is a real variable for the horizontal emittance (not normalized by gamma beta). It is automatically set by the keyword E |
| `EMITXE` | EMITXE specifies the minimum horizontal emittance for EMITTANCE(EMIT) and Emittance[] calculations. It is useful to give emitt |
| `EMITY` | EMITY is a real variable for the vertical emittance (not normalized by gamma beta). It is automatically set by the keyword EMI |
| `EMITYE` | EMITYE specifies the minimum vertical emittance for EMITTANCE(EMIT) and Emittance[] calculations. It is useful to give emittan |
| `EMITZ` | EMITZ is a real variable for the longitudinal emittance (not normalized by gamma beta). It is automatically set by the keyword |
| `EMITZE` | EMITZE specifies the minimum longitudinal emittance for EMITTANCE(EMIT) and Emittance[] calculations. It is useful to give emi |
| `ExponentOfResidual` | The exponent to calculate MatchinResidual. The default is 2. |
| `FFS$NumericalDerivative` | If False (default), the calculation of response matrix for each matching variables uses analytical expressions as much as poss |
| `FitFunction` | FitFunction is a symbol to assign user-defined functions for matching with the GO command. |
| `FSHIFT` | FSHIFT is the relative shift df/f0 of the revolution (or rf) frequency in a ring. This is valid in EMITTANCE(EMIT), the partic |
| `GCUT` | GCUT specifies the cut-off value of Gaussian distribution in unit of the standard deviation. Accessible in MAIN. |
| `InitialOrbits` | Usage: InitialOrbits = { {x1, px1, y1, py1, z1, dp1}, ...}; |
| `LOSSAMPL` | LOSSAMPL is the transverse amplitude beyond which a particle is judged to have been lost. The default is 1 m. Accessible in MAIN. |
| `LOSSDZ` | LOSSDZ is the longitudinal position z beyond which a particle is judged to have been lost. The default is 100 m. Accessible in |
| `MatchingAmplitude` | MatchingAmplitude is a list of amplitudes for the finite-amplitude matching. |
| `MatchingResidual` | MatchingResidual holds the convergence in the last GO or CALCULATE(CAL) commands. It is calculated by |
| `MASS` | MASS is the rest mass of the particle in eV. The default is the electron mass. |
| `MINCOUP` | MINCOUP is the minimum emittance ratio to be assumed in the calculation of intrabeam effects in EMITTANCE(EMIT) and Emittance[ |
| `MOMENTUM` | MOMENTUM is the nominal momentum of the beam line at the entrance in eV/c. Accessible in MAIN. |
| `NBUNCH` | NBUNCH is the number of bunches for the calculation of WakeFunction. Accessible in MAIN. |
| `NetResidual` | The net residual of convergence except the penalty for unstable optics. |
| `NP` | NP is the number of particles in the tracking. Accessible in MAIN. |
| `NPARA` | NPARA specifies the maximum number of parallel processes for various calculations such as CALC, GO, TrackParticles, DynamicApe |
| `OffMomentumWeight` | Relative weight of the off-momentum deviation for the off-momentum matching. The default is 1. |
| `OMEGA0` | OMEGA0 is 2*Pi*SpeedOfLight/LINE["s","$$$"] . Accessible in MAIN. |
| `OpticsEpilog` | OpticsEpilog is a variable to assign a user-defined function which is to be executed every time after an optics calculation is |
| `OpticsProlog` | OpticsProlog is a variable to assign a user-defined function which is to be executed every time before an optics calculation i |
| `PBUNCH` | PBUNCH is the number of particles/bunch for the calculation of the intra-beam and space charge effects in EMITTANCE(EMIT), and |
| `PhotonList` | When PHOTONS is ON (default is OFF), with RAD and FLUC, TrackParticles generates a list of all photons radiated through the tr |
| `PHICAV` | The equilibrium RF phase (in radian, = 2 Pi EFFRFFREQ DTSYNCH/c) where the synchrotron radiation loss and acceleration balance |
| `SIGE` | SIGE is the equilibrium momentum spread at the entrance of the beam line calculated by EMIT/Emittance[SaveEMIT->True]. |
| `SIGZ` | SIGZ is the equilibrium bunch length at the entrance of the beam line calculated by EMIT/Emittance[SaveEMIT->True]. |
| `SpeedOfLight` | SpeedOfLight is 299792458. |
| `TITLE` | TITLE is a character-string to make the title of the plot in DRAW or GEO commands. |
| `Track$Excitation` | An artificial random beam exciation is added in the tracking, by setting |
| `UnstableLevel` | Number of unstable planes in the optics calculations. Equal to zero if all x and y optics are stable for on/off-momentum and f |

## constants (0)

There are pre-defined special symbols for constants in FFS:

```
There are pre-defined special symbols for constants in FFS:

symbol         value
True           1
False          0
Infinity       INF
INF            INF
NaN            NaN
Pi             ArcSin[1]*2
E              Exp[1]
I              Complex[0,1]
Degree         Pi/180
GoldenRatio    (1+Sqrt[5])/2
EulerGamma     0.57721566490153286061
```

## physical-constants (0)

EquAlign["Text[Physical constants available in FFS are:\n\\begin{align*}\n{\\tt SpeedOfLight} &\\qquad c\\equiv 299792458 {\\\

```
EquAlign["Text[Physical constants available in FFS are:\n\\begin{align*}\n{\\tt SpeedOfLight} &\\qquad c\\equiv 299792458 {\\\
rm\\ m/s}\\\\\n{\\tt PlanckConstant} &\\qquad h \\equiv 6.62607015\\times 10^{-34} {\\rm\\ Js}\\\\\n{\\tt PlanckHbar} &\\qquad \\h\
bar\\equiv h/(2\\pi)\\\\\n{\\tt ElectronCharge} &\\qquad  e_e \\equiv 1.602176634\\times 10^{-19} {\\rm\\ C}\\\\\n{\\tt FineStruct\
ureConstant}&\\qquad  \\alpha = 1/137.035999084 \\\\\n{\\tt ElectronMass}&\\qquad  m_e = 0.51099895000\\times 10^6 {\\rm\\ eV}\\\\\
\n{\\tt ElectronRadius}&\\qquad  {\\rm classical\\ radius\\ of\\ electron\\ in\\ m,} r_e \\equiv \\alpha \\hbar c / (e_e m_e)\\\\\\
n{\\tt ProtonMass}&\\qquad  m_p  = 938.27208816\\times 10^6 {\\rm\\ eV}\\\\\n{\\tt ProtonRadius}&\\qquad  {\\rm classical\\ radius\
\\ of\\ proton\\ in\\ m}\\\\\n{\\tt SIMu0}&\\qquad \\mu_0\\equiv 2\\alpha h/(c e_e^2)\\\\\n{\\tt SIEpsilon0}&\\qquad    \\varepsil\
on_0\\equiv 1/(\\mu_0 c^2)\\\\\n{\\tt ElectronGminus2over2}&\\qquad   (g-2)/2 {\\rm\\ of\\ electron} = 0.001159652181280002\\\\\n{\
\\tt BoltzmannConstant}&\\qquad   1.380649\\times 10^{-23} {\\rm\\ J/K .}\n\\end{align*};;;0]",
"Physical constants available in FFS are:\n\nSpeedOfLight           c = 299792458 m/s\nPlanckConstant         h = 6.62607015d-34 J\
s\nPlanckHbar             hbar = h/2/Pi\nElectronCharge         e_e = 1.602176634d-19 C\nFineStructureConstant  alpha = 1.d0 / 137\
.035999084d0\nElectronMass           m_e in eV = 0.51099895000d6\nElectronRadius         classical radius of electron = alpha hbar\
 c / (e_e m_e), in m\nProtonMass             m_p in eV = 938.27208816d6\nProtonRadius           classical radius of Proton, in m\n\
SIMu0                  2 alpha h / e_e^2 / c\nSIEpsilon0             1/SIMu0/c^2\nElectronGminus2over2   (g-2)/2 of electrons = 0.\
001159652181280002\nBoltzmannConstant      1.380649d-23 J/K"
]
```

## その他の章 (19)

| 章 | 下位項目数 | 説明 (最初の行) |
|---|---|---|
| `beam-line` | 7 | A beam line is defined in the MAIN level by LINE command as: |
| `character-string` | 13 | A character-string is expressed by enclosing in "". Special characters are expressed using \: |
| `command-syntax` | 0 | The command syntax in FFS is |
| `components` | 0 | Components are the objects which consist the beam line. A component simulates an individual magnet, drift space, or rf-cavity. |
| `data-structure` | 5 | All data and "programs" in SAD Script are expressed either by an atom or a list-structure: |
| `defining-functions` | 0 | A function is defined by one of the following forms: |
| `dynamics` | 13 |  |
| `elements` | 170 | An element in FFS represents an object which has a unique name and several keywords with values. This simulates a power supply |
| `expression` | 46 | An expression in FFS consists of a symbol, constants, and operators. |
| `flags` | 60 | Usage: [NO]flag |
| `geometric-functions` | 0 | Available geometric-functions are: |
| `machine-error-commands` | 0 | Usage: machine-error-command [options] amount component-pattern .. |
| `matching-function-commands` | 0 | Usage: (1) matching-function goal-value [off-momentum-points] |
| `multi-turn-tracking` | 0 | The multi-turn tracking can be done by TrackParticles or DynamicApertureSurvey[] function in FFS. |
| `off-momentum-matching` | 0 | FFS matches the optical functions for an orbit with finite momentum deviation. |
| `optical-functions` | 0 | Available optical functions for matching are: |
| `pattern` | 1 | Pattern is a special expression for mathing arguments in function definitions and rules with several forms: |
| `set-value-of-element` | 2 | Usage: element-pattern [keyword] [{MIN \| MAX \| MINMAX \| MAXMIN}] value |
| `wildcards` | 0 | Many commands and functions accept the wildcards as a specification for the name of elements or components. The valid wildc |

### beam-line の下位項目 (7)

| 名前 | 説明 (最初の行) |
|---|---|
| `beam-line-functions` | Functions/objects to construct/edit beam lines and elements in FFS. |
|   `BeamLine` | Usage: BeamLine[e1, e2, ...]; |
|   `BeamLineName` | BeamLineName[] returns the name of the current beam line. If a BeamLine object is used by USE or VISIT, the new beam line beco |
|   `ExtractBeamLine` | Usage: ExtractBeamLine[line] |
|   `PrintBeamLine` | Usage: PrintBeamLine[b1,.. ,option] |
|   `WriteBeamLine` | Usage: WriteBeamLine[f, b1,.. ,option] |
| `orientation-of-an-element` | An element with negative orientation means a reversal of the element along the z-axis. Thus all magnets except for a solenoid |

### character-string の下位項目 (13)

| 名前 | 説明 (最初の行) |
|---|---|
| `FromCharacterCode` | FromCharacterCode[r_Real] returns a character whose character code is r. |
| `StringFill` | StringFill[s, sf, n] with strings s and sf, n > 0, returns (s//sf//sf...)[1,n] . |
| `StringJoin` | StringJoin[s1, s2, [,s3...]] ===> s1 // s2 [//s3...] concatenates strings s1, s2 [,s3...]. |
| `StringLength` | StringLength[s] returns the length of string s. |
| `StringMatchQ` | StringMatchQ[s, spat] returns True/False whether string s matches string-patten spat. |
| `StringPart` | s_String[n] returns the n-th character in s. |
| `StringPosition` | StringPosition[s, subs] returns a list of positions of subs in string s. |
| `StringReplace` | StringReplace[s, rules] replaces the parts of string s accoding to rules, which is a Rule or a list of Rules: |
| `StringTrim` | StringTrim[s] removes the leading and trailing spaces and tabs from s. |
| `Symbol` | Symbol[s] returns a Symbol whose name is character-string s. |
| `ToCharacterCode` | ToCharacterCode[s] returns the list of character codes of character-string s. |
| `ToExpression` | ToExpression[s] converts a character-string s to an expression and evaluate it. |
| `ToString` | ToString[expr] evaluates an expression expr, then converts to a character-string. |

### command-syntax (本文)

```
The command syntax in FFS is

     expression1 [param1..] [;] expression2..	

(1) The input is first evaluated as an expression. If the expression returns a Symbol with the same name as the expression itself, 
it is interpreted as an FFS command, otherwise the returned value is printed out unless it is Null.
(2) Each command takes succeeding its parameters if necessary. A command with indefinite number of parameters can be terminated by 
semicolon. Most commands terminate itself at the end of line.
(3) A line can be continued to the next line if a backslash is placed at the end of the line.
(4) An expression continues to the next line if it is not closed in the line.
(5) An exclamation mark comments out the rest of the line.	

Example: A command line

          QF* .1

means the set-value-of-element command as unless the symbol QF has been defined otherwise. If QF has been defines as a number, such
 as QF=2.5, the above command line is interpreted as Times[QF,.1] then returns .25 .
```

### components (本文)

```
Components are the objects which consist the beam line. A component simulates an individual magnet, drift space, or rf-cavity.
 The parameters of a component is specified the values in the corresponding element with the same name as the component, which simu
lates a power supply. Many components can be attached to the same element. Parameters of each component may deviate from the corres
ponding element if machine errors are given.
   A component is specified with the form name[.order][{+-}offset], where name is the name of the component. The number order means
 the order-th component which belongs to name element, counted from the beginning of the line starting from 1. Offset is a positive
 or negative number to specify the downstream or upstream components from the given component. If order is omitted, the first eleme
nt is assumed, and if offset is omitted, zero is assumed. The order can be renumbererd by RENUMBER(RENUM).
   The end of line is specified by $$$. The first component can be specified by ^^^.
```

### data-structure の下位項目 (5)

| 名前 | 説明 (最初の行) |
|---|---|
| `Extract` | Extract[f, part [,head]] |
| `Head` | Head[f] takes the head of an expression f. |
| `Length` | Length[f] returns the number of elements in the body of a structure f. |
| `List` | List is a special symbol to be the head of generic list-structure. |
| `Part` | Part[f, a [,b ,...]] ===> f[[a, [,b ...]]] |

### defining-functions (本文)

```
 A function is defined by one of the following forms:

f[pat1 [,pat2...]] (:)= body;
f[pat1 [,pat2...]] ^(:)= body;
g/:f[pat1 [,pat2...]] (:)= body;

where pat1 [,pat2...] are patterns (including expressions).
 If UpSet(^=) or UpSetDelayed (^:=) is used, the definition is associated with the symbol in the first level of l.h.s.
 If TagSet(/:) is used, the definition is associated with the symbol on the left of /: .

 The patters can be an expression including constants. The definition with constant arguments can be accessed faster than searching
 a list, in general, so they are suitable for a data base. Definitions with constant arguments have higher priorities than with pat
terns.
```

### dynamics の下位項目 (13)

| 名前 | 説明 (最初の行) |
|---|---|
| `independent-variable` | EquAlign["Text[SAD uses $s$, the distance along a reference line as the independent variable. The reference line is either a \ |
| `Lagrangean` | EquAlign["Text[The primary position variables are $(x,y,s)$, where $x$ and $y$ are the displacements along the normal and bin\ |
| `Hamiltonian` | EquAlign["Text[The Lagrangean $L$ defines the {\\it canonical} momenta as\n\\begin{align}\np_x=&\\frac{\\partial L}{\\partial\ |
|   `2nd-order-Hamiltonian` | EquAlign["Text[The Hamiltonian $H$ can be expanded up to the second order of the transverse coordinates as:\n\\begin{align}\n\ |
|     `solution-H2` | EquAlign["Text[We assume the focusing is uniform and $K=B'/(B\\rho)>0$, which can be always obtained by an appropriate rotati\ |
|     `solution-dH` | EquAlign["Text[The transformation for the correction term $\\varDelta H=H-H_2$ for $1/\\rho=0$:\\begin{align}\\varDelta H&=-\\ |
| `remarks-on-dynamics` | EquAlign["Text[There are several remarks on the dynamics in SAD:\n\\begin{itemize}\n\\item All transformations in SAD are sym\ |
| `x-y-coupling` | EquAlign["Text[The transformation matrix from the physical coordinate $(x,p_x,y,p_y)$ to the $x$-$y$ decoupled coordinate $(X\ |
| `extended-Twiss-parameters` | A symplectic matrix such as the normal mode matrix can be expressed in terms of the extended Twiss parameters. In 6 by 6 case, |
|   `definitions` | EquAlign["Text[Let V denote the matrix to define the normal mode, i.e.,\n\\begin{align}\n {\\bm U} = {\\rm V}{\\bm u}\\,,\n\ |
| `synchrotron-radiation` | EquAlign["Text[The evaluation of synchrotron radiation in SAD is done usingbased on ``kinematical method\":\\\\\n\\\\\nLet ${\ |
| `equilibrium-beam-envelope` | EquAlign["Text[The equilibrium beam envelope (second order moment) $\\left\\langle X_iX_j\\right\\rangle$ is calculated by th\ |
| `Bmag-coefficient` | EquAlign["Text[The $B_{\\rm mag}$ coeeficient represents the mismatch between two optics:\n\\begin{align}B_{\\rm mag}\\equiv\\ |

### elements の下位項目 (170)

| 名前 | 説明 (最初の行) |
|---|---|
| `APERT` | An aperture. Only valid in tracking. A particle with |
| `BEND` | A bending magnet. |
|   `AE1` | The absolute face angle at the entrance. The effective face angle is E1 * ANGLE + AE1, and a positive angle at the entrance co |
|   `AE2` | The absolute face-angle at the exit to the bending angle. The effective face angle is E2 * ANGLE + AE2, and a positive angle a |
|   `ANGLE` | The bending angle. If positive, it bends the orbit in x-s plane toward negative-x-direction. ANGLE determines the geometry of |
|   `DISFRIN` | If nonzero, the nonlinear Maxwellian fringe is suppressed. |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DROTATE` | Additional rotation in x-y plane to simulate a rotation error. DROTATE does not affect the geometry of the ring. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `E1` | The ratio of the face-angle at the entrance to the bending angle. The effective face angle is E1 * ANGLE + AE1, and a positive |
|   `E2` | The ratio of the face-angle at the exit to the bending angle. The effective face angle is E2 * ANGLE + AE2, and a positive ang |
|   `F1` | Length of the slope of the field at the edge as: |
|   `FB1` | F1 at the entrance. Actually F1 + FB1 is used at the entrance. |
|   `FB2` | F1 at the exit. Actually F1 + FB2 is used at the exit. |
|   `FRINGE` | When FRINGE is non-zero, the effect of the linear fringe F1 is taken into account both at the entrance and the exit. |
|   `K0` | The normal dipole magnetic field component (times the length L). |
|   `K1` | The normal quadrupole magnetic field component (times the length L). |
|   `L` | The effective length along the arc of the orbit. |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `transformation:BEND` | The transformation of a bend depends on the value of K1. If K1 is zero, it is a series of transformations: |
| `CAVI` | Accelerating structure. |
|   `DISFRIN` | If nonzero, the Maxwellian fringe is suppressed. The effects of DISFRIN and FRINGE are summarized as |
|   `DPHI` | Relative phase offset. The stable synchrotron phase above the transition is near PHI = 0. The acceleration is given by |
|   `DVOLT` | Additional accelerating voltage to be added to VOLT without affecting the design momentum. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `FREQ` | Rf frequency. If this keyword is nonzero, the keyword HARM is ignored. |
|   `HARM` | A harmonic number. This is valid only when FREQ is zero. |
|   `L` | The effective length. |
|   `PHI` | Relative phase offset. The stable synchrotron phase above the transition is near PHI = 0. The acceleration is given by |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `V02` | The y^2-dependence of the acceleration. Tracking only. |
|   `V1` | The linear x-dependence of the acceleration. Tracking only. |
|   `V11` | The xy-dependence of the acceleration. Tracking only. |
|   `VOLT` | Accelerating peak voltage in Volt. |
| `COORD` | An element for an arbitrary coordinate transformation. This element can be used to express an off-axis element. |
| `default-keyword` | The default and available non-default variable keywords are: |
| `DECA` | A decapole magnet. |
|   `DISFRIN` | If nonzero, the nonlinear Maxwellian fringe is suppressed. |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `K4` | The normal decapole magnetic field component (times the length L). |
|   `L` | The effective length. |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `transformation:THIN` | EquAlign["Text[The transformation of a {\\tt DECA} is given as\n\\begin{align}\n &\\exp(:F_{\\rm in}:)\\exp(:a {\\tt L}:)\\\ |
| `DODECA` | A dodecapole magnet. |
|   `DISFRIN` | If nonzero, the nonlinear Maxwellian fringe is suppressed. |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `K5` | The normal dodecapole magnetic field component (times the length L). |
|   `L` | The effective length. |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `transformation:THIN` | EquAlign["Text[The transformation of a {\\tt DODECA} is given as\n\\begin{align}\n &\\exp(:F_{\\rm in}:)\\exp(:a {\\tt L}:)\ |
| `DRIFT` | A drift space. |
|   `L` | The length, can be negative. |
|   `RADIUS` | Radius of the vacuum chamber. Effective when SPAC is ON. |
|   `transformation:DRIFT` | The transformation of a drift is written as |
| `keywords` | Available keywords are: |
| `MARK` | MARK elements play special roles in FFS: |
|   `OFFSET` | OFFSET is a relative position from the current position. A fraction is allowed to specify a location within an element. |
| `MULT` | A magnet with multipoles. Note that the reference plane is defined so that the skew quadrupole component becomes zero. |
|   `K0` | The normal dipole magnetic field component (times the length L). |
|   `SK0` | The skew dipole magnetic field component (times the length L). |
|   `K1` | The normal quadrupole magnetic field component (times the length L). |
|   `SK1` | The skew quadrupole magnetic field component (times the length L). |
|   `K2` | The normal sextupole magnetic field component (times the length L). |
|   `SK2` | The skew sextupole magnetic field component (times the length L). |
|   `K3` | The normal octupole magnetic field component (times the length L). |
|   `SK3` | The skew octupole magnetic field component (times the length L). |
|   `K4` | The normal decapole magnetic field component (times the length L). |
|   `SK4` | The skew decapole magnetic field component (times the length L). |
|   `K5` | The normal dodecapole magnetic field component (times the length L). |
|   `SK5` | The skew dodecapole magnetic field component (times the length L). |
|   `K6` | The normal 14-pole magnetic field component (times the length L). |
|   `SK6` | The skew 14-pole magnetic field component (times the length L). |
|   `K7` | The normal 16-pole magnetic field component (times the length L). |
|   `SK7` | The skew 16-pole magnetic field component (times the length L). |
|   `K8` | The normal 18-pole magnetic field component (times the length L). |
|   `SK8` | The skew 18-pole magnetic field component (times the length L). |
|   `K9` | The normal 20-pole magnetic field component (times the length L). |
|   `SK9` | The skew 20-pole magnetic field component (times the length L). |
|   `K10` | The normal 22-pole magnetic field component (times the length L). |
|   `SK10` | The skew 22-pole magnetic field component (times the length L). |
|   `K11` | The normal 24-pole magnetic field component (times the length L). |
|   `SK11` | The skew 24-pole magnetic field component (times the length L). |
|   `K12` | The normal 26-pole magnetic field component (times the length L). |
|   `SK12` | The skew 26-pole magnetic field component (times the length L). |
|   `K13` | The normal 28-pole magnetic field component (times the length L). |
|   `SK13` | The skew 28-pole magnetic field component (times the length L). |
|   `K14` | The normal 30-pole magnetic field component (times the length L). |
|   `SK14` | The skew 30-pole magnetic field component (times the length L). |
|   `K15` | The normal 32-pole magnetic field component (times the length L). |
|   `SK15` | The skew 32-pole magnetic field component (times the length L). |
|   `K16` | The normal 34-pole magnetic field component (times the length L). |
|   `SK16` | The skew 34-pole magnetic field component (times the length L). |
|   `K17` | The normal 36-pole magnetic field component (times the length L). |
|   `SK17` | The skew 36-pole magnetic field component (times the length L). |
|   `K18` | The normal 38-pole magnetic field component (times the length L). |
|   `SK18` | The skew 38-pole magnetic field component (times the length L). |
|   `K19` | The normal 40-pole magnetic field component (times the length L). |
|   `SK19` | The skew 40-pole magnetic field component (times the length L). |
|   `K20` | The normal 42-pole magnetic field component (times the length L). |
|   `SK20` | The skew 42-pole magnetic field component (times the length L). |
|   `K21` | The normal 44-pole magnetic field component (times the length L). |
|   `SK21` | The skew 44-pole magnetic field component (times the length L). |
|   `AE1` | The absolute face angle at the entrance. The effective face angle is E1 * ANGLE + AE1, and a positive angle at the entrance co |
|   `AE2` | The absolute face-angle at the exit to the bending angle. The effective face angle is E2 * ANGLE + AE2, and a positive angle a |
|   `ANGLE` | The bending angle. If positive, it bends the orbit in x-s plane toward negative-x-direction. ANGLE determines the geometry of |
|   `DISFRIN` | If nonzero, the nonlinear maxwellian fringe is suppressed. The effects of DISFRIN and FRINGE are summarized as |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DPHI` | Relative phase offset. The stable synchrotron phase above the transition is near PHI = 0. The acceleration is given as |
|   `DVOLT` | Additional accelerating peak voltage to be added to Volt, without affecting the design momentum p0(s). |
|   `E1` | The ratio of the face-angle at the entrance to the bending angle. The effective face angle is E1 * ANGLE + AE1, and a positive |
|   `E2` | The ratio of the face-angle at the exit to the bending angle. The effective face angle is E2 * ANGLE + AE2, and a positive ang |
|   `F1` | EquAlign["Text[{\\tt F1} and {\\tt F2} are parameters to characterize the slope of the field at the edges defined as:\n\\begin |
|   `F2` | EquAlign["Text[{\\tt F1} and {\\tt F2} are parameters to characterize the slope of the field at the edges defined as:\n\\begin |
|   `FB1` | Linear Fringe length F1 for the K0 component at the entrance. |
|   `FB2` | Linear Fringe length F1 for the K0 component at the exit. |
|   `FREQ` | Rf frequency. If this keyword is nonzero, the keyword HARM is ignored. |
|   `FRINGE` | The effects of the linear fringe (characterized by F1 and F2), and the nonlinear Mexwellian fringe are controled as: |
|   `HARM` | A harmonic number. This is valid only when FREQ is zero. |
|   `L` | The effective length. |
|   `misalignments` | Misalignments of a MULT element are expressed by the keywords DX, DY, DZ, CHI1, CHI2, and ROTATE(=CHI3). They specify all misa |
|   `multipole_with_nonzero_ANGLE` | The multipoles in MULT with nonzero ANGLE are defined by |
|   `PHI` | Relative phase offset. The stable synchrotron phase above the transition is near PHI = 0. |
|   `RADIUS` | Radius of the vacuum chamber. Effective when SPAC is ON. |
|   `VOLT` | Accelerating peak voltage in Volt. |
| `OCT` | A octupole magnet. |
|   `DISFRIN` | If nonzero, the nonlinear Maxwellian fringe is suppressed. |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `K3` | The normal octupole magnetic field component (times the length L). |
|   `L` | The effective length. |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `transformation:THIN` | EquAlign["Text[The transformation of a {\\tt OCT} is given as\n\\begin{align}\n &\\exp(:F_{\\rm in}:)\\exp(:a {\\tt L}:)\\e\ |
| `QUAD` | A quadrupole magnet. |
|   `DISFRIN` | If nonzero, the nonlinear maxwellian fringe is suppressed. The effects of DISFRIN and FRINGE are summarized as |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `F1` | EquAlign["Text[{\\tt F1} and {\\tt F2} are parameters to characterize the slope of the field at the edges defined as:\n\\begin |
|   `F2` | EquAlign["Text[{\\tt F1} and {\\tt F2} are parameters to characterize the slope of the field at the edges defined as:\n\\begin |
|   `FRINGE` | The effects of the linear fringe (characterized by F1 and F2), and the nonlinear Mexwellian fringe are controled as: |
|   `K1` | The normal quadrupole magnetic field component (times the length L). |
|   `L` | The effective length. |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `transformation:QUAD` | The transformation in a QUAD is a sequence of: |
| `SEXT` | A sextupole magnet. |
|   `DISFRIN` | If nonzero, the nonlinear Maxwellian fringe is suppressed. |
|   `DISRAD` | If nonzero, the synchrotron radiation in the particle-tracking is inhibited. |
|   `DX` | Horizontal displacement of magnet. This applied before the rotation by ROTATE. |
|   `DY` | Vertical displacement of magnet. This applied before the rotation by ROTATE. |
|   `K2` | The normal sextupole magnetic field component (times the length L). |
|   `L` | The effective length. |
|   `ROTATE` | Rotation in x-y plane. After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by -(amount given b |
|   `transformation:THIN` | EquAlign["Text[The transformation of a {\\tt SEXT} is given as\n\\begin{align}\n &\\exp(:F_{\\rm in}:)\\exp(:a {\\tt L}:)\\\ |
| `SOL` | A solenoid. Unlike other elements, SOL elements inserted at boundaries or of a solenoid or at where the field changes. Between |
|   `BOUND` | BOUND = 1 must be given at both sides of the boundaries of a solenoid, otherwise SOL only specifies the change of BZ |
|   `BZ` | The longitudinal field of a solenoid. If a SOL is used with BOUND = 0 (default), only BZ is used to change the field, and no c |
|   `DISFRIN` | EquAlign["Text[Disables the nonlinear fringe of solenoid if nonzero. The default is 0. The transformation for the nonlinear f\ |
|   `DPX` | An offset of the design ORBIT angle CHI1 relative to the solenoid axis at SOL with GEO = 1. |
|   `DPY` | An offset of the design ORBIT angle CHI2 relative to the solenoid axis at SOL with GEO = 1. |
|   `DX` | An x-offset of the design ORBIT relative to the solenoid center at SOL with GEO = 1. |
|   `DY` | A y-offset of the design ORBIT relative to the solenoid center at SOL with GEO = 1. |
|   `F1` | The length of fringe of the solenoid field. It affects only the EMITTANCE calculation. If F1 = 0 (default), no radiation arise |
|   `GEO` | One of boundaries (with GEO = 1) of a solenoid must have GEO = 1 to specify the alignment of the design orbit. At a SOL elemen |

### expression の下位項目 (46)

| 名前 | 説明 (最初の行) |
|---|---|
| `(-)` | operator for subtraction or unary minus. |
| `(/)` | operator for division. |
| `AddTo(+=)` | a+=b is equivalent to a=a+b . |
| `Alternatives(|)` | a \| b \| ... represents a pattern which matches one of patterns a, b, ... |
| `And(&&)` | a && b returns True(==1) when both a and b are nonzero real, False(==0) otherwise. b is not evaluated when a is zero. |
| `Apply (@@)` | f@@a applies function f to subexpressions of a. f@@[a,level] specifies a levelspec to apply by level. |
| `CompoundExpression(;)` | a ; b evaluates a, then evaluates b and returns its result. |
| `Decrement(--)` | a-- decrements a by 1, returning the old value of a. --a decrements a by 1, returning the new value of a. |
| `DivideBy(/=)` | a/=b is equivalent to a=a/b. |
| `Dot(.)` | a . b returns the inner product of a and b. |
| `Equal(==)` | a == b returns True(==1) if a and b are both Real with the same value or a String with same length and value. If one of a or b |
| `Function(&)` | a & is a pure-function whose argument is specified #, #n, ##, ##n. |
| `Greater(>)` | If both a and b are real, a > b returns True if a is greater than b, False otherwise. If one of a or b is a list, or both are |
| `GreaterEqual(>= or =>)` | If both a and b are real, a => b returns True if a is greater than or equal to b, False otherwise. If one of a or b is a list, |
| `Increment(++)` | a++ increments a by 1, returning the old value of a. ++a increments a by 1, returning the new value of a. |
| `Less(<)` | If both a and b are real, a < b returns True if a is less than b, False otherwise. If one of a or b is a list, or both are lis |
| `LessEqual(<= or =<)` | If both a and b are real, a <= b returns True if a is less than or equal to b, False otherwise. If one of a or b is a list, or |
| `List({})` | {a,b,c...} is a list structure. |
| `Map (/@)` | f/@a maps function f to subexpressions of a. f/@[a,level] specifies a levelspec of map by level. |
| `MapAll(//@)` | f//@a maps function f to all subexpressions of a. f//@[a,Heads->True] maps including the heads of a and its subexpressions. |
| `Member(@)` | f@a refers the member a of an instance or a class f. Otherwise it is same as f[a]. f@g@h means (f@g)@h, and f@g[h] (f@g)[h]. |
| `MessageName(::)` | symbol::tag returns a message associated with symbol and tag |
| `Not(~)` | ~a returns True(==1) when a is zero, False(==0) when a is a nonzero real. |
| `Or(||)` | a \|\| b returns True(==1) when a is nonzero real or b is nonzero real, False(==0) otherwise. b is not evaluated when a is nonzero. |
| `Part([[]])` | a[[b,..]] is a subexpression of an expression a. |
| `PatternTest(?)` | pattern?test matches to an object which matches pattern then test[object] gives True. |
| `Plus(+)` | a + b returns the sum of a and b. If one of a or b is a list, or both are lists with same size, returns a list of element-wise |
| `Power(^)` | a ^ b returns the power of a to b. If one of a or b is a list, or both are lists with same size, returns a list of element-wis |
| `Repeated(..)` | p.. matches sequence of one ore more expressions, each matching p. |
| `RepeatedNull(...)` | p... matches sequence of zero ore more expressions, each matching p. |
| `ReplaceAll(/.)` | expr/.rule replaces all subexpressions of expr using rule. |
| `ReplaceRepeated(//.)` | expr//.rule replaces all subexpressionso of expr using rule, while a replacement is performed. |
| `Rule(->)` | pattern->expr represents a rule for ReplaceAll. |
| `RuleDelayed(:>)` | pattern:>expr represents a rule for ReplaceAll, where expr is kept unevaluated until the replacement. |
| `SameQ(===)` | a === b returns True(==1) if a and b have the same type and same value, False(==0) otherwise. |
| `Sequence([])` | a[b,c,..] means a list of b, c,.. with the head a. It is evaluated as a function-reference when a is a function or a defined-f |
| `Set(=)` | a = value sets the value b to the symbol a. {a,b,..}={v1,v2,..} sets a,b,c simultaneously. a[b,c,..]=v1 sets the part of a[b,c |
| `SetDelayed(:=)` | same as Set but the right hand side is not evaluated when it is set. |
| `StringJoin (//)` | a // b converts a and b to character-strings, then join them. |
| `SubtractFrom(-=)` | a-=b is equivalent to a=a-b . |
| `TagSet(/:)` | symb/:lhs = rhs sets rhs to lhs, associated with symbol symb. symb/:lhs := rhs sets rhs to lhs unevaluated, associated with sy |
| `Times(*)` | a * b returns the product of a and b. If one of a or b is a list, or both are lists with same size, returns a list of element- |
| `TimesBy(*=)` | a*=b is equivalent to a=a*b . |
| `Unequal(<>)` | a <> b returns True(==1) if a and b are both Real with different values or a String with different length or values. If one of |
| `UnsameQ(<=>)` | a <=> b returns True(==1) if a and b have the different types or different values. False(==0) otherwise. |
| `Unset(=.)` | a=. clears the definition assigned to a. |

### flags の下位項目 (60)

| 名前 | 説明 (最初の行) |
|---|---|
| `ABSW` | ABSW or NORELW sets the weights of variable elements independent from their values in the matching. Otherwise they are weighte |
| `BIPOL` | BIPOL or NOUNIPOL allows the change of sign of the value of the element during the matching. It affects the default keywords o |
| `CALC4D` | If CALC4D is on, the optics calculation in CALCULATE(CAL) and GO performs a 4x5 calculation (4D + dispersion). The antonyms is |
| `CALC6D` | If CALC6D is on, the optics calculation in CALCULATE(CAL) and GO performs full 6D calculation, which may takes RADCOD into acc |
| `CALEXP` | The default behavior on EXPAND for CALCULATE(CAL) and GO, when unspecified. If it is ON, the values of elements are expanded t |
| `CELL` | CELL or NOINS sets the periodic condition in calculating the optical-functions. |
| `CMPLOT` | CMPLOT enables the shift of the center of mass at the beginning of the tracking. This is almost obsolete. |
| `COD` | COD turns on finding the closed-orbit in the emittance calculation. Accessible in MAIN level. |
| `CODPLOT` | CODPLOT lets the emittance calculation return the information on the closed-orbit, extended Twiss parameters, and the beam siz |
| `CONV` | CONV is a flag set by the CALCULATE(CAL) or GO commands. It becomes True when MatchingResidual is less than CONVERGENCE. |
| `CONVCASE` | When CONVCASE is on, FFS command line parser converts the input characters to the upper case.(Default on) CONVCASE actions for |
| `DAMPONLY` | DAMPONLY is the antonym of FLUC. |
| `DAPERT` | DAPERT enables the DAPERT procedure in the multi-turn tracking to obtain the dynamic aperture diagram. Accessible in the MAIN |
| `DIFFRES` | In the case of CELL, if the total tune deviates across a difference resonance, it is counted as unstable. The default is NODIFFRES. |
| `ECHO` | ECHO enables the echo of the main input in the MAIN level. |
| `EMIOUT` | EMIOUT turns on the extended output of emittance calculation. |
| `FFSPRMPT` | When FFSPRMPT is off(default) the input prompt is In[n]:= , where n is $Line+1. Otherwise the prompt is the traditional FFS pr |
| `FIXSEED` | FIXSEED or NOMOVESEED disables the change of the seed of the random number generator after the particle tracking. |
| `FLUC` | FLUC or NODAMPONLY enables the diffusion due to synchrotron radiation in the particle tracking. Otherwise only the damping is |
| `GAUSS` | GAUSS or NOUNIFORM sets the momentum distribution of the incoming beam to be Gaussian, otherwise uniform(square) distribution |
| `GEOCAL` | When GEOCAL is on(default), the geometry of the beam line is always updated by CALCULATE(CAL) or GO commands using the current |
| `GEOFIX` | GEOFIX is the antonym of GEOCAL. |
| `HALFRES` | In the case of CELL, if the total tune deviates across a half-integer resonance, or the traces is less than -2, it is counted |
| `IDEAL` | IDEAL or NOREAL inhibits to use the component-specific deviations in the optics calculation. |
| `INS` | INS is the antonym of CELL. |
| `INTRA` | INTRA turns on the calculation of intra-beam scattering in the emittance calculation. Accessible in MAIN level. |
| `INTRES` | In the case of CELL, if the total tune deviates across an integer resonance, or the trace is larger than 2, it is counted as u |
| `JITTER` | JITTER or NOQUIET allows jitter of the center-of-mass of the incoming beam in the case of TRPT. Otherwise the center-of-mass s |
| `KEEPEXP` | If it is ON, CALCULATE(CAL) preserves the behavior on EXPAND of the last CALCULATE(CAL) or GO. The default is ON. |
| `LOG` | LOG enables the echo of all inputs in the MAIN level. Also suppresses the prompt in FFS. |
| `LOSSMAP` | If LOSSMAP is on, TrackParticles returns the location and turn where the particle lost is detected, in the 8/10th and 9/11th e |
| `LWAKE` | LWAKE turns on optics calculation with Longitudinal WakeFunction |
| `MOVESEED` | MOVESEED is the antonym of FIXSEED. |
| `PHOTONS` | When PHOTONS is ON (default is OFF), with RAD and FLUC, TrackParticles generates a list of all photons radiated through the tr |
| `POL` | POL turns on spin tracking and the calculation of equilibrium polarization in EMIT/Emittance[]. |
| `PRSVCASE` | When PRSVCASE is on, FFS command line parser preserves the input characters for the element names and patterns.(Default off) |
| `PSPAC` | When on, performs space-charge simulation in a "Particle-In-Cell" method. PSPAC is effective in tracking only. |
| `QUIET` | QUIET is the antonym of JITEER. |
| `RAD` | RAD turns on the synchrotron radiation in the particle-tracking. Accessible in the MAIN level. |
| `RADCOD` | RADCOD turns on the energy loss due to synchrotron radiation at the closed-orbit in the emittance calculation. Also turns off |
| `RADLIGHT` | When RADLIGHT is on, the function TrackParticles returns a list of trajectories which are used to calculate the synchrotron ra |
| `RADPOL` | If POL is on, RADPOL turns on Sokolov-Ternov effect in tracking. |
| `RADTAPER` | Scales all magnets except for solenoids according to the local momentum (DDP) of the closed orbit. It uses the average of DDPs |
| `REAL` | REAL is the antonym of IDEAL. |
| `RELW` | RELW is the antonym of ABSW. |
| `RFSW` | RFSW turns on the acceleration by CAVI and TCAVI element in the particle-tracking and the emittance calculation. Accessible in |
| `RING` | RING is the antonym of TRPT. |
| `SELFCOD` | When on with WSPAC, in tracking, the space charge force is calculated relative to the center of mass of the current set of par |
| `SORG` | SORG sets the origin of S (design orbit length) at the location set by ORG. |
| `SPAC` | When SPAC is on, tracking is done with space charge effect. The actual number of particles in the beam and the number of macro |
| `STABLE` | STABLE is a flag set by the CALCULATE(CAL) or GO commands in the case of CELL. It becomes True when the closed orbit is found |
| `SUMRES` | In the case of CELL, if the total tune deviates across a sum resonance, it is counted as unstable. The default is SUMRES. |
| `SUS` | If OFF, SUSPEND(SUSP) and END commands have no action. The default is ON. |
| `TRPT` | TRPT or NORING declares that the beam line is a transport line, not a part of a storage ring. The nominal momentum be changed |
| `TWAKE` | TWAKE turns on optics calculation with Transverse WakeFunction |
| `UNIFORM` | UNIFORM is the antonym of GAUSS. It assumes the momentum distribution to be a uniform(square) within +-DP. |
| `UNIPOL` | UNIPOL is the antonym of BIPOL. |
| `UNSTABLE` | UNSTABLE is the antonym of STABLE. |
| `WAKEOPT` | WAKEOPT enables an optics calculation with LWAKE or TWAKE under TRPT. The default is OFF. |
| `WSPAC` | When on, performs space-charge simulation in a "strong-weak" mode. The beam size through the beam line is to be calculated by |

### geometric-functions (本文)

```
Available geometric-functions are:

GX      geometrical coordinate xi
GY      geometrical coordinate eta
GZ      geometrical coordinate zeta
CHI1    geometrical rotation angle ch1_1
CHI2    geometrical rotation angle ch1_2
CHI3    geometrical rotation angle ch1_3

The geometrical coordinate {xi, eta, zeta} is set by the ORG command, and its default origin is at the entrance of the beam line, a
nd the default directions are xi in s-direction, eta in -(x-direction), and zeta in -(y-direction) at the entrance. The set {xi, et
a, zeta} forms a right-hand system.
   The rotation angles are defined so as to give the local {x,y,s} is written

    {x, y, s}_local
   = rotate[s, chi3] rotate[x, chi2] rotate[y, -chi1]{x, y, s}_origin,

where rotate[a, b] reads "rotate around the new-a vector by b right-handedly.

geometric functions disignate the geometry of the coordinate. If the geometry of orbit is needed, use LINE["OGEO"], etc., or DISP O
G.
```

### machine-error-commands (本文)

```
Usage: machine-error-command [options] amount component-pattern ..

where machine-error-command is one of

command    keyword affected     applicable types
DELX       DX                   BEND QUAD SEXT OCT DECA DODECA SOL CAV
DELY       DY                   BEND QUAD SEXT OCT DECA DODECA SOL CAV
DL         L                    DRIFT SOL
DTHETA     ROTATE               QUAD SEXT OCT DECA DODECA CAV
DTHETA     DROTATE              BEND
DK         default-keyword      DRIFT BEND QUAD SEXT OCT DECA DODECA MULT SOL CAV
DDK        K0 or DBZ            BEND SOL

   amount is the amount of the error,
   component-pattern is the pattern to specify the components to be applied.

   Options are	

RANDOM(R)        Set amount*GaussRandom[] to the keyword. 
UNIFORM(U)       Set the specified amount to the keyword without random number.
INCOHERENT(INC)  GaussRandom[] is called for each component. Default.
COHERENT(C)      GaussRandom[] is called once for each component-pattern.
PUT(P)           Set the error to the keyword. Default.
ADD(A)           Add the error to the keyword.	
```

### matching-function-commands (本文)

```
Usage: (1) matching-function  goal-value [off-momentum-points]
       (2) matching-functionM goal-value [off-momentum-points]
       (3) matching-functionI incoming-value
       (4) matching-functionSCALE scale

(1) sets the matching condition for matching-function at the current fitting point or region with the goal-value and the off-moment
um-points (see off-momentum-matching).
   If off-momentum-points is omitted, the previous value for this matching-function at this fitting-point is assumed. If the previo
us value is not defined, 1 is assumed. If -1 is given for off-momentum-points, the matching-function is rejected from the matching 
(see REJECT(REJ)).
   If "*" is given for goal-value, the previous value is used if exists.

Example:   BX 10 3    (beta_x to be 10 at 3 momenta) 
           BX 20      (now beta_x to be 20 at 3 momenta (previous setting))
           BX *  5    (now beta_x to be 20 (previous setting) at 5 momenta)

(2) If the letter "M" is appended to matching-function, it means the maximum-fitting for the function. The maximum of either the va
lue (for positive-definite functions) or the absolute value (for bipolar functions) are to be limited in the matching.

(3) If SCALE is appended to matching-function, it sets the scale of the input/output of the function to scale. This scale is used i
n the matching-function commands, DISPLAY(DISP), SHOW, etc.

(4) If the current fit location is at a MARK, @ refers the saved value at the MARK. @- refers -(saved value), and @* refers (the or
ientation of MARK)*(saved value).  These are useful to match using the already-saved values, for instance to match between two beam
 lines having the common MARK.

   Available matching-functions are:
optical-functions (see optical-functions):
AX BX GMX NX AY BY GMY NY EX EPX EY EPY R1 R2 R3 R4 DETR DX DPX DY DPY DZ DDP PEX PEPX PEY PEPY TRX TRY LENG
geometric-functions (see geometric-functions):
GX GY GZ CHI1 CHI2 CHI3
```

### multi-turn-tracking (本文)

```
The multi-turn tracking can be done by TrackParticles or DynamicApertureSurvey[] function in FFS.

The multi-turn tracking uses the closed orbit, normal coordinate, and the equilibrium emittances. Therefore One of EMITTANCE(EMIT),
 the Emittance[] function, or the EMIT command in the MAIN level are necessary to be done once in prior to the multi-turn-tracking.
 The values of EMITX, EMITY, EMITZ, SIGE can be changed between EMITs and the multi-turn-tracking.
```

### off-momentum-matching (本文)

```
FFS matches the optical functions for an orbit with finite momentum deviation.

Example: 
   DP=0.01;      sets the range of momentum to DP0-0.01 < dp/p0 < DP0+0.01 .
   BX 10 9;      sets the goal of betax to 10 m, at 9 points.
                 in the range above, i.e.,
                 dp/p0 = {-0.01,-0.0075,-0.005,-0.0025,0,
                          0.0025,0.005,0.0075,0.01} + DP0 .
   GO starts the matching.

As this example, the off-momentum points are chosen with equal separation. If the off-momentum point n is an even number, the chose
n points are same as the case of n+1, but the on-momentum point (==DP0) is excluded.
 If a list {min, max} is set to DP, the off-momentum matching is performed at dp/p(k) = min*Sin[k Pi/2/n], max*Sin[k Pi/2/n], {k,0,
fitp/2}, where fitp is the number of off-momentum points given by FITP.
```

### optical-functions (本文)

```
Available optical functions for matching are:	

AX      alpha_X
BX      beta_X
GMX     gamma_X
NX      psi_X, the default scale is 1/(2Pi)
AY      alpha_Y
BY      beta_Y
GMY     gamma_Y
NY      psi_Y, the default scale is 1/(2Pi)
EX      eta_X   (dispersion_X)
EPX     eta_Px  (dispersion_PX) 
EY      eta_Y   (dispersion_Y)
EPY     eta_Py  (dispersion_PY)
R1      R_1     (see x-y-coupling)
R2      R_2     (see x-y-coupling)
R3      R_3     (see x-y-coupling)
R4      R_4     (see x-y-coupling)
DETR    R_1*R_4 - R_2*R_3 (see x-y-coupling)
DX      dx
DPX     dpx
DY      dy
DPY     dpy
DZ      dz
DDP     delta=dp/p0
AZ      alpha_Z
BZ      beta_Z
GMZ     gamma_Z
NZ      psi_Z, the default scale is 1/(2Pi)
ZX      zeta_X  (z-dispersion_X)
ZPX     zeta_Px (z-dispersion_PX) 
ZY      zeta_Y  (z-dispersion_Y)
ZPY     zeta_Py (z-dispersion_PY)
PEX     eta_x   (dispersion_x)
PEPX    eta_px  (dispersion_px) 
PEY     eta_y   (dispersion_y)
PEPY    eta_yy  (dispersion_py)
TRX     trace(T_X), only defined at the end of the beam line.
TRY     trace(T_Y), only defined at the end of the beam line.
BMAGX   BmagX, against the reference optics.
BMAGY   BmagY, against the reference optics.
BMAGZ   BmagZ, against the reference optics.
LENG    length of the design orbit	

In the above, upper case X, Px, Y, Py represents the x-y decoupled coordinate. EX, EPX, EY, EPY refer the decoupled coordinate, whi
le PEX, PEPX, PEY, PEPY are in the physical coordinate. On the other hand,  DX, DPX, DY, DPY refer the physical coordinate.
```

### pattern の下位項目 (1)

| 名前 | 説明 (最初の行) |
|---|---|
| `MatchQ` | MatchQ[x, pat] returns True if x matches pat |

### set-value-of-element の下位項目 (2)

| 名前 | 説明 (最初の行) |
|---|---|
| `keywords` | Available keywords are: |
| `default-keyword` | The default and available non-default variable keywords are: |

### wildcards (本文)

```
   Many commands and functions accept the wildcards as a specification for the name of elements or components. The valid wildc
ards are:

       *         matches any zero or more characters.
       %         matches one character.
       {..}      matches any character enclosed.
       {^..}     matches any character not enclosed.
       ..|..     alternative pattern.
```

