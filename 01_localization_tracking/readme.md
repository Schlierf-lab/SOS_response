## Single-molecule Localization and Tracking

[SMAP](https://github.com/jries/SMAP) is used on the raw movie files (.nd2) with parameters specified in Supplementary Table S2. Cell outlines are determined from brightfield images using [oufti](https://oufti.org/). For PALM movies with fiducial markers, the `Drift_Correction/drift_correction.m` script is used to shift localizations. The localization table is then converted using the `SMAP2swift.m` script, using SMAP and oufti files as input. Next, the tracking software swift (Endesfelder, M., Schießl, C., Turkowyd, B., Lechner, T., and Endesfelder, U., 2020, *Manuscript in prep.*) is called.

<sub>-- Leonard Schärfen</sub>
<sub>-- Andreas Hartmann</sub>
