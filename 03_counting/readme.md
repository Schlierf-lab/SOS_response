## Single-molecule Counting

This folder contains scripts to perform single-molecule counting, a modified approach after Lee *et al.* (2012). Counting single photoactivatable fluorescent molecules by photoactivated localization microscopy (PALM). *PNAS* 109, 17436–17441. First, `S1_swift_Pbleach0p0001_Pblink0.m` is used to run swift with minimal parameter assumptions. Scripts S1b-S1f are used to estimate photophysical parameters of the FP; this needs to be done only once. Scripts S2-S4 are then ran to perform single-molecule counting with the determined parameters on individual localization files.

<sub>-- Andreas Hartmann</sub>
<sub>-- Leonard Schärfen</sub>
