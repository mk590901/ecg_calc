# ECG waveform generator on dart

Ported version of ECGSYN of ECG waveform generator

## Description
__ECGSYN__ is a library for generating synthetic __ECG__ signals developed by _Patrick McSharry_ and _Gary Clifford_ in 2003. It allows you to specify parameters such as mean heart rate, beat count, sampling frequency, and wave morphology (P, Q, R, S, T).
The source code in __C__ is available on the __PhysioNet__ platform as an archive at https://physionet.org/files/ecgsyn/1.0.0/ecgsyn.tar.gz. The archive size is 161.4 KB, last updated on April 12, 2019.
The __ecgsyn.c__ file is the main source of the code, which can be found in the __C__ directory inside the archive. Additional files __dfour1.c__ and __ran1.c__from _Numerical Recipes in C_ are required for compilation, which aren't included in the archive.

## Notes
* Those who wish can follow the link https://www.physionet.org/content/ecgsyn/1.0.0/ and get information about the package __ECGSYN__ first-hand.
* Also, for those who want to refresh their knowledge in the field of computational mathematics, a pdf copy of the above-mentioned book __"Numerical Recipes in C"__ is attached to the project.

