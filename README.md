# ECG waveform generator on dart

Ported version of ECGSYN of ECG waveform generator

## Description
__ECGSYN__ is a library for generating synthetic __ECG__ signals developed by _Patrick McSharry_ and _Gary Clifford_ in 2003. It allows you to specify parameters such as mean heart rate, beat count, sampling frequency, and wave morphology (P, Q, R, S, T).
The source code in __C__ is available on the __PhysioNet__ platform as an archive at https://physionet.org/files/ecgsyn/1.0.0/ecgsyn.tar.gz. The archive size is 161.4 KB, last updated on April 12, 2019.
The __ecgsyn.c__ file is the main source of the code, which can be found in the __C__ directory inside the archive. Additional files __dfour1.c__ and __ran1.c__ from _Numerical Recipes in C_ are required for compilation, which aren't included in the archive.

In addition to the __C__ version, the archive contains a __Java__ version in the __java__ folder.

Finally, a separate Java version of the __Java ECG Generator__ is available on the __MIT__ website at http://www.mit.edu/~gari/CODE/ECGSYN/JAVA/APPLET2/ecgsyn/ecg-java/source.html, where can download the source code. __It was this project that became the prototype for the dart version__.

## Notes
* Those who wish can follow the link https://www.physionet.org/content/ecgsyn/1.0.0/ and get information about the package __ECGSYN__ first-hand.
* Also, for those who want to refresh their knowledge in the field of computational mathematics, a pdf copy of the above-mentioned book __"Numerical Recipes in C"__ is attached to the project.

## Porting
After a day spent trying to port the C version, finding and porting additional functions __ran1.c__ and __dfour1.c__ using the additional __fftea__ dart package to implement the Fast Fourier Transform (FFT), I realized that this was a dead end. I was able to generate an ECG signal, but I realized that this is far from a universal standalone solution.

As an alternative, I tried using the previously mentioned standalone __MIT__ Java version. The only external function __Math.IEEEremainder(double x, double y)__ was ported and included into the project without much effort.

## Results
__ecg_calc__ is a console __dart__ application that generates a text file ecgsyn.dat containing a time series of a synthetic ECG signal. It consists of three columns separated by a space:
* Timestamp in seconds, starting from 0, with a step determined by the sampling frequency (sfecg). For example, for sfecg=256 Hz, the step is 1/256 ≈ 0.00390625 s.
* ECG voltage value in millivolts (mV), corresponding to the signal amplitude at each point.
* ECG signal peaks
  
__The file contains a fixed number of lines: 64K__.
