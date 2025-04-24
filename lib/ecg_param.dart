import 'ieee_remainder.dart';

class EcgParam {
  //////////////////////////////////////////////////////////////////
  //  GLOBAL ECG PARAMETERS:
  //////////////////////////////////////////////////////////////////
  late int _n; // Number of heart beats
  late double _hrstd; // Heart rate std
  late double _hrmean; // Heart rate mean
  late double _lfhfratio; // LF/HF ratio
  late int _sfecg; // ECG sampling frequency
  late int _sf; // Internal sampling frequency
  late double _amplitude; // Amplitude for the plot area
  late int _seed; // Seed
  late double _anoise; // Amplitude of additive uniform noise
  late int _period;

  // Define frequency parameters for rr process
  // flo and fhi correspond to the Mayer waves and respiratory rate respectively
  late double _flo; // Low frequency
  late double _fhi; // High frequency
  late double _flostd; // Low frequency std
  late double _fhistd; // High frequency std

  // Order of extrema: [P Q R S T]
  late List<double> _theta; // ti not in radians
  late List<double> _a;
  late List<double> _b;

  // Animation variables
  late int _ecgAnimateInterval;

  // Flag to know if all parameters are valid
  late bool _allParametersValid;

  //////////////////////////////////////////////////////////////////
  // Creates a new instance of EcgParam
  //////////////////////////////////////////////////////////////////
  EcgParam() {
    resetParameters();
  }

  //////////////////////////////////////////////////////////////////
  // GLOBAL Set/Get Parameter Functions:
  //////////////////////////////////////////////////////////////////
  void setN(int value) {
    _n = value;
    _allParametersValid = false;
  }

  int getN() {
    return _n;
  }

  void setHrStd(double value) {
    _hrstd = value;
    _allParametersValid = false;
  }

  double getHrStd() {
    return _hrstd;
  }

  void setHrMean(double value) {
    _hrmean = value;
    _allParametersValid = false;
  }

  double getHrMean() {
    return _hrmean;
  }

  void setLfHfRatio(double value) {
    _lfhfratio = value;
    _allParametersValid = false;
  }

  double getLfHfRatio() {
    return _lfhfratio;
  }

  void setSfEcg(int value) {
    _sfecg = value;
    _ecgAnimateInterval = (1000 ~/ _sfecg).toInt();
    _allParametersValid = false;
  }

  int getSfEcg() {
    return _sfecg;
  }

  void setSf(int value) {
    _sf = value;
    _allParametersValid = false;
  }

  int getSf() {
    return _sf;
  }

  void setAmplitude(double value) {
    _amplitude = value;
    _allParametersValid = false;
  }

  double getAmplitude() {
    return _amplitude;
  }

  void setSeed(int value) {
    _seed = value;
    _allParametersValid = false;
  }

  int getSeed() {
    return _seed;
  }

  void setANoise(double value) {
    _anoise = value;
    _allParametersValid = false;
  }

  double getANoise() {
    return _anoise;
  }

  void setPeriod(int value) {
    _period = value;
    _allParametersValid = false;
  }

  int getPeriod() {
    return _period;
  }

  void setFLo(double value) {
    _flo = value;
    _allParametersValid = false;
  }

  double getFLo() {
    return _flo;
  }

  void setFHi(double value) {
    _fhi = value;
    _allParametersValid = false;
  }

  double getFHi() {
    return _fhi;
  }

  void setFLoStd(double value) {
    _flostd = value;
    _allParametersValid = false;
  }

  double getFLoStd() {
    return _flostd;
  }

  void setFHiStd(double value) {
    _fhistd = value;
    _allParametersValid = false;
  }

  double getFHiStd() {
    return _fhistd;
  }

  void setTheta(int index, double value) {
    _theta[index] = value;
    _allParametersValid = false;
  }

  double getTheta(int index) {
    return _theta[index];
  }

  void setA(int index, double value) {
    _a[index] = value;
    _allParametersValid = false;
  }

  double getA(int index) {
    return _a[index];
  }

  void setB(int index, double value) {
    _b[index] = value;
    _allParametersValid = false;
  }

  double getB(int index) {
    return _b[index];
  }

  void setEcgAnimateInterval(int value) {
    _ecgAnimateInterval = value;
    _allParametersValid = false;
  }

  int getEcgAnimateInterval() {
    return _ecgAnimateInterval;
  }

  //////////////////////////////////////////////////////////////////
  // Functions:
  //////////////////////////////////////////////////////////////////
  /*
   * Check to see if all parameters are valid
   */
  bool isValid() {
    return _allParametersValid;
  }

  /*
   * In this function, it can be enforced
   * additional rules.
   */
  bool checkParameters() {
    bool retValue = true;
    _allParametersValid = true;

    // Check the Internal frequency respect to ECG frequency
    if (ieeeRemainder(_sf.toDouble(), _sfecg.toDouble()).toInt() != 0){
    retValue = false;
      _allParametersValid = false;
    }

    return retValue;
  }


  /*
   * ReInit the Button Parameters' values
   */
  void resetParameters() {
    /* General Interface parameters */
    _n = 256;
    _sfecg = 256;
    _sf = 512;
    _anoise = 0.04;  //  0.1;
    _hrmean = 60.0;
    _hrstd = 1.0;
    _seed = 1;
    _amplitude = 1.4;
    _flo = 0.1;
    _fhi = 0.25;
    _flostd = 0.01;
    _fhistd = 0.01;
    _lfhfratio = 0.5;

    /* ECG morphology: Order of extrema: [P Q R S T] */
    _theta = [-60.0, -15.0, 0.0, 15.0, 90.0];
    _a = [1.2, -5.0, 30.0, -7.5, 0.75];
    _b = [0.25, 0.1, 0.1, 0.1, 0.4];

    _ecgAnimateInterval = (1000 ~/ _sfecg).toInt();

    _allParametersValid = true;
  }
}
