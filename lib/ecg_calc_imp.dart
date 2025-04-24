import 'dart:math';
import 'dart:io';
import 'ecg_log.dart';
import 'ecg_param.dart';
import 'ieee_remainder.dart';

class EcgCalc {

  final String outfile = "ecgsyn.dat";
  //  Creates a new instance of EcgCalc
  EcgCalc(EcgParam parameters, EcgLogWindow logOb) {
    paramOb = parameters;
    ecgLog = logOb;

    /* variables for static function ran1() */
    iy = 0;
    iv = List<int>.filled(NTAB, 0);
  }

  bool calculateEcg() {
    bool retValue;

    ecgLog.println("Starting to calculate ECG....");

    retValue = dorun();

    ecgLog.println("Finished calculating ECG table data.\n");

    return retValue;
  }

  int getEcgResultNumRows() {
    return ecgResultNumRows;
  }

  double getEcgResultTime(int index) {
    return ecgResultTime[index];
  }

  double getEcgResultVoltage(int index) {
    return ecgResultVoltage[index];
  }

  int getEcgResultPeak(int index) {
    return ecgResultPeak[index];
  }

  /* C defines */
  static const double PI = 3.141592653589793; //2.0 * asin(1.0);
  static const int NR_END = 1;
  static const int IA = 16807;
  static const int IM = 2147483647;
  static const double AM = (1.0 / IM);
  static const int IQ = 127773;
  static const int IR = 2836;
  static const int NTAB = 32;
  static const double NDIV = (1 + (IM - 1) / NTAB);
  static const double EPS = 1.2e-7;
  static const double RNMX = (1.0 - EPS);

  ///////////////////////////////////////////////////////////////////////////////
  //  DEFINE PARAMETERS AS GLOBAL VARIABLES
  ///////////////////////////////////////////////////////////////////////////////
  // Order of extrema: [P Q R S T]
  List<double> ti = List<double>.filled(6, 0.0); // ti converted in radians
  List<double> ai = List<double>.filled(6, 0.0); // new calculated a
  List<double> bi = List<double>.filled(6, 0.0); // new calculated b

  int necg = 0; // Number of ECG outputs
  int mstate = 3; // System state space dimension
  double xinitial = 1.0; // Initial x co-ordinate value
  double yinitial = 0.0; // Initial y co-ordinate value
  double zinitial = 0.04; // Initial z co-ordinate value
  int rseed = 0;
  double h = 0.0;
  List<double> rr = [];
  List<double> rrpc = [];

  /*
   * Variables for static function ran1()
   */
  int iy = 0;
  List<int> iv = [];

  /*
   * ECG Result Variables
   */
  /* Result Vectors*/
  List<double> ecgResultTime = [];
  List<double> ecgResultVoltage = [];
  List<int> ecgResultPeak = [];

  int ecgResultNumRows = 0;

  /* Object Variables */
  late EcgParam paramOb;
  late EcgLogWindow ecgLog;

  /*--------------------------------------------------------------------------*/
  /*    UNIFORM DEVIATES                                                      */
  /*--------------------------------------------------------------------------*/

  double ran1() {
    int j;
    int k;
    double temp;
    bool flg;

    if (iy == 0) {
      flg = false;
    } else {
      flg = true;
    }

    if (rseed <= 0 || !flg) {
      if (-rseed < 1) {
        rseed = 1;
      } else {
        rseed = -rseed;
      }

      for (j = NTAB + 7; j >= 0; j--) {
        k = (rseed ~/ IQ);
        rseed = IA * (rseed - k * IQ) - IR * k;
        if (rseed < 0) rseed += IM;
        if (j < NTAB) iv[j] = rseed;
      }
      iy = iv[0];
    }

    k = (rseed ~/ IQ);
    rseed = IA * (rseed - k * IQ) - IR * k;
    if (rseed < 0) rseed += IM;

    j = (iy ~/ NDIV).toInt();
    iy = iv[j];
    iv[j] = rseed;

    if ((temp = AM * iy) > RNMX) {
      return RNMX;
    } else {
      return temp;
    }
  }

  /*
   * FFT
   */
  void ifft(List<double> data, int nn, int isign) {
    int n, mmax, m, istep, i, j;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    double swap;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
      if (j > i) {
        swap = data[j];
        data[j] = data[i];
        data[i] = swap;
        swap = data[j + 1];
        data[j + 1] = data[i + 1];
        data[i + 1] = swap;
      }
      m = n >> 1;
      while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
    mmax = 2;
    while (n > mmax) {
      istep = mmax << 1;
      theta = isign * (6.28318530717959 / mmax);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (m = 1; m < mmax; m += 2) {
        for (i = m; i <= n; i += istep) {
          j = i + mmax;
          tempr = wr * data[j] - wi * data[j + 1];
          tempi = wr * data[j + 1] + wi * data[j];
          data[j] = data[i] - tempr;
          data[j + 1] = data[i + 1] - tempi;
          data[i] += tempr;
          data[i + 1] += tempi;
        }
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }
      mmax = istep;
    }
  }

  /*
   * STANDARD DEVIATION CALCULATOR
   */
  /* n-by-1 vector, calculate standard deviation */
  double stdev(List<double> x, int n) {
    int j;
    double add, mean, diff, total;

    add = 0.0;
    for (j = 1; j <= n; j++) {
      add += x[j];
    }

    mean = add / n;

    total = 0.0;
    for (j = 1; j <= n; j++) {
      diff = x[j] - mean;
      total += diff * diff;
    }
    return sqrt(total / (n - 1));
  }

  /*
   * THE ANGULAR FREQUENCY
   */
  double angfreq(double t) {
    int i = 1 + (t / h).floor();
    return (2.0 * PI / rrpc[i]);
  }

  /*--------------------------------------------------------------------------*/
  /*    THE EXACT NONLINEAR DERIVATIVES                                       */
  /*--------------------------------------------------------------------------*/
  void derivspqrst(double t0, List<double> x, List<double> dxdt) {
    int i, k;
    double a0, w0, r0, x0, y0;
    double t, dt, dt2, zbase;
    List<double> xi = List<double>.filled(6, 0.0);
    List<double> yi = List<double>.filled(6, 0.0);

    k = 5;
    w0 = angfreq(t0);
    r0 = 1.0;
    x0 = 0.0;
    y0 = 0.0;
    a0 = 1.0 - sqrt((x[1] - x0) * (x[1] - x0) + (x[2] - y0) * (x[2] - y0)) / r0;

    for (i = 1; i <= k; i++) {
      xi[i] = cos(ti[i]);
    }
    for (i = 1; i <= k; i++) {
      yi[i] = sin(ti[i]);
    }

    zbase = 0.005 * sin(2.0 * PI * paramOb.getFHi() * t0);

    t = atan2(x[2], x[1]);
    dxdt[1] = a0 * (x[1] - x0) - w0 * (x[2] - y0);
    dxdt[2] = a0 * (x[2] - y0) + w0 * (x[1] - x0);
    dxdt[3] = 0.0;

    for (i = 1; i <= k; i++) {
      dt = ieeeRemainder(t - ti[i], 2.0 * PI);
      dt2 = dt * dt;
      dxdt[3] += -ai[i] * dt * exp(-0.5 * dt2 / (bi[i] * bi[i]));
    }
    dxdt[3] += -1.0 * (x[3] - zbase);
  }

  /*
   * RUNGA-KUTTA FOURTH ORDER INTEGRATION
   */
  void rk4(List<double> y, int n, double x, double h, List<double> yout) {
    int i;
    double xh, hh, h6;
    List<double> dydx = List<double>.filled(n + 1, 0.0);
    List<double> dym = List<double>.filled(n + 1, 0.0);
    List<double> dyt = List<double>.filled(n + 1, 0.0);
    List<double> yt = List<double>.filled(n + 1, 0.0);

    hh = h * 0.5;
    h6 = h / 6.0;
    xh = x + hh;

    derivspqrst(x, y, dydx);
    for (i = 1; i <= n; i++) {
      yt[i] = y[i] + hh * dydx[i];
    }

    derivspqrst(xh, yt, dyt);
    for (i = 1; i <= n; i++) {
      yt[i] = y[i] + hh * dyt[i];
    }

    derivspqrst(xh, yt, dym);
    for (i = 1; i <= n; i++) {
      yt[i] = y[i] + h * dym[i];
      dym[i] += dyt[i];
    }

    derivspqrst(x + h, yt, dyt);
    for (i = 1; i <= n; i++) {
      yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    }
  }

  /*
   * GENERATE RR PROCESS
   */
  void rrprocess(List<double> rr, double flo, double fhi, double flostd, double fhistd, double lfhfratio, double hrmean, double hrstd, double sf, int n) {
    int i;
    double c1, c2, w1, w2, sig1, sig2, rrmean, rrstd, xstd, ratio;
    double df;
    List<double> w = List<double>.filled(n + 1, 0.0);
    List<double> hw = List<double>.filled(n + 1, 0.0);
    List<double> sw = List<double>.filled(n + 1, 0.0);
    List<double> ph0 = List<double>.filled((n ~/ 2 - 1 + 1), 0.0);
    List<double> ph = List<double>.filled(n + 1, 0.0);
    List<double> swC = List<double>.filled((2 * n) + 1, 0.0);

    w1 = 2.0 * PI * flo;
    w2 = 2.0 * PI * fhi;
    c1 = 2.0 * PI * flostd;
    c2 = 2.0 * PI * fhistd;
    sig2 = 1.0;
    sig1 = lfhfratio;
    rrmean = 60.0 / hrmean;
    rrstd = 60.0 * hrstd / (hrmean * hrmean);

    df = sf / n;
    for (i = 1; i <= n; i++) {
      w[i] = (i - 1) * 2.0 * PI * df;
    }

    for (i = 1; i <= n; i++) {
      hw[i] = (sig1 * exp(-0.5 * pow(w[i] - w1, 2) / pow(c1, 2)) / sqrt(2 * PI * c1 * c1)) + (sig2 * exp(-0.5 * pow(w[i] - w2, 2) / pow(c2, 2)) / sqrt(2 * PI * c2 * c2));
    }

    for (i = 1; i <= n ~/ 2; i++) {
      sw[i] = (sf / 2.0) * sqrt(hw[i]);
    }

    for (i = n ~/ 2 + 1; i <= n; i++) {
      sw[i] = (sf / 2.0) * sqrt(hw[n - i + 1]);
    }

    // Randomize the phases
    for (i = 1; i <= n ~/ 2 - 1; i++) {
      ph0[i] = 2.0 * PI * ran1();
    }

    ph[1] = 0.0;
    for (i = 1; i <= n ~/ 2 - 1; i++) {
      ph[i + 1] = ph0[i];
    }

    ph[n ~/ 2 + 1] = 0.0;
    for (i = 1; i <= n ~/ 2 - 1; i++) {
      ph[n - i + 1] = -ph0[i];
    }

    // Make complex spectrum
    for (i = 1; i <= n; i++) {
      swC[2 * i - 1] = sw[i] * cos(ph[i]);
    }

    for (i = 1; i <= n; i++) {
      swC[2 * i] = sw[i] * sin(ph[i]);
    }

    // Calculate inverse fft
    ifft(swC, n, -1);

    // Extract real part
    for (i = 1; i <= n; i++) {
      rr[i] = (1.0 / n) * swC[2 * i - 1];
    }

    xstd = stdev(rr, n);
    ratio = rrstd / xstd;

    for (i = 1; i <= n; i++) {
      rr[i] *= ratio;
    }

    for (i = 1; i <= n; i++) {
      rr[i] += rrmean;
    }
  }

  /*
   * DETECT PEAKS
   */
  void detectpeaks(List<double> ipeak, List<double> x, List<double> y, List<double> z, int n) {
    int i, j, j1, j2, jmin, jmax, d;
    double thetap1, thetap2, thetap3, thetap4, thetap5;
    double theta1, theta2, d1, d2, zmin, zmax;

    thetap1 = ti[1];
    thetap2 = ti[2];
    thetap3 = ti[3];
    thetap4 = ti[4];
    thetap5 = ti[5];

    for (i = 1; i <= n; i++) {
      ipeak[i] = 0.0;
    }

    theta1 = atan2(y[1], x[1]);

    for (i = 1; i < n; i++) {
      theta2 = atan2(y[i + 1], x[i + 1]);

      if ((theta1 <= thetap1) && (thetap1 <= theta2)) {
        d1 = thetap1 - theta1;
        d2 = theta2 - thetap1;
        if (d1 < d2) {
          ipeak[i] = 1.0;
        } else {
          ipeak[i + 1] = 1.0;
        }
      } else if ((theta1 <= thetap2) && (thetap2 <= theta2)) {
        d1 = thetap2 - theta1;
        d2 = theta2 - thetap2;
        if (d1 < d2) {
          ipeak[i] = 2.0;
        } else {
          ipeak[i + 1] = 2.0;
        }
      } else if ((theta1 <= thetap3) && (thetap3 <= theta2)) {
        d1 = thetap3 - theta1;
        d2 = theta2 - thetap3;
        if (d1 < d2) {
          ipeak[i] = 3.0;
        } else {
          ipeak[i + 1] = 3.0;
        }
      } else if ((theta1 <= thetap4) && (thetap4 <= theta2)) {
        d1 = thetap4 - theta1;
        d2 = theta2 - thetap4;
        if (d1 < d2) {
          ipeak[i] = 4.0;
        } else {
          ipeak[i + 1] = 4.0;
        }
      } else if ((theta1 <= thetap5) && (thetap5 <= theta2)) {
        d1 = thetap5 - theta1;
        d2 = theta2 - thetap5;
        if (d1 < d2) {
          ipeak[i] = 5.0;
        } else {
          ipeak[i + 1] = 5.0;
        }
      }
      theta1 = theta2;
    }

    // Correct the peaks
    d = (paramOb.getSfEcg() / 64).ceil();
    for (i = 1; i <= n; i++) {
      if (ipeak[i] == 1 || ipeak[i] == 3 || ipeak[i] == 5) {
        j1 = max(1, i - d);
        j2 = min(n, i + d);
        jmax = j1;
        zmax = z[j1];
        for (j = j1 + 1; j <= j2; j++) {
          if (z[j] > zmax) {
            jmax = j;
            zmax = z[j];
          }
        }
        if (jmax != i) {
          ipeak[jmax] = ipeak[i];
          ipeak[i] = 0;
        }
      } else if (ipeak[i] == 2 || ipeak[i] == 4) {
        j1 = max(1, i - d);
        j2 = min(n, i + d);
        jmin = j1;
        zmin = z[j1];
        for (j = j1 + 1; j <= j2; j++) {
          if (z[j] < zmin) {
            jmin = j;
            zmin = z[j];
          }
        }
        if (jmin != i) {
          ipeak[jmin] = ipeak[i];
          ipeak[i] = 0;
        }
      }
    }
  }

  /*
   * DORUN PART OF PROGRAM
   */
  bool dorun() {
    bool retValue = true;

    int i, j, k, nrr, nt, nts;
    int q;
    List<double> x;
    double tstep, tecg, rrmean, hrfact, hrfact2;
    double qd;
    List<double> xt, yt, zt, xts, yts, zts;
    double timev, zmin, zmax, zrange;
    List<double> ipeak;

    // Perform some checks on input values
    q = (paramOb.getSf() ~/ paramOb.getSfEcg()).round();
    qd = paramOb.getSf() / paramOb.getSfEcg();

    /* Convert angles from degrees to radians and copy a vector to ai*/
    for (i = 1; i <= 5; i++) {
      ti[i] = paramOb.getTheta(i - 1) * PI / 180.0;
      ai[i] = paramOb.getA(i - 1);
    }

    /* Adjust extrema parameters for mean heart rate */
    hrfact = sqrt(paramOb.getHrMean() / 60);
    hrfact2 = sqrt(hrfact);

    for (i = 1; i <= 5; i++) {
      bi[i] = paramOb.getB(i - 1) * hrfact;
    }

    ti[1] *= hrfact2;
    ti[2] *= hrfact;
    ti[3] *= 1.0;
    ti[4] *= hrfact;
    ti[5] *= 1.0;

    /* Declare state vector */
    x = List<double>.filled(4, 0.0);

    ecgLog.println("Approximate number of heart beats: ${paramOb.getN()}");
    ecgLog.println("ECG sampling frequency: ${paramOb.getSfEcg()} Hertz");
    ecgLog.println("Internal sampling frequency: ${paramOb.getSf()} Hertz");
    ecgLog.println("Amplitude of additive uniformly distributed noise: ${paramOb.getANoise()} mV");
    ecgLog.println("Heart rate mean: ${paramOb.getHrMean()} beats per minute");
    ecgLog.println("Heart rate std: ${paramOb.getHrStd()} beats per minute");
    ecgLog.println("Low frequency: ${paramOb.getFLo()} Hertz");
    ecgLog.println("High frequency std: ${paramOb.getFHiStd()} Hertz");
    ecgLog.println("Low frequency std: ${paramOb.getFLoStd()} Hertz");
    ecgLog.println("High frequency: ${paramOb.getFHi()} Hertz");
    ecgLog.println("LF/HF ratio: ${paramOb.getLfHfRatio()}");
    ecgLog.println("time step milliseconds: ${paramOb.getEcgAnimateInterval()}\n");
    ecgLog.println("Order of Extrema:");
    ecgLog.println("      theta(radians)");
    ecgLog.println("P: [${ti[1]}\t]");
    ecgLog.println("Q: [${ti[2]}\t]");
    ecgLog.println("R: [${ti[3]}\t]");
    ecgLog.println("S: [${ti[4]}\t]");
    ecgLog.println("T: [${ti[5]}\t]\n");
    ecgLog.println("      a(calculated)");
    ecgLog.println("P: [${ai[1]}\t]");
    ecgLog.println("Q: [${ai[2]}\t]");
    ecgLog.println("R: [${ai[3]}\t]");
    ecgLog.println("S: [${ai[4]}\t]");
    ecgLog.println("T: [${ai[5]}\t]\n");
    ecgLog.println("      b(calculated)");
    ecgLog.println("P: [${bi[1]}\t]");
    ecgLog.println("Q: [${bi[2]}\t]");
    ecgLog.println("R: [${bi[3]}\t]");
    ecgLog.println("S: [${bi[4]}\t]");
    ecgLog.println("T: [${bi[5]}\t]\n");

    /* Initialize the vector */
    x[1] = xinitial;
    x[2] = yinitial;
    x[3] = zinitial;

    /* Initialize seed */
    rseed = -paramOb.getSeed();

    /* Calculate time scales */
    h = 1.0 / paramOb.getSf();
    tstep = 1.0 / paramOb.getSfEcg();

    /* Calculate length of RR time series */
    rrmean = (60.0 / paramOb.getHrMean());
    nrr = pow(2.0, (log(paramOb.getN() * rrmean * paramOb.getSf()) / log(2.0)).ceil()).toInt();

    ecgLog.println("Using $nrr = 2^ ${(log(1.0 * nrr) / log(2.0)).round()} samples for calculating RR intervals");

    /* Create rrprocess with required spectrum */
    rr = List<double>.filled(nrr + 1, 0.0);
    rrprocess(rr, paramOb.getFLo(), paramOb.getFHi(), paramOb.getFLoStd(), paramOb.getFHiStd(), paramOb.getLfHfRatio(), paramOb.getHrMean(), paramOb.getHrStd(), paramOb.getSf().toDouble(), nrr);

    /* Create piecewise constant rr */
    rrpc = List<double>.filled((2 * nrr) + 1, 0.0);
    tecg = 0.0;
    i = 1;
    j = 1;
    while (i <= nrr) {
      tecg += rr[j];
      j = (tecg / h).round();
      for (k = i; k <= j; k++) {
        rrpc[k] = rr[i];
      }
      i = j + 1;
    }
    nt = j;

    /* Integrate dynamical system using fourth order Runge-Kutta*/
    xt = List<double>.filled(nt + 1, 0.0);
    yt = List<double>.filled(nt + 1, 0.0);
    zt = List<double>.filled(nt + 1, 0.0);
    timev = 0.0;
    for (i = 1; i <= nt; i++) {
      xt[i] = x[1];
      yt[i] = x[2];
      zt[i] = x[3];
      rk4(x, mstate, timev, h, x);
      timev += h;
    }

    /* Down sample to ECG sampling frequency */
    xts = List<double>.filled(nt + 1, 0.0);
    yts = List<double>.filled(nt + 1, 0.0);
    zts = List<double>.filled(nt + 1, 0.0);

    j = 0;
    for (i = 1; i <= nt; i += q) {
      j++;
      xts[j] = xt[i];
      yts[j] = yt[i];
      zts[j] = zt[i];
    }
    nts = j;

    /* Do peak detection using angle */
    ipeak = List<double>.filled(nts + 1, 0.0);
    detectpeaks(ipeak, xts, yts, zts, nts);

    /* Scale signal to lie between -0.4 and 1.2 mV */
    zmin = zts[1];
    zmax = zts[1];
    for (i = 2; i <= nts; i++) {
      if (zts[i] < zmin) {
        zmin = zts[i];
      } else if (zts[i] > zmax) zmax = zts[i];
    }
    zrange = zmax - zmin;
    for (i = 1; i <= nts; i++) {
      zts[i] = (zts[i] - zmin) * (1.6) / zrange - 0.4;
    }

    /* Include additive uniformly distributed measurement noise */
    for (i = 1; i <= nts; i++) {
      zts[i] += paramOb.getANoise() * (2.0 * ran1() - 1.0);
    }

    /*
     * Insert into the ECG data table
     */
    ecgLog.println("Generating result matrix...");

    ecgResultNumRows = nts;

    ecgResultTime = List<double>.filled(ecgResultNumRows, 0.0);
    ecgResultVoltage = List<double>.filled(ecgResultNumRows, 0.0);
    ecgResultPeak = List<int>.filled(ecgResultNumRows, 0);

    for (i = 1; i <= nts; i++) {
      ecgResultTime[i - 1] = (i - 1) * tstep;
      ecgResultVoltage[i - 1] = zts[i];
      ecgResultPeak[i - 1] = ipeak[i].toInt();
    }

    print ('ecgResultNumRows->$ecgResultNumRows');

    File file = File(outfile);

    IOSink sink = file.openWrite();
    for (i = 0; i < nts; i++) {
      sink.write('${(i) * tstep} ${zts[i]} ${ipeak[i].toInt()}\n');
    }

    sink.close();

    ecgLog.println("Finished generating result matrix.");

    return retValue;
  }
}

