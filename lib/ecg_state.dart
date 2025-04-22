import 'dart:math';

class ECGState {
  double hrmean, hrstd, lfhfratio, sfecg, sfint, anoise;
  int n;
  List<double> ti, ai, bi;
  double t, theta, x, y, z;
  List<double> rrBuffer;
  int rrIndex, rrCount;
  int idum;
  Random rng;

  ECGState({
    required this.hrmean,
    required this.hrstd,
    required this.lfhfratio,
    required this.sfecg,
    required this.sfint,
    required this.anoise,
    required this.n,
    required this.ti,
    required this.ai,
    required this.bi,
    required int seed,
  })  : t = 0.0,
        theta = 0.0,
        x = 0.0,
        y = 0.0,
        z = 0.0,
        rrBuffer = List.filled(10, 0.0),
        rrIndex = 0,
        rrCount = 0,
        idum = seed < 0 ? seed : -seed,
        rng = Random(seed.abs()) {
    _generateInitialRR();
  }

  void _generateInitialRR() {
    // Упрощённая генерация RR (заменить на FFT-логику)
    double rrmean = 60.0 / hrmean;
    for (int i = 0; i < 10; i++) {
      rrBuffer[i] = rrmean + hrstd * (rng.nextDouble() - 0.5);
      if (rrBuffer[i] < 0.2) rrBuffer[i] = 0.2;
    }
    rrCount = 10;
  }

  void generateRR() {
    if (rrIndex >= rrCount) {
      _generateInitialRR(); // Заменить на полноценный rrprocess с FFT
      rrIndex = 0;
    }
  }

  List<double> generateECGBuffer(int bufSize) {
    List<double> buffer = List.filled(bufSize, 0.0);
    double dt = 1.0 / sfecg;

    for (int i = 0; i < bufSize; i++) {
      if (t >= t + rrBuffer[rrIndex]) {
        generateRR();
        t += rrBuffer[rrIndex];
      }

      theta += 2.0 * pi * dt / rrBuffer[rrIndex];
      if (theta > 2.0 * pi) theta -= 2.0 * pi;

      double dx = 0.0, dy = 0.0, dz = 0.0;
      for (int j = 0; j < 5; j++) {
        double dtheta = (theta - ti[j]) % (2.0 * pi);
        if (dtheta > pi) dtheta -= 2.0 * pi;
        else if (dtheta < -pi) dtheta += 2.0 * pi;
        double r = ai[j] * exp(-dtheta * dtheta / (2.0 * bi[j] * bi[j]));
        dx -= r * sin(dtheta);
        dy += r * cos(dtheta);
        dz += r;
      }

      x += dt * (-sfecg * x + dx);
      y += dt * (-sfecg * y + dy);
      z += dt * (-sfecg * z + dz);

      buffer[i] = z + anoise * rng.nextDouble();
      t += dt;
    }

    return buffer;
  }
}

void main() async {
  var state = ECGState(
    hrmean: 60.0,
    hrstd: 1.0,
    lfhfratio: 0.5,
    sfecg: 256.0,
    sfint: 512.0,
    anoise: 0.0,
    n: 1000,
    ti: [-70.0, -15.0, 0.0, 15.0, 100.0],
    ai: [1.2, -5.0, 30.0, -7.5, 0.75],
    bi: [0.25, 0.1, 0.1, 0.1, 0.4],
    seed: -1,
  );

  for (int i = 0; i < 10; i++) {
    var buffer = state.generateECGBuffer(256);
    print(buffer);
    await Future.delayed(Duration(seconds: 1));
  }
}