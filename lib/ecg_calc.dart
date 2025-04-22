import 'dart:io';

import 'package:ecg_calc/ecg_log.dart';

import 'ecg_calc_imp.dart';
import 'ecg_param.dart';
import 'ecg_generator.dart';
import 'ecg_param.dart';
import 'ecg_log.dart';
import 'ecg_param2.dart';

int calculate() {
  EcgCalc calc = EcgCalc(EcgParam(), EcgLogWindow());
  bool rc = calc.calculateEcg();
  print ('rc=$rc');
  return 0;
}

void generate () async {
  // Initialize parameters and logger (replace with your actual EcgParam implementation)
  EcgParam2 param = EcgParam2()
    ..setHrMean(60.0)
    ..setHrStd(1.0)
    ..setLfHfRatio(0.5)
    ..setSfEcg(256)
    ..setSf(512)
    ..setANoise(0.0)
    ..setFLo(0.1)
    ..setFHi(0.25)
    ..setFLoStd(0.01)
    ..setFHiStd(0.01)
    ..setN(1000)
    ..setSeed(-1);
    // ..setTheta([0, -70.0, -15.0, 0.0, 15.0, 100.0])
    // ..setA([0, 1.2, -5.0, 30.0, -7.5, 0.75])
    // ..setB([0, 0.25, 0.1, 0.1, 0.1, 0.4]);

  param.setThetaValue(0, -70.0);
  param.setThetaValue(1, -15.0);
  param.setThetaValue(2, 0.0);
  param.setThetaValue(3, 15.0);
  param.setThetaValue(4, 100.0);
  param.setAValue(0, 1.2);
  param.setAValue(1, -5.0);
  param.setAValue(2, 30.0);
  param.setAValue(3, -7.5);
  param.setAValue(4, 0.75);
  param.setBValue(0, 0.25);
  param.setBValue(1, 0.1);
  param.setBValue(2, 0.1);
  param.setBValue(3, 0.1);
  param.setBValue(4, 0.4);
  EcgLogWindow log = EcgLogWindow();

  // Initialize generator
  final generator = ECGGenerator(param: param, log: log);

  // Example 1: Generate a single buffer
  var [time, voltage, peak] = generator.generateBuffer(256);
  print('Time: ${time.take(10)}...');
  print('Voltage: ${voltage.take(10)}...');
  print('Peak: ${peak.take(10)}...');

  // Example 2: Stream buffers every second
  final stream = generator.streamBuffers(256, Duration(seconds: 1));
  int count = 0;
  await for (var [time, voltage, peak] in stream) {
    print('Buffer $count: Voltage=${voltage.take(10)}..., Peak=${peak.take(10)}...');
    count++;
    if (count >= 5) break; // Stop after 5 buffers
  }

  // Optionally, write to file incrementally
  IOSink sink = File('ecgsyn.dat').openWrite();
  for (int i = 0; i < 10; i++) {
    var [time, voltage, peak] = generator.generateBuffer(256);
    for (int j = 0; j < time.length; j++) {
      //sink.write('${time[j]} ${voltage[j]} ${peak[j]}\n');
      double doubleValue = voltage[j] as double;
      int intValue = (doubleValue * 1000).toInt();
      sink.write('$intValue, ');
    }
  }
  await sink.close();
}

