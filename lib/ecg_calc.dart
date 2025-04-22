import 'package:ecg_calc/ecg_log.dart';

import 'ecg_calc_imp.dart';
import 'ecg_param.dart';

int calculate() {
  EcgCalc calc = EcgCalc(EcgParam(), EcgLogWindow());
  bool rc = calc.calculateEcg();
  print ('rc=$rc');
  return 0;
}
