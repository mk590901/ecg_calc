import 'ecg_calc_imp.dart';
import 'ecg_param.dart';
import 'ecg_log.dart';

int calculate() {
  EcgCalc calc = EcgCalc(EcgParam(), EcgLogWindow());
  bool rc = calc.calculateEcg();
  print ('rc=$rc');
  return 0;
}

