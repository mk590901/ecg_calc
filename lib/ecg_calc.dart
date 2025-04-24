import 'ecg_calc_imp.dart';
import 'ecg_param.dart';
import 'ecg_log.dart';

bool calculate() {
  EcgCalc calc = EcgCalc(EcgParam(), EcgLogWindow());
  bool rc = calc.calculateEcg();
  return rc;
}

