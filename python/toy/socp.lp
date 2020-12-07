\ Model model_tp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  z + 0 I1 + 0 I2
Subject To
 R0: - x1 + x2 + m = 0
 R1: - y1 + y2 + n = 0
 cst_u1: - x2 - u1 + [ - x1 * a1 + x2 * a1 ] = 0
 cst_v1: - y2 - v1 + [ - y1 * a1 + y2 * a1 ] = 0
 cst_u2: - x2 - u2 + [ - x1 * a2 + x2 * a2 ] = -4
 cst_v2: - y2 - v2 + [ - y1 * a2 + y2 * a2 ] = 0
 cst_I1: [ u1 ^2 + v1 ^2 ] <= 1
 cst_I2: [ u2 ^2 + v2 ^2 ] <= 1
 qc6: [ - z ^2 + m ^2 + n ^2 ] <= 0
Bounds
 x1 free
 y1 free
 x2 free
 y2 free
 a1 <= 1
 a2 <= 1
 u1 free
 v1 free
 u2 free
 v2 free
 m free
 n free
Binaries
 I1 I2
End
