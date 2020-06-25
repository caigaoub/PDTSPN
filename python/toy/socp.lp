\ Model SOCP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  Z[0] + Z[1] + Z[2] + Z[3]
Subject To
 cw0: X[2] - W[0] = 0
 cw1: - X[0] - W[1] = 0
 cw2: X[0] - X[1] - W[2] = 0
 cw3: X[1] - X[2] - W[3] = 0
 cu0: Y[2] - U[0] = 0
 cu1: - Y[0] - U[1] = 0
 cu2: Y[0] - Y[1] - U[2] = 0
 cu3: Y[1] - Y[2] - U[3] = 0
 Dx_0: - X[0] - S[0] = 0
 Dx_1: - X[1] - S[1] = -40
 Dx_2: - X[2] - S[2] = -40
 Dy_0: - Y[0] - T[0] = -40
 Dy_1: - Y[1] - T[1] = -40
 Dy_2: - Y[2] - T[2] = 0
 Dist_0: [ - Z[0] ^2 + W[0] ^2 + U[0] ^2 ] <= 0
 Dist_1: [ - Z[1] ^2 + W[1] ^2 + U[1] ^2 ] <= 0
 Dist_2: [ - Z[2] ^2 + W[2] ^2 + U[2] ^2 ] <= 0
 Dist_3: [ - Z[3] ^2 + W[3] ^2 + U[3] ^2 ] <= 0
 InCir_0: [ S[0] ^2 + T[0] ^2 ] <= 1
 InCir_1: [ S[1] ^2 + T[1] ^2 ] <= 1
 InCir_2: [ S[2] ^2 + T[2] ^2 ] <= 1
Bounds
 X[0] free
 X[1] free
 X[2] free
 Y[0] free
 Y[1] free
 Y[2] free
 W[0] free
 W[1] free
 W[2] free
 W[3] free
 U[0] free
 U[1] free
 U[2] free
 U[3] free
 S[0] free
 S[1] free
 S[2] free
 T[0] free
 T[1] free
 T[2] free
End
