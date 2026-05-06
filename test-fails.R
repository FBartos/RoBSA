Failure (test-5-methods.R:33:3): Summary functions work
Expected `capture_output_lines(...)` to equal `c(...)`.
Differences:
  8/35 mismatches
x[18]: "lnorm-aft      2/8       0.167       0.000            10.031"
y[18]: "lnorm-aft      2/8       0.167       0.000            10.029"

x[23]: "x_fac3    4/8       0.500       0.189             1.454"
y[23]: "x_fac3    4/8       0.500       0.189             1.453"

x[27]: "x_cont    0.055  0.055 -0.035 0.055 0.143"
y[27]: "x_cont    0.056  0.056 -0.033 0.056 0.146"

x[28]: "x_fac3[B] 0.053  0.016 -0.076 0.016 0.282"
y[28]: "x_fac3[B] 0.053  0.020 -0.076 0.020 0.275"

x[29]: "x_fac3[C] 0.094  0.050 -0.046 0.050 0.378"
y[29]: "x_fac3[C] 0.093  0.052 -0.051 0.052 0.376"

Failure (test-5-methods.R:75:3): Summary functions work
Expected `capture_output_lines(...)` to equal `c(...)`.
Differences:
  3/20 mismatches
x[8]: "exp-aft        1/5       0.200       0.852       23.103"
y[8]: "exp-aft        1/5       0.200       0.852       23.061"

x[9]: "weibull-aft    1/5       0.200       0.062        0.265"
y[9]: "weibull-aft    1/5       0.200       0.062        0.266"

x[12]: "gamma-aft      1/5       0.200       0.084        0.369"
y[12]: "gamma-aft      1/5       0.200       0.085        0.370"

Failure (test-5-methods.R:101:3): Summary functions work
Expected `capture_output_lines(...)` to equal `c(...)`.
Differences:
  14/35 mismatches
x[9]: "exp-aft        1/5       0.200       0.842       21.294"
y[9]: "exp-aft        1/5       0.200       0.842       21.310"

x[10]: "weibull-aft    1/5       0.200       0.067        0.287"
y[10]: "weibull-aft    1/5       0.200       0.067        0.286"

x[17]: "x_cont     0.059  0.059 -0.078 0.196     1.061       1.061      0.925      1.216"
y[17]: "x_cont     0.058  0.059 -0.076 0.194     1.060       1.061      0.926      1.214"

x[18]: "x_bin[1]  -0.173 -0.175 -0.447 0.097     0.841       0.840      0.640      1.102"
y[18]: "x_bin[1]  -0.171 -0.173 -0.446 0.098     0.843       0.841      0.640      1.103"

x[19]: "x_fac3[B]  0.099  0.099 -0.242 0.433     1.104       1.104      0.785      1.542"
y[19]: "x_fac3[B]  0.099  0.099 -0.237 0.431     1.104       1.104      0.789      1.539"

Failure (test-5-methods.R:153:3): Models summary functions work
Expected `capture_output_lines(...)` to equal `c(...)`.
Differences:
  3/23 mismatches
x[16]: "     1    gamma-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.250      -320.81       0.739        8.514"
y[16]: "     1    gamma-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.250      -320.81       0.739        8.511"

x[19]: "     4   llogis-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.083      -326.19       0.001        0.012"
y[19]: "     4   llogis-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.083      -326.19       0.001        0.013"

x[21]: "     6  weibull-aft      N(0, 1)      treatment contrast: B(3, 3)       0.083      -323.69       0.014        0.154"
y[21]: "     6  weibull-aft      N(0, 1)      treatment contrast: B(3, 3)       0.083      -323.70       0.014        0.154"

Failure (test-5-methods.R:181:3): Models summary functions work
Expected `capture_output_lines(...)` to equal `c(...)`.
Differences:
  5/19 mismatches
x[10]: "     1      exp-aft                                              0.100      -317.64       0.601       13.570"
y[10]: "     1      exp-aft                                              0.100      -317.64       0.601       13.580"

x[11]: "     2  weibull-aft                                              0.100      -320.25       0.044        0.412"
y[11]: "     2  weibull-aft                                              0.100      -320.26       0.044        0.412"

x[14]: "     5    gamma-aft                                              0.100      -319.95       0.060        0.570"
y[14]: "     5    gamma-aft                                              0.100      -319.95       0.060        0.571"

x[15]: "     6      exp-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -318.51       0.250        3.001"
y[15]: "     6      exp-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -318.51       0.250        2.996"

x[19]: "    10    gamma-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -320.82       0.025        0.231"
y[19]: "    10    gamma-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -320.81       0.025        0.232"


Failure (test-8-predict.R:66:3): Predict survival works
Expected `as.matrix(prediction_models[[3]])` to equal `as.matrix(...)`.
Differences:
  14/15 mismatches (average diff: 0.114)
[1]   8.97 -  8.91 ==  0.0675
[2]   9.62 -  9.57 ==  0.0543
[3]  12.87 - 12.76 ==  0.1018
[4]   3.01 -  2.98 ==  0.0313
[5]   3.26 -  3.21 ==  0.0524
[6]   4.51 -  4.47 ==  0.0479
[8]   5.07 -  5.13 == -0.0561
[9]   6.70 -  6.68 ==  0.0185
[10]  8.40 -  8.35 ==  0.0476
