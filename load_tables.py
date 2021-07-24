import numpy as np

odir = 'tables'
d    = 5
C    = 25

try:
    A = np.fromfile(odir + '/A{}-{}.dat'.format(C,d), dtype=np.float64)
    A = A.reshape((C,d), order="F")
    print(A)

try:
    Q_dir = np.fromfile(odir + '/Q{}.dat'.format(d), dtype=np.float64)
    Q_dir = Qdir.reshape((d,d), order="F")
    print(Q_dir)


try:
    Q_neu  = np.fromfile(odir + '/Qn{}.dat'.format(d), dtype=np.float64)
    dx_neu = Qneu[0]
    Q_neu  = Qneu[1:].reshape((d,d), order="F")
    print(dx_neu)
    print(Q_neu)


try:
    Q_neu2  = np.fromfile(odir + '/Q2n{}.dat'.format(d), dtype=np.float64)
    dx_neu2 = Qneu2[0]
    Q_neu2  = Qneu2[1:].reshape((d,d), order="F")
    print(dx_neu2)
    print(Q_neu2)
