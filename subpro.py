#!/usr/bin/python
import subprocess
from subprocess import Popen, PIPE
import numpy as np

bMode = 3
oMode = 1
h = 0.002
gama = 0.1
v = 0.0501
A = 1
omega = 10
f = 0.25
phi = 0.5
numdt = 20000000
W = 0
for tN in range(1, 3):
    cmd = ["./a.out", "%d" % tN, "%d" % bMode, "%d" % oMode, "%f" % h, "%f" %
           gama, "%f" % v, "%f" % A, "%f" % omega, "%f" % f, "%f" % phi, "%d" %
           numdt, "%f" % W]
    result = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    outfile = open('p%d.txt' % tN, 'a+')
    out = result.stdout.read()

    outfile.write(out)
    outfile.close()
