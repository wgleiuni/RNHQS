#!/usr/bin/python
import subprocess
from subprocess import Popen, PIPE
import numpy as np

h = 0.002
h = h*np.pi
sigma = 2.1
ri = 0.5
omega = 1
phi = 0.0
numdt = 20000
step = 11
omega = np.linspace(0, 10, step)
phi = np.linspace(0, 2, step)
[omegav, phiv] = np.meshgrid(omega, phi)
for i in range(step):
    for j in range(step):
        tN = (i)*step+j+1
        print tN
        cmd = ["./python_job.sh", "%d" % tN, "%f" % h, "%f" % sigma, "%f" % ri,
               "%f" % omegav[i, j], "%f" % phiv[i, j], "%d" % numdt]
        result = subprocess.Popen(cmd, stdout = subprocess.PIPE)
        outfile = open('p%d.txt' % tN, 'a+')
        out = result.stdout.read()
        outfile.write(out)
        outfile.close()
