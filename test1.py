import numpy
import math
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
import matplotlib.pyplot as plt

def rhs(t,y):
    k = 70 
    L = k*(math.sqrt(y[0]**2+y[1]**2)-1)/math.sqrt(y[0]**2+y[1]**2)
    result = numpy.array([y[2],y[3],-y[0]*L,-y[1]*L-1])
    return result
t0 = 0
y0 = numpy.array([0.8,0.8,0,0])

model = Explicit_Problem(rhs,y0,t0)
model.name = 'task1'

sim = CVode(model)

tfinal = 70

(t,y) = sim.simulate(tfinal)

sim.plot()

plt.plot(y[:,0],y[:,1])
plt.show()
