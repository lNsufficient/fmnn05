import scipy
import scipy.optimize as so
import scipy.linalg as sl
import numpy
import math
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
import matplotlib.pyplot as plt

class Woodpecker(Implicit_Problem):
    def __init__(self, t, y0, yp0, sw):
        super().__init__(self.motion_state1, y0, yp0, t, sw0=sw)
        self.variables(y0)
        self.lambda_1_old = y0[6] 

            
    def variables(self, y):
        self.m_s = 3.0e-4 # Mass of sleeve [kg]
        self.j_s = 5.0e-9 # Moment of inertia of the sleeve [kgm]
        self.m_b = 4.5e-3 # Mass of bird [kg]
        #masstotal=mS+mB # total mass
        self.j_b = 7.0e-7 # Moment of inertia of bird [kgm]
        self.r_0 = 2.5e-3 # Radius of the bar [m]
        self.r_s = 3.1e-3 # Inner Radius of sleeve [m]
        self.h_s = 5.8e-3 # 1/2 height of sleeve [m]
        self.l_s = 1.0e-2 # verical distance sleeve origin to spring origin [m]
        self.l_g = 1.5e-2 # vertical distance spring origin to bird origin [m]
        self.h_b = 2.0e-2 # y coordinate beak (in bird coordinate system) [m]
        self.l_b = 2.01e-2 # -x coordinate beak (in bird coordinate system) [m]
        self.c_p = 5.6e-3 # rotational spring constant [N/rad]
        self.g  = 9.81 #  [m/s^2] 
#    zpp, phipp_s, phipp_b = yp[3], yp[4], yp[5]
#    zp, phip_s, phip_b = yp[0], yp[1], yp[2]
#    lambdap_1, lambdap_2  = yp[6], yp[7]
#    z, phi_s, phi_b = y[0], y[1], y[2]
#    v, u_s, u_b = y[3], y[4], y[5]
#    lambda_1, lambda_2  = y[6], y[7]

    def state_event(self,t, y, yp, sw):
        print('In state_event') 
        phi_b=y[2]
        phi_s=y[1]
        lambda_1=y[6]
        tol = 0.00001 
        lambda_1_old = self.lambda_1_old
        
        if self.current_state == 1: 
            
            if phip_b < 0 and self.h_s*phi_s < -(self.r_s - self.r_0):
                new_state = 2
            if phip_b > 0 and self.h_s*phi_s > (self.r_s - self.r_0):
                new_state = 3

        elif self.current_state == 2 and lambda_1_old*lambda_1 < 0:
            new_state = 1

        elif self.current_state == 3:
            if phip_b < 0 and lambda_1_old*lambda_1 < tol: #maybe use tol? 

                new_state = 1
            if phip_b > 0 and h_b*phi_b == l_s+l_g-l_b-r_0:
                new_state = 4

        elif self.current_state == 4:
            new_state = 3
          
        self.lambda_1_old = lambda_1 #vi kräver att detta körs varje varv, annars får vi problem
        
        new_sw = numpy.zeros(4)        

        new_sw[new_state-1] = 1
        return [list(new_sw), 0]

    def moi(self,solver):
        y = solver.y
        yp = solver.yp
        a = (self.j_b+self.m_b*self.l_g*self.l_g)
        I = self.m_b*self.l_g*yp[0] + (self.m_b*self.l_s*self.l_g)*yp[1] + a*yp[2]
        phip_b = I/a;
        yp[2] = phip_b
        y[5] = phip_b
        yp[0] = 0
        yp[1] = 0
        y[3] = 0
        y[4] = 0
        solver.y = y
        solver.yp = yp
        return solver

    def handle_event(self,solver,event_info):
        print('In handle_event')
        print(event_info) #vi misstänker att detta är vad som returneras från state_event


        print(solver)
        sw = event_info(0)

        if sw[0]==1:
            self.res = self.motion_state1
        elif sw[1]==1:
            self.res = self.motion_state2
            solver = self.moi(solver)
        elif sw[2]==1:
            self.res = self.motion_state3
            solver = self.moi(solver)
        elif sw[3]==1:
            #solver.y[4] = -solver.y[4] these are the lines that Claus suggested, we changed them. 
            solver.y[5] = -solver.y[5]
            #solver.yp[1] = -solver.yp[1]
            solver.yp[2] = -solver.yp[2]
            print("Changed values inside solver (state 4)")
            #new_state=3, but that is not necessary to change here, already set. . . . . . 

        self.current_state = new_state

        
            

    def motion_state1(self,t,y, yp, sw):
        #print('In state 1 res')
        self.state_event(t,y,yp,sw)
        lambda_1, lambda_2  = y[6], y[7]
        
        r = self.default_residual(y, yp)
        r[0] = r[0]
        r[1] = r[1] + lambda_1
        r[2] = r[2] + lambda_2
        r[3] = lambda_1
        r[4] = lambda_2
        
        return r
    def motion_state2(self, t,y, yp, sw):
        print('In state 2 res')
        zp, phi_s, phip_s = yp[0], y[1], yp[1]
        lambda_1, lambda_2  = y[6], y[7]

        
        r = self.default_residual(y, yp)
        r[0] = r[0] + lambda_2
        r[1] = r[1] + self.h_s*lambda_1+self.r_s*lambda_2
        r[2] = r[2]
        r[3] = (self.r_s - self.r_0) + self.h_s*phi_s
        r[4] = zp + self.r_s*phip_s

        return r

    def motion_state3(self,t,y,yp, sw):
        print('In state 3 res')

        zp, phi_s, phip_s = yp[0], y[1], yp[1]
        lambda_1, lambda_2  = y[6], y[7]

        r = self.default_residual(y, yp)
        r[0] = r[0] + lambda_2
        r[1] = r[1] - self.h_s*lambda_1+self.r_s*lambda_2
        r[2] = r[2]
        r[3] = (self.r_s - self.r_0) - self.h_s*phi_s
        r[4] = zp + self.r_s*phip_s

        return r

    def default_residual(self,y, yp):
        zpp, phipp_s, phipp_b = yp[3], yp[4], yp[5]
        zp, phip_s, phip_b = yp[0], yp[1], yp[2]
        z, phi_s, phi_b = y[0], y[1], y[2]
        v, u_s, u_b = y[3], y[4], y[5]

        r = numpy.zeros(numpy.size(y))
        r[0] = (self.m_s + self.m_b)*zpp + self.m_b*self.l_s*phipp_s+self.m_b*self.l_g*phipp_b+(self.m_s+self.m_b)*self.g
        r[1] = (self.m_b*self.l_s)*zpp+(self.j_s+self.m_b*self.l_s*self.l_s)*phipp_s + (self.m_b*self.l_s*self.l_g)*phipp_b -self.c_p*(phi_b-phi_s) + self.m_b*self.l_s*self.g
        r[2]  = self.m_b*self.l_g*zpp + (self.m_b*self.l_s*self.l_g)*phipp_s+ (self.j_b + self.m_b*self.l_g*self.l_g)*phipp_b - self.c_p*(phi_s-phi_b)+self.m_b*self.l_g*self.g
        r[3] = 63 
        r[4] = 63
        r[5] = v - zp
        r[6] = u_s - phip_s
        r[7] = u_b - phip_b
        
        return r
        
t0 = 0;
sw0 = [[1, 0, 0 ,0],0]
y0 = numpy.array([5, 0,0, 0, 0, 1,0,0])
yd0 =  numpy.array([0, 0, 1,-9.82, 1e-12, 0, 0, 0])
problem = Woodpecker(t0, y0, yd0, sw0)
sim = IDA(problem)
problem.name = 'Woodpecker'
ncp = 500
tfinal = 2
t, y, yd = sim.simulate(tfinal, ncp)
plt.plot(t, y)
plt.legend(["z", "phi_s", "phi_b", "zp", "phip_s", "phip_b", "lambda_1", "lambda_2"], loc = 'lower left')
plt.show()

