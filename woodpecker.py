import scipy
import scipy.optimize as so
import scipy.linalg as sl
import numpy
import math
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
import matplotlib.pyplot as plt

class Woodpecker(Implicit_Problem):

def variables():
    global m_b, m_s
    m_b = 
    m_s =

    
#    zpp, phipp_s, phipp_b = yp[3], yp[4], yp[5]
#    zp, phip_s, phip_b = yp[0], yp[1], yp[2]
#    lambdap_1, lambdap_2  = yp[6], yp[7]
#    z, phi_s, phi_b = y[0], y[1], y[2]
#    v, u_s, u_b = y[3], y[4], y[5]
#    lambda_1, lambda_2  = y[6], y[7]
def state_event(self,y, yp, lambda_1_old):
    phi_b=y[2]
    phi_s=y[1]
    lambda_1=y[6]


    if self.current_state == 1:
        
        if phip_b < 0:
            new_state = 2
        if phip_b > 0:
            new_state = 3

    elif self.current_state == 2 and lambda_1_old*lambda_1 < 0:
        new_state = 1

    elif self.current_state == 3:
        if phip_b < 0 and lambda_1_old*lambda_1 < 0: #maybe use tol? 
            new_state = 1
        if phip_b > 0 and h_b*phi_b == l_s+l_g-l_b-r_0:
            new_state = 4

    elif self.current_state == 4:
        new_state = 3
        phi_s = phi_s*(-1)

    return phi_s, new_state

def handle_event(self,y,yp,new_state):

    if new_state==1:
        self.res = self.motion_state1
    elif new_state==2:
        self.res = self.motion_state2
    elif new_state==3:
        self.res = self.motion_state3
    elif new_state==4:
        y[4] = -y[4]
        yp[1] = -yp[1]
        new_state=3
    self.current_state = new_state

    return y, yp
        

def motion_state1(t,y, yp):

    lambda_1, lambda_2  = y[6], y[7]
    
    r = default_residual(y, yp)
    r[0] = r[0]
    r[1] = r[1] + lambda_1
    r[2] = r[2] + lambda_2
    r[3] = lambda_1
    r[4] = lambda_2
    
    return r
def motion_state2(t,y, yp):

    zp, phi_s, phip_s = yp[0], y[1], yp[1]
    lambda_1, lambda_2  = y[6], y[7]

    
    r = default_residual(y, yp)
    r[0] = r[0] + lambda_2
    r[1] = r[1] + h_s*lambda_1+r_s*lambda_2
    r[2] = r[2]
    r[3] = (r_s - r_0) + h_s*phi_s
    r[4] = zp + r_s*phip_s

    return r

def motion_state3(t,y,yp):


    zp, phi_s, phip_s = yp[0], y[1], yp[1]
    lambda_1, lambda_2  = y[6], y[7]

    r = default_residual(y, yp)
    r[0] = r[0] + lambda_2
    r[1] = r[1] - h_s*lambda_1+r_s*lambda_2
    r[2] = r[2]
    r[3] = (r_s - r_0) - h_s*phi_s
    r[4] = zp + r_s*phip_s

    return r

def default_residual(y, yp):
    
    zpp, phipp_s, phipp_b = yp[3], yp[4], yp[5]
    zp, phip_s, phip_b = yp[0], yp[1], yp[2]
    z, phi_s, phi_b = y[0], y[1], y[2]
    v, u_s, u_b = y[3], y[4], y[5]

    r = numpy.zeros(size(y))
    r[0] = (m_s + m_b)*zpp + m_b*l_s*phipp_s+m_b*l_g*phipp_b+(m_s+m_b)*g
    r[1] = (m_b*l_s)*zpp+(j_s+m_b*l_s*l_s)*phipp_s + (m_b*l_s*l_g)*phipp_b -c_p*(phi_b-phi_s) + m_b*l_s*g
    r[2]  = m_b*l_g*zpp + (m_b*l_s*l_g)*phipp_s+ (j_b + m_b*l_g*l_g)*phipp_b - c_p*(phi_s-phi_b)+m_b*l_g*g
    r[3] = 63 
    r[4] = 63
    r[5] = v - zp
    r[6] = u_s - phip_s
    r[7] = u_b - phip_b
    
    return r
    
    
