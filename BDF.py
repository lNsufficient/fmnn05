from assimulo.explicit_ode import *
import numpy as N
import scipy.linalg as SL
import matplotlib.pyplot as plt
import math

class BDF(Explicit_ODE):
    """
    Explicit Euler.
    """
    Tol=1.e-8
    maxit=100
    maxsteps = 100000
    order=3
    
    def integrate(self, t0, y0, tf, opts):
        """
        _integrates (t,y) values until t > tf
        """

        self.jCalc = 1
        if opts["output_list"] == None:
            raise Explicit_ODE_Exception('BDF-2 is a fixed step-size method. Provide' \
                                         ' the number of communication points.')
        
        self.h = N.diff(opts["output_list"])[0]
       
        
        t = N.array([t0])
        y = N.array([y0]).T
        print(N.shape(y))
        for i in range(self.maxsteps):
            print('t och y ' , N.shape(t), N.shape(y))
            if t[-1] >= tf:
                break
            if i==0:  # initial step
                t,y = self.step_EE(t,y)
            elif i==1 or self.order==2:   
                t,y = self.step_BDF(t,y, 2) #tilldelningssatsen behövs förmodligen inte.
            elif i==2 or self.order==3:
                t,y = self.step_BDF(t,y,3)
            else:
                t,y = self.step_BDF(t,y,4)

            self.h=min(self.h,N.abs(tf-t[-1]))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        self.h = N.diff(opts["output_list"])[0]
        y = y.T
        return 3, list(t), list(y)
    
    def step_EE(self, T, Y):
        """
        This calculates the next step in the integration with explicit Euler.
        """
        f = self.problem.rhs
        h = self.h
        print(N.shape(Y))
        T = N.append(T,T[-1]+h)
        #append = Y[:,-1].reshape(2,-1)+h*(f(T[-2],Y[:,-1]).reshape(2,-1))
        append = Y[:,-1]+h*(f(T[-2],Y[:,-1])).T
        
        print('size to append', N.shape(append))
        Y = N.append(Y,append.T, axis=1)
        
        return T, Y
        
    def step_BDF(self,T,Y,bdf_order):
        """
        B        append = Y[:,-1].reshape(2,-1)+h*(f(T[-2],Y[:,-1]).reshape(2,-1))
DF-2 with Fixed Point Iteration and Zero order predictor
        
        alpha_0*y_np1+alpha_1*y_n+alpha_2*y_nm1=h f(t_np1,y_np1)
        alpha=[3/2,-2,1/2]
        """
        if bdf_order==2:
            alpha=[3./2.,-2.,1./2]
        elif bdf_order==3:
            alpha=[11./6,-3,3./2,-1./3]
        elif bdf_order==4:
            alpha=[25./12.,-4,3,-4./3,1./4.]
        f=self.problem.rhs
        h=self.h
        # predictor
        T = N.append(T,T[-1]+h)
        appendix = Y[:,-1]
        print('Size of appendix', N.shape(appendix))
        print('Size of vector to append to', N.shape(Y))
        
        #Y = N.append(Y,Y[:,-1].reshape(2,-1), axis=1)
        Y = N.append(Y,N.array([Y[:,-1]]).T, axis=1)
        
        # # corrector with fixed point iteration
        # for i in range(self.maxit):
            # y_np1_ip1=h*f(T[-1],Y[:,-1])			
            # for j in range(bdf_order):
                # print('y_np1_ip1 ',N.shape(y_np1_ip1))
                # #y_np1_ip1=y_np1_ip1.reshape(2,1) - (alpha[j+1]*Y[:,-j-2]).reshape(2,1)
                # new = y_np1_ip1 - N.array([N.dot(alpha[j+1],Y[:,-j-2])]).T
                # #print(y)
                # y_np1_ip1= new 
                # print(N.shape(alpha[j+1]*Y[:,-j-2]))
            # y_np1_ip1=y_np1_ip1/alpha[0]
            # if SL.norm(y_np1_ip1[:,0]-Y[:,-1]) < self.Tol:
                # Y[:,-1]=y_np1_ip1[:,0]
                # return T,Y
            # print('Shapes ',N.shape(Y[:,-1]), N.shape(y_np1_ip1)) 
            # Y[:,-1]=y_np1_ip1[:,0]
        if self.jCalc == 1:	
            J = N.eye(N.size(Y[:,-1])) + h/alpha[0]*self.fder(f,T,Y)
        else:
            J = self.J
        #print("Linneas test", N.array([Y[:,-1]]).T)
        Y[:,-1]=h*f(T[-1],Y[:,-1])[:,0] +Y[:,-1]
        for i in range(self.maxit):
            F = h*f(T[-1], Y[:,-1])[:,0]
            for j in range(bdf_order):
                F = F - (alpha[j+1]*Y[:,-j-2])
            F = F/alpha[0]
            F = F - Y[:,-1]
            if SL.norm(F) < self.Tol:
                self.J = J
                self.jCalc = 0
                return T, Y
            else:
                #print("F ", F)
                dx = SL.solve(J,F)
                #print("dx ", dx)
                Y[:,-1] = Y[:,-1] + dx
        else:
            if self.jCalc == 1:
                raise Explicit_ODE_Exception('Corrector could not converge within % iterations'%i)
            else:
                self.jCalc = 1
                return step_BDF(T[0:-1],Y[:,:-1], bdf_order)
                
    def fder(self,f,T,Y):
        ySize = N.size(Y[:,-1])
        #print(ySize)
        J = N.zeros((ySize, ySize))
        h = self.h*1e-1
        I = N.eye(ySize)
        for j in range(ySize):
            d = h*I[:,j]
            J[:, j] = (f(T[-1],Y[:,-1]+d)[:,0]-f(T[-1],Y[:,-1]-d)[:,0])/(2*h)
			
        #print(J)
        return J
	
    def print_statistics(self, verbose=NORMAL):
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        self.log_message(' Step-length          : %s '%(self.h), verbose)
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : BDF%s'%(self.order),                       verbose)
        self.log_message(' Solver type       : Fixed step\n',                      verbose)
            

#TEST EXAMPLE
def rhs(t,y):
    k = 10 
    L = k*(math.sqrt(y[0]**2+y[1]**2)-1)/math.sqrt(y[0]**2+y[1]**2)
    result = N.array([[y[2],y[3],-y[0]*L,-y[1]*L-1]]).T
    return result
def pend(t,y):
    #g=9.81    l=0.7134354980239037
    gl=13.7503671
    print(N.size(y))
    result = N.array([[y[1]],[-gl*N.sin(y[0])]])
    print('result ',N.shape(result))
    return result

#Initial conditions
y0=N.array([0.6, -0.7,0.,0.])    

#Specify an explicit problem
pend_mod=Explicit_Problem(rhs, y0)
pend_mod.name='Nonlinear Pendulum'

#Define an explicit solver
pend_sim = BDF(pend_mod) #Create a BDF solver

#Simulate the problem
c = 500
t,y = pend_sim.simulate(1*c, 100*c)

#Plot the result
print('I slutet', N.shape(y))

print(N.shape(t))
#P.plot(t,y)
#P.show()
print(y)
plt.plot(y[:,0],y[:,1])
plt.axis('equal')
#plt.plot(t,y)
plt.show()



