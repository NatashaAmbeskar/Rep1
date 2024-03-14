"""
File Name: Euler's Method: error convergence
Author: Natasha Ambeskar
Date:3/6/2024

Description:Calculates Euler's Method given first order differential equation, initial values, increment size, and final x-value.
"""
from sympy import *
from sympy import symbols

x, y=symbols('x y')

#Parameters to be changed by users
diffEQ=2*x*y
x0=1.0
y0=1.0
step=0.1
xf=1.5
sol=exp((x**2)-1)

def Eulers(diffEQ,x0,y0,xf,step,sol):             #uses Euler's method to approximate value of function at xF 
    if(step>(xf-x0)):                         #prevents steps that go past the final value
        print("step too big")
    elif(step<10**(-6)):                      #prevents step size of 0
        print("step too small")
    else:
        print('for step size:', step)
        iterations=int((xf-x0)/step)
        print(f'{x0:.3f},{y0:.3f},{1:.3f},{0:.3f}')
        for i in range(iterations):
            m=diffEQ.subs(x,x0).subs(y,y0)
            y1=y0+m*step
            x0+=step
            y0=y1
            actVal=float(sol.subs(x,x0).evalf())
            error=float(abs(actVal-y0))
            print(f'{x0:.3f},{y0:.3f},{actVal:.3f},{error:.3f}') #prints x and y value at each step, the actual value of the function, and the error(actual value minus calculated value)
        #return y0
yf=Eulers(diffEQ,x0,y0,xf,step,sol)
print('')
yf=Eulers(diffEQ,x0,y0,xf,0.05,sol)

#The final error value for a step size of 0.05 is half of the value of the final error value listed for a step size of 0.05. 
#Given that the global error should be proportional to step size, deceasing the step size by a factor of one-half should also decrease the error by a factor of one-half,
#which is demonstrated in the output, thus proving that error converges. 

