"""
File Name: Euler's Method
Author: Natasha Ambeskar
Date:3/6/2024
File Name: Euler's Method
Author: Natasha Ambeskar
Date:3/6/2024

Description:Calculates Euler's Method given first order differential equation, initial values, increment size, and final x-value.
"""
from sympy import *
from sympy import symbols
x, y=symbols('x y')

#Parameters to be changed by users
diffEQ=2*x+y
x0=0
y0=1
step=0.1
xf=5

def ImprovedEulers(diffEQ,x0,y0,xf,step):                  #uses Improved Euler's method to approximate value of function at xf
    if(step>10**-6 and step<(xf-x0)):                      #prevents steps that are too small or to big(past the final value)
        iterations=int((xf-x0)/step)                       #calculates the number of times Improved Euler's method must be performed to reach xf
        for i in range(iterations):                        #Loop that performs the Improved Euler's method the number of times necessary to reach xF
            m=diffEQ.subs(x,x0).subs(y,y0)
            y0+=(step/2)*(m+diffEQ.subs(x,x0+step).subs(y,y0+m*step))
            x0+=step
        return y0 
    else:
        print('invalid step size')
    
   


yf=ImprovedEulers(diffEQ,x0,y0,xf,step)

print(yf)