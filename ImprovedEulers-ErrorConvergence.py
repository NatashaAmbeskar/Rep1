"""
File Name: Improved Euler's Method: error convergence
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
diffEQ=2*x*y
x0=1
y0=1
step=0.1
xf=1.5
sol=exp((x**2)-1)
def ImprovedEulers(diffEQ,x0,y0,xf,step,sol):                  #uses Improved Euler's method to approximate value of function at xf
    if(step>10**-6 and step<(xf-x0)):                      #prevents steps that are too small or to big(past the final value)
        iterations=int((xf-x0)/step) 
        print('for step',step)
        
        for i in range(iterations):                        #Loop that performs the Improved Euler's method the number of times necessary to reach xF
            m=diffEQ.subs(x,x0).subs(y,y0)
            y0+=(step/2)*(m+diffEQ.subs(x,x0+step).subs(y,y0+m*step)) #Recalculates y value using Impoved Euler's method by taking the aveage of the left and right end of the slopes
            x0+=step
            actVal=float(sol.subs(x,x0).evalf())
            error=float(abs(actVal-y0))
            print(f'{x0:.3f},{y0:.3f},{actVal:.3f},{error:.3f}')
        #return y0
    else:
        print('invalid step size')
    
   


yf=ImprovedEulers(diffEQ,x0,y0,xf,step,sol)
print('')
yf=ImprovedEulers(diffEQ,x0,y0,xf,0.05,sol)
#The final error value for a step size of 0.05 is one-fourth of the value of the final error value listed for a step size of 0.05. Given that the global error should be proportional to step size squared for improved euler's method, deceasing the step size by a factor of one-half should decrease the error by a factor of one-fourth, which is demonstrated in the output, thus proving that error converges.
