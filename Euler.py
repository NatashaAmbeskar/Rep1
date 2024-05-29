"""
File Name: Euler's Method
Author: Natasha Ambeskar
Date:3/6/2024

Description:Uses Euler's Method to approximate value of solution to a given differntial equation at a desired point.
Parameters: 
Differential Equation(diffEQ)
Initial x(x0) and y(y0) values for the solution of the given differential equation
final x(xf) value for which the method should find the value of the solution of the diffEQ
step size(step) desired by user when travesring the interval from initial to final x values(small step size results in higher accuracy)

"""
from sympy import *
from sympy import symbols
x, y=symbols('x y')
#comment
#Paameters to be changed by users
diffEQ=y
x0=0
y0=1
step=0.1
xf=5

def Eulers(diffEQ,x0,y0,xf,step):             #uses Euler's method to approximate value of function at xF 
    if(step>(xf-x0)):                         #prevents steps that go past the final value
        print("step too big")
    elif(step<10**(-6)):                      #prevents step size of 0
        print("step too small")
    else:
        iterations=int((xf-x0)/step)          #calculates the number of times Euler's method must be performed to reach approximation for value of function at xF
        for i in range(iterations):           #Loop that performs the Euler's method the number of times necessary to reach xF
            m=diffEQ.subs(x,x0).subs(y,y0)
            y1=y0+m*step
            x0+=step
            y0=y1
        return y0
    
   


yf=Eulers(diffEQ,x0,y0,xf,step)

print(yf)
