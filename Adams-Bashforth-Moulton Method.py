from sympy import *
import time
from math import *

"""
Name: Natasha Ambeskar
Date: 04/01/2024
Description: Adams-Bashforth-Moulton Method for approximating the solution to first order differential equations at a desired point.
Parameters:
Differential Equation
initial x and y values for the function that is a solution to the given differential equation
Final x value at which point the method should find the value of the solution to the differntial equation
Step size desired by user when traversing the interval from initial to final x value

"""
# Parameters to be changed by user
x = Symbol("x")
y = Symbol("y")
diffEq = x+y-1
x0 = 0
y0 = 1
step = 0.2
xf = 0.8
 
def RK4(diffEq, x0, y0, step, xf):                                         #rk4 method; finds intial points for Adam Bashforth Method
    x = Symbol("x")
    y = Symbol("y")
    iterations = int((xf - x0) / step)                                     # indicates how many attempts it will take to get from x0 to xf with a given step size
    start_time = time.time()
    InitialConditions = [y0]
    for i in range(3):
        k1 = diffEq.subs({x: x0, y: y0})                                   # k1-slope at x0,y0
        k2 = diffEq.subs({x: x0 + (step / 2) , y: y0 + step * (k1 / 2)})   # k2
        k3 = diffEq.subs({x: x0 + (step / 2) , y: y0 + step * (k2 / 2)})   # k3
        k4 = diffEq.subs({x: x0 + step, y: y0 + step * (k3)})              # k4
        y0 += (step/6) * (k1 + 2*k2 + 2*k3 + k4)                           # update y0
        InitialConditions.append(y0)
        x0 += step                                                         # increment x0 by step to go to next iteration
    return InitialConditions


def ABM(diffEQ,x0,y0,step, xf):                                            #ABM method
    r=(xf-(x0+3*step))
    iterations=int(round(r,2)/step)                                        #calculates number of iterations necessary to reach final step value from last rk4 x-value
    InitialCondition = RK4(diffEq, x0, y0, step, xf)                       #calls RK4 method to determine other three initial conditions(that will be stored in a list)
    y1 = InitialCondition[1];                                              #calls initial conditions in from list and stores in variables
    y2 = InitialCondition[2];
    y3 = InitialCondition[3];
    x1=x0+step
    x2=x1+step
    x3=round(x2+step,2)
    for i in range(iterations):
        f0 = diffEQ.subs({x:x0,y:y0})                                      #calculates slopes at initial conditions
        f1 = diffEQ.subs({x:x1,y:y1})
        f2 = diffEQ.subs({x:x2,y:y2})
        f3 =diffEQ.subs({x:x3,y:y3})
        yPredicted = y3 + (step / 24) * (55*f3-59*f2+37*f1-9*f0)           #calculates initial prediction at x3+step
        f4 = diffEQ.subs({x:x0 + 4 * step,y:yPredicted})                   #calculates slope at using predicted value
        yPrev = yPredicted
        yCurr = y3 + (step/24) * (9*f4+19*f3-5*f2+f1)                      #corrects predicted value using slope calculated
        while(yPrev.round(11)!=yCurr.round(11)):                           #reapplies corrector until values are the same before and after using corrector method
            yPrev=yCurr
            f4=diffEQ.subs({x:x0 + 4 * step,y:yCurr})
            yCurr=y3+(step/24)*(9*f4+19*f3-5*f2+f1)
        y0=y1                                                             #shifts values to create new set of 'four initial conditions' using previously calculated value for the next step.
        y1=y2
        y2=y3
        y3=yCurr
        x0=x1
        x1=x2
        x2=x3
        x3+=step
    return yCurr

print(ABM(diffEq, x0, y0, step, xf))
