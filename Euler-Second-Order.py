'''
File Name: Euler's Method-Second Order
Author: Natasha Ambeskar
Date: 3/6/2024

Description: Performs Euler's method to approximate solution of second order differential equation at given point. 
Parameters: 
Second Oder Differntial Equation as a string
Inital x, y, and y prime values for solution of the differential equation
Final x value at for which the method should approximate the solution
Size of interval(step size) desired by user when traversing interval from initial to final x value(smaller step size leads to increased accuracy)
'''

from sympy import *
from math import *
init_printing(use_unicode=True)
x = Symbol("x")
y = Symbol("y")
yp = Symbol("yp") #represents y prime

#parameters to be changed by user

eqn2='-1*y+-x*yp' #Second order differential equation, can be changed by user.
x0=0
y0=1
yp0=2
xf=0.2
step=0.1




def EulersSecondOrder(eqn2,x0,y0,yp0,xf,step):#Second order Euler's method, accepts second order differential(eqn2) as a string along with x0,y0,yp0,xf,step size 
    M=Matrix([[0,1]])
    CoeffList=[]
    for i in range(2):                        #Extracts coefficients from the equation and adds to a matrix
        var=eqn2.find("y")
        appending=sympify(eqn2[:var-1])
        CoeffList.append(appending)
        eqn2=eqn2[var+2:]
    M=M.row_insert(1,Matrix([CoeffList]))

    a=Matrix([y,yp])                         #Creates a new matrix to store y and y prime values

    iterations=int((xf-x0)/step)
    for i in range(iterations):             #Performs Euler's method on Matrix equation: a(next)=a(previous)+step*(M(with previous values subbed in)*a(previous)
        f=M.subs({x:x0})*a.subs({y:y0,yp:yp0})
        a=a.subs({y:y0,yp:yp0})+step*f
        y0=a.row(0)[0]
        yp0=a.row(1)[0]
        x0+=step
    return a                                 #returns a 2x1 matrix with the final y value in the first row and the final y prime value in the second row


yfMatrix=EulersSecondOrder(eqn2,x0,y0,yp0,xf,step) #parameters to be changed by user

print('yf=',yfMatrix.row(0)[0])
print('yfprime=',yfMatrix.row(1)[0])

