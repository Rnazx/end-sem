import math
import library as lib
L=1
g=9.8
thetamax=math.pi/4
a=math.sin(thetamax/2)
def func(x):
    return 4*math.sqrt(L/g)*(1/math.sqrt(1-(a*math.sin(x))**2))
y=lib.simpson(0,(math.pi/2),10,func)
print("The value of the time period is found to be"+str(y)+"seconds")
