import math
import library as n
#defining the constants
h=6.626*10**(-34)
k=1.381*10**(-23)
c=3*10**(8)
#Defining the function
def func(x):
    return (x-5)*math.exp(x)+5
print("Newton-Raphson method : ")
y=n.newtonraphson(1, func, 200,"Q1",10**(-4))
print("The root of the function is :",y )
print('\n')
print("The estimated value of b is ",(h*c)/(k*y) )

