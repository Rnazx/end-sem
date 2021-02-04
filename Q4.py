import library as diff
g=9.8
def f(z,y,x):#The second order differential equation
    return -(g**2)
def j(z,x):#Splitting into two first order equations
    return z
ydash=diff.boundary(0,5,2,45,1.5,j,f,0.01,0.001)
print(ydash)