import math
import random
#function for storing the augmented matrix matrix(A) from a file
def store(st,rows):
    X=[]
    T=[]
    rows=int(rows)
    with open(st,'r+') as file:
        i=0
        while i<rows:
            s=file.readline()
            p=s.split()
            j=0
            while (j<len(p)):
                T.append(float(p[j]))
                j=j+1
            X.append(T)
            T=[]
            i=i+1 
        return X
#function for partial pivoting
def partialpivot(A,n,b):
    count=0
    for sum in range(n-1):
        if A[sum][sum]==0:
            for l in range(sum+1,n):
                if abs(A[l][sum])>abs(A[sum][sum]):
                    A[sum],A[l]=A[l],A[sum]
                    b[sum],b[l]=b[l],b[sum]
                    count+=1
                else : continue
        else : continue
    return count
#function for Gauss jordan method
def gaussjordan(Aug):
    n=len(Aug)
    m=len(Aug[1])
    sum=0
    flag=True
    while (sum<n) and flag:
        partialpivot(Aug,n,[0 for i in range(n)])
        piv = Aug[sum][sum]
        for l in range(sum,m):
            Aug[sum][l]/=piv
        for i in range(n):
            if i==sum or Aug[i][sum]==0: continue
            else: 
                factor=Aug[i][sum]
                for j in range(sum,m):
                    Aug[i][j]-=factor*Aug[sum][j]
        sum+=1
        #to checsum if unique solution exists
        l=0
        for i in range(n):
            if Aug[n-1][i]==0:
                l+=1
        if (l==n):
            flag=False
            if (Aug[n-1][n]!=0): print("The system has no solution")
            else: 
                print("The system has no unique solution")
                print("One of the infinite solutions is")
#function for matrix multiplication
def matrixmultiply(A,B):
    n=len(A)
    prod=[[] for x in range(n)]
    t=0
    for i in range(n):
        for j in range(n):
            for sum in range(n):
                t+=A[i][sum]*B[sum][j]
            prod[i].append(t)
            t=0
        j=0
    return prod
#function for identity matrix of order n 
def identity(n):
    I=[[int(x==y) for x in range(n)] for y in range(n)]
    return I
#function for augmenting a matrix upon the other
def augment(A,B):
    n=len(A)
    for i in range(n):
        A[i]+=B[i]
#function to print matrix
def printmatrix(A):
    n=len(A)
    m=len(A[0])
    for i in range(n):
        print()
        for j in range(m):
            print(A[i][j],end="   ")
    print('\n\n')
#function for Crout's L, U decomposition
def croutLU(A,b):
    n=len(A)
    c=0
    augment(A,A)#this will be a 3x6 matrix with LHS as L and RHS as U
    for j in range(n):
        c+=partialpivot(A,n,b)
        for i in range(n):
            l=sum=0
            while sum<i and i<=j:
                A[i][n+j] -= A[i][sum]*A[sum][n+j]
                sum+=1 
            while l<j and i>j:
                A[i][j] -= A[i][l]*A[l][n+j]
                l+=1
            if (i>j): 
                A[i][j]/=A[j][n+j]
                A[i][j+n]=0.0 #lower triangular entries in U to be 0
            else : A[i][j]=float(i==j)  #setting upper triangular entries in U to be 0 and diagnol entriesto be 1
    L=[[A[i][j]for j in range(n)] for i in range(n)]
    U=[[A[i][j+n]for j in range(n)] for i in range(n)]
    return L,U,c
#Function for forward bacsumward substitution Ax=b
def forwardbacsumward(A,b):
    n=len(A)
    sum=0
    #to checsum if lower or upper triangular
    for l in range(n):
        for m in range(n):
            if (l<m):
                sum+=A[l][m]
    if(sum==0):
        p=1 #it is lower triangular
    else: 
        p=-1 #it is upper triangular
        n+=1
    s=int(sum!=0)
    #substitution
    for i in range(s,n):
        for j in range(s,i):
            b[p*i]-=A[p*i][p*j]*b[p*j]
        b[p*i]/=A[p*i][p*i]
#function for finding the inverse of a decomposed matrix
def inverseLU(L,U,P,c):
    n=len(L)
    det=1
    for i in range(n):
        det*=U[i][i]
    det*=(-1)**c
    if (det!=0):
        X=fbm(U,fbm(L,P))
    return X,det
#forward bacsumward for matrices LUxAinv=P
def fbm(A,P):
    n=len(A)
    sum=0
    y=list(P)
    #to checsum if lower or upper triangular
    for l in range(n):
        for m in range(n):
            if (l<m):
                sum+=A[l][m]
    if(sum==0):
        p=1 #it is lower triangular
    else: 
        p=-1 #it is upper triangular
        n+=1
    s=int(sum!=0)
    #substitution
    for i in range(s,n):
        for sum in range(s,n):
            for j in range(s,i):
                y[p*i][p*sum]-=A[p*i][p*j]*P[p*j][p*sum]
            y[p*i][p*sum]/=A[p*i][p*i]
    return y
#ASSIGNMENT_05 CODE
#function for bracsumeting the root
def bracsumet(a,b,f,N=200,factor=1.5):
        n=0
        a=a
        b=b
        while (f(a)*f(b) > 0) and (n<N):
            if abs(f(a))>abs(f(b)):
                b+=factor*(b-a)
            else:
                a-=factor*(b-a)     
        if(n==199): 
            print("The root cannot be bracsumeted")
        else:
            return a,b
#function for bisection with error printing
def bisection(a, b, f, N,question,tolerance=10**(-6)):
    if (f(a)*f(b) > 0): (a,b)=bracsumet(a,b,f)#bracsumeting the root
    n = 1
    file = open("bis"+ str(question)+".txt", "w+")
    avg=a
    avgnext=b
    #print("Iteration     Error")
    while (abs(b-a) >= tolerance) and (abs(avgnext-avg) >= tolerance) and (n<N):
        avg = (a + b)/2 
        if f(a)*f(avg) < 0: b = avg
        elif f(b)*f(avg) < 0: a = avg
        elif f(avg) == 0:
            print("Found exact solution.")
            return avg
        else:
            print("Bisection method fails.")
            return None
        avgnext = (a + b)/2
        #print(n,"           ",abs(avgnext-avg))
        file.writelines([str(n)+"   ",str(abs(avgnext-avg))+"\n"]) 
        n+=1
    file.close()
    return avg

#function for regula falsi method
def regulafalsi(a, b, f, N,question,tolerance=10**(-6)):
    if (f(a)*f(b) > 0): (a,b)=bracsumet(a,b,f)#bracsumeting the root
    n = 1
    file = open("reg"+ str(question)+".txt", "w+")
    approx=a
    approxnext=b
    #print("Iteration     Error")
    while (abs(b-a) >= tolerance) and (abs(approxnext-approx) >= tolerance) and (n<N):
        approx = (a*f(b) - b*f(a))/(f(b) - f(a))
        if f(a)*f(approx) < 0: b = approx
        elif f(b)*f(approx) < 0: a = approx
        elif f(approx) == 0:
            print("Found exact solution.")
            return approx
        else:
            print("Bisection method fails.")
            return None
        approxnext = (a*f(b) - b*f(a))/(f(b) - f(a))
        #print(n,"           ",abs(approxnext-approx))
        file.writelines([str(n)+"   ",str(abs(approxnext-approx))+"\n"]) 
        n+=1
    file.close()
    return approx

#function for symmetric derivative of a function
def derivative(f, x):
    h = 10**(-6)
    return (f(x+h)-f(x-h))/(2*h)

#function for newton-raphson method
def newtonraphson(x,f,N,question,tolerance=10**(-6)):
    file = open("newt"+ str(question)+".txt", "w+")
    n = 1
    y = derivative(f, x)
    h = f(x)/y
    #print("Iteration  Error")
    while (abs(h) >= tolerance) and n<N:
        #print(n,"        ", abs(h))
        file.writelines([str(n)+"   ",str(abs(h))+"\n"]) 
        x -= h
        h = f(x)/derivative(f, x)
        n+=1
    file.close()
    return(x)
#For Laquerres method to find roots of a polynomial
#derivative of a polynomial
def polynomialderivative(x,order,coeffs):
    pol=0
    deg=len(coeffs)-1
    for i in coeffs:
        j=deg
        f=1
        sum=0
        while(sum<order):
            f*=j
            j-=1
            sum+=1
        pol+=i*f*(x**(deg-order))
        deg-=1
    return pol
#Laguerres method of estimating the root
def laguerres(a,coeffs,N=200):
    n = len(coeffs)
    x=a
    sum=0
    if(n==2):return -coeffs[1]/coeffs[0]
    else:
        if polynomialderivative(a,0,coeffs)==0:
            return a
        else:
            while(abs(a) > 10**(-6)) and sum<N:
                G = polynomialderivative(x,1,coeffs)/polynomialderivative(x,0,coeffs)
                H = G**2 - polynomialderivative(x,2,coeffs)/polynomialderivative(x,0,coeffs)
                if (abs(G + math.sqrt(abs((n-1)*(n*H - G*G))))>=abs(G - math.sqrt(abs((n-1)*(n*H - G*G))))):p=1
                else:p=-1
                a= n/(G +p*math.sqrt(abs((n-1)*(n*H - G*G))))
                x-= a
                sum+=1
            return x
#function for Synthetic division  
def syntheticdivision(est,coeffs):
    n=len(coeffs)
    B=[0 for i in range(n)]
    B[0]=coeffs[0]
    for i in range(n-1):
        B[i+1]+=coeffs[i+1]+est*(B[i])
    if(abs(B[-1])<10**(-6)): B.pop()
    else: est=laguerres(est,coeffs)
    return B
#function to find the roots of the polynomial
def roots(initialguess,coeffs):
    roots=[]
    a=initialguess
    n=len(coeffs)-1
    while (n>1):
        n=len(coeffs)-1
        a=laguerres(a,coeffs)
        roots.append(a)
        coeffs=syntheticdivision(a,coeffs)
    return roots
#*****************************************************************************************************
#INTEGRATION
# function for integrating using midpoint method
def midpoint(a,b,N,f):
    h=(b-a)/N
    y=0
    for i in range(1,N+1):
        x=a+h
        m=(x+a)/2
        y+=h*f(m)                     
        a=x
    return y
# integration using Trapezoid method 
def trapezoid(a,b,N,f):
    h=(b-a)/N
    y=h/2*(f(a)+f(b))
    for i in range(1,N):
        x=a+h
        y+=(h/2)*2*f(x)          
        a=x
    return y
# integration using Simpson method 
def simpson(a,b,N,f):
    h=abs(b-a)/N
    y=h/3*(f(a)+f(b))
    for i in range(1,N):
        x=a+h
        if i%2==0:                       
            y+=h/3*f(x)*2
        else:                            
            y+=h/3*f(x)*4
        a=x
    return y
# integration using Montecarlo method 
def montecarlo(a,b,N,f):
    sum=0
    for i in range(1,N+1):
        x = random.uniform(a,b)                
        sum+=sum+f(x)
    return sum/N
#**************************************************************************************************************
#DIFFERENTIAL EQUATIONS
#function for euler method
def euler(x,y,f,h,xend):
    file = open("Euler"+ str(h)+".txt", "w+")
    while (x<=xend):
        file.writelines([str(x)+"   ",str(y)+"\n"])
        y+=h*f(y,x)
        x+=h
    file.close()
'''#function for  RK4 method first order
def rungekutta41(x,y,f,h,xend):
    file = open("RK4"+ str(h)+".txt", "w+")
    while (x<=xend):
        file.writelines([str(x)+"   ",str(y)+"\n"])
        k1=h*f(y,x)
        k2=h*f((y+(k1/2)),(x+(h/2)))
        k3=h*f((y+(k2/2)),(x+(h/2)))
        k4=h*f((y+k3),(x+h))
        y+=(k1+2*k2+2*k3+k4)/6
        x+=h
    file.close()'''
#function for  RK4 method second order
def rungekutta4(x,y,z,dydx,dzdx,h,xstart,xend,flag=True,name="RK4"):
    if (flag==True):  file = open(str(name)+ str(h)+".txt", "w+")
    a=x
    b=y
    c=z
    g=-h
    #here from the initial point, x goes forward (in this case from 0 to 5) whereas a goes backwards (in this case from 0 to -5)
    while (x<=xend) or (a>=xstart):
        if (flag==True): file.writelines([str(x)+"   ",str(y)+"\n"])

        k1y=h*dydx(z,x)
        k1z=h*dzdx(z,y,x)
        k1b=g*dydx(c,a)
        k1c=g*dzdx(c,b,a)

        k2y=h*dydx(z+k1z/2,x+h/2)
        k2z=h*dzdx(z,y+k1y/2,x+h/2)
        k2b=g*dydx(c+k1c/2,a+g/2)
        k2c=g*dzdx(c,b+k1b/2,a+g/2)

        k3y=h*dydx(z+k2z/2,x+h/2)
        k3z=h*dzdx(z,y+k2y/2,x+h/2)
        k3b=g*dydx(c+k2c/2,a+g/2)
        k3c=g*dzdx(c,b+k2b/2,a+g/2)

        k4y=h*dydx(z+k3z/2,x+h)
        k4z=h*dzdx(z,y+k3y/2,x+h)
        k4b=g*dydx(c+k3c/2,a+g)
        k4c=g*dzdx(c,b+k3b/2,a+g)

        y+=(k1y + 2*k2y + 2*k3y + k4y)/6
        b+=(k1b + 2*k2b + 2*k3b + k4b)/6
        z+=(k1z + 2*k2z + 2*k3z + k4z)/6
        c+=(k1c + 2*k2c + 2*k3c + k4c)/6
        x+=h
        a+=g
        if (flag==True): file.writelines([str(a)+"   ",str(b)+"\n"])
    if (flag==True): file.close()
    return y
#function for shooting method using RK4
def boundary(a,b,ya,yb,guess,dydx,dzdx,h,tol):
    #bracketing
    x=rungekutta4(a,ya,guess,dydx,dzdx,h,a,b,False)
    Ch=guess
    i=1
    if (x==yb): rungekutta4(a,ya,guess,dydx,dzdx,h,a,b,False)
    else:
        while(x>yb):#Bracketing y(b)
            guess-= i
            x=rungekutta4(a,ya,guess,dydx,dzdx,h,a,b,False)
        #Converging to y(b) using lagrange extrapolation
        Cl=guess
        ych=rungekutta4(a,ya,Ch,dydx,dzdx,h,a,b,False)
        ycl=rungekutta4(a,ya,Cl,dydx,dzdx,h,a,b,False)
        yc=yb+10*abs(tol)
        while (abs(yc-yb)>=tol):
            C=Cl+(Ch-Cl)*(yb-ycl)/(ych-ycl)
            yc=rungekutta4(a,ya,C,dydx,dzdx,h,a,b,False)
            if (yc>yb) :
                Ch=C
                ych=rungekutta4(a,ya,Ch,dydx,dzdx,h,a,b,False)
            else :
                Cl=C
                ycl=rungekutta4(a,ya,Cl,dydx,dzdx,h,a,b,False)
    return C
    






