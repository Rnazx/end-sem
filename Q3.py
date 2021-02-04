#first part
import math
T=[0.00,0.30,0.60,0.90,1.20,1.50,1.80,2.10,2.40,2.70,3.00,3.30]
W=[2.2,1.96,1.72,1.53,1.36,1.22,1.10,1.00,0.86,0.75,0.65,0.60]
T2=[]
TW=[]
N=len(T)
#define sum
def sum(X):
    s=0
    for i in range(N):
        s+=X[i]
    return s

for i in range(N):
    T2.append((T[i])**2)
    TW.append(T[i]*W[i])
#slope using linear fit
slope=(sum(T)*sum(W)-N*sum(TW))/((sum(T))**2-N*sum(T2))
intercept=(sum(T)*sum(TW)-(sum(T)*sum(W)))/((sum(T)**2)-N*sum(T2))
print("The slope and the intercept obtained are",slope,intercept)
#second part
logw=[]
tlogw=[]
for i in range(N):
    logw.append(math.log(W[i]))
    tlogw.append(T[i]*math.log(W[i]))
wc=(N*sum(tlogw)-sum(T)*sum(logw))/((sum(T)**2)-N*sum(T2)) 
s=(sum(X2)*sum(logw))

print("The wc and the A obtained are",wc)