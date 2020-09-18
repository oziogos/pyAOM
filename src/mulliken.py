import math
import numpy as np
def equal(A,B,equal_diff=1.0e-6):
    if abs(A-B)<equal_diff:
        value=1
    else:
        value=0
    return value
def Afunction(k,p):
    prefac=math.exp(-p)
    mysum=0
    for mu in range(1,k+1+1):
        nom=math.factorial(k)
        denom=pow(p,mu)
        n=k-mu+1
        denom=denom*math.factorial(n)
        mysum+=nom/denom
    res=prefac*mysum
    return res
def Bfunction(k,p,t):
    prefac1=-math.exp(-p*t)
    prefac2=-math.exp(p*t)
    sum1=0.0;sum2=0.0
    for mu in range(1,k+1+1):
        nom=math.factorial(k)
        u=p*t
        denom=pow(u,mu)
        n=k-mu+1
        denom=denom*math.factorial(n)
        u=nom/denom
        sum1=sum1+u
        n=k-mu
        sum2=sum2+u*pow(-1,n)
    res=prefac1*sum1+prefac2*sum2
    return res
def pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2):
    norm=(X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)
    norm=math.sqrt(norm)
    res=0.5*norm*(mu1+mu2)
    return res
def tvalue(mu1,mu2):
    res=mu1-mu2
    res=res/(mu1+mu2)
    return res
def Smulliken_s_s(p,t,type1,type2):
    res=0
    if type1==1 and type2==1: #1s-1s
        if equal(t,0.0)==1:
            res=(1.0/6.0)*(p*p*p)*(3.0*Afunction(2,p)-Afunction(0,p))
        else:
            res=(1.0/4.0)*(p*p*p)*pow((1.0-t*t),1.5)*(Afunction(2,p)*Bfunction(0,p,t)-Afunction(0,p)*Bfunction(2,p,t))
    if type1==2 and type2==2: #2s-2s
        if equal(t,0.0)==1:
            res=(1.0/360.0)*(p*p*p*p*p)*(15.0*Afunction(4,p)-10.0*Afunction(2,p)+3.0*Afunction(0,p))
        else:
            res=(1.0/48.0)*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Afunction(4,p)*Bfunction(0,p,t)-2.0*Afunction(2,p)*Bfunction(2,p,t)+Afunction(0,p)*Bfunction(4,p,t))
    if (type1==1 and type2==2) or (type1==2 and type2==1): #1s-2s
        if equal(t,0.0)==1:
            res=(1.0/12.0)*pow(3.0,-0.5)*(p*p*p*p)*(3.0*Afunction(3,p)-Afunction(1,p))
        else:
            res=(1.0/8.0)*pow(3.0,-0.5)*(p*p*p*p)*pow((1.0+t),1.5)*pow((1.0-t),2.5)*(Afunction(3,p)*Bfunction(0,p,t)-Afunction(2,p)*Bfunction(1,p,t)-Afunction(1,p)*Bfunction(2,p,t)+Afunction(0,p)*Bfunction(3,p,t))
    if (type1==1 and type2==6) or (type1==6 and type2==1): #1s-3s
        if equal(t,0.0)==1:
            res=(1.0/60.0)*(1.0/math.sqrt(10.0))*(p*p*p*p*p)*(5.0*Afunction(4,p)-Afunction(0,p))
        else:
            res=(1.0/24.0)*(1.0/math.sqrt(10.0))*(p*p*p*p*p)*pow((1.0+t),1.5)*pow((1.0-t),3.5)*(Afunction(4,p)*Bfunction(0,p,t)-2.0*Afunction(3,p)*Bfunction(1,p,t)+2.0*Afunction(1,p)*Bfunction(3,p,t)-Afunction(0,p)*Bfunction(4,p,t))
    if (type1==2 and type2==6) or (type1==6 and type2==2): #2s-3s
        if equal(t,0.0)==1:
            res=(1.0/360.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p*p)*(15.0*Afunction(5,p)-10.0*Afunction(3,p)+3.0*Afunction(1,p))
        else:
            res=(1.0/48.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(5,p)*Bfunction(0,p,t)-Afunction(4,p)*Bfunction(1,p,t)-2.0*Afunction(3,p)*Bfunction(2,p,t)+2.0*Afunction(2,p)*Bfunction(3,p,t)+Afunction(1,p)*Bfunction(4,p,t)-Afunction(0,p)*Bfunction(5,p,t))
    if type1==6 or type2==6: #3s-3s
        if equal(t,0.0)==1:
            res=(1.0/25200.0)*(p*p*p*p*p*p*p)*(35.0*Afunction(6,p)-35.0*Afunction(4,p)+21.0*Afunction(2,p)-5.0*Afunction(0,p))
        else:
            res=(1.0/1440.0)*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(Afunction(6,p)*Bfunction(0,p,t)-3.0*Afunction(4,p)*Bfunction(2,p,t)+3.0*Afunction(2,p)*Bfunction(4,p,t)-Afunction(0,p)*Bfunction(6,p,t))
    return res;
def Smulliken_s_psigma(p,t,type1,type2):
    res=0
    if (type1==1 and (type2>=3 and type2<=5)) or (type2==1 and (type1>=3 and type1<=5)): #1s-2psigma
        if equal(t,0.0)==1:
            res=(1.0/12.0)*(p*p*p*p)*(3.0*Afunction(2,p)-Afunction(0,p))
        else:
            res=(1.0/8.0)*(p*p*p*p)*pow((1.0+t),1.5)*pow((1.0-t),2.5)*(-1.0*Afunction(3,p)*Bfunction(1,p,t)+Afunction(2,p)*Bfunction(0,p,t)+Afunction(1,p)*Bfunction(3,p,t)-Afunction(0,p)*Bfunction(2,p,t))
    if (type1==2 and (type2>=3 and type2<=5)) or (type2==2 and (type1>=3 and type1<=5)): #2s-2psigma
        if equal(t,0.0)==1:
            res=(1.0/60.0)*(1.0/math.sqrt(3.0))*(p*p*p*p*p)*(5.0*Afunction(3,p)-Afunction(1,p))
        else:
            res=(1.0/16.0)*(1.0/math.sqrt(3.0))*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Afunction(3,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(1,p)*(Bfunction(4,p,t)-Bfunction(2,p,t))+Bfunction(1,p,t)*(Afunction(2,p)-Afunction(4,p))+Bfunction(3,p,t)*(Afunction(2,p)-Afunction(0,p)))
    if (type1==1 and (type2>=7 and type2<=9)) or (type2==1 and (type1>=7 and type1<=9)): #1s-3psigma
        if equal(t,0.0)==1:
            res=(1.0/15.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p)*(5.0*Afunction(3,p)-2.0*Afunction(1,p))
        else:
            res=(1.0/8.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p)*pow((1.0+t),1.5)*pow((1-t),3.5)*(Afunction(3,p)*(Bfunction(0,p,t)+Bfunction(2,p,t))-Afunction(1,p)*(Bfunction(2,p,t)+Bfunction(4,p,t))-Bfunction(1,p,t)*(Afunction(2,p)+Afunction(4,p))+Bfunction(3,p,t)*(Afunction(0,p)+Afunction(2,p)))
    if (type1==2 and (type2>=7 and type2<=9)) or (type2==2 and (type1>=7 and type1<=9)): #2s-3psigma
        if equal(t,0.0)==1:
            res=(1.0/360.0)*(1.0/math.sqrt(10.0))*(p*p*p*p*p*p)*(15.0*Afunction(4,p)-10.0*Afunction(2,p)+3.0*Afunction(0,p))
        else:
            res=(1.0/48.0)*(1.0/math.sqrt(10.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1-t),3.5)*(-1.0*Afunction(5,p)*Bfunction(1,p,t)+Afunction(4,p)*Bfunction(0,p,t)+2.0*Afunction(3,p)*Bfunction(3,p,t)-2.0*Afunction(2,p)*Bfunction(2,p,t)-Afunction(1,p)*Bfunction(5,p,t)+Afunction(0,p)*Bfunction(4,p,t))
    if (type1==6 and (type2>=3 and type2<=5)) or (type2==6 and (type1>=3 and type1<=5)): #3s-2psigma
        if equal(t,0.0)==1:
            res=(1.0/360.0)*(1.0/math.sqrt(10.0))*(p*p*p*p*p*p)*(5.0*Afunction(4,p)+6.0*Afunction(2,p)-3.0*Afunction(0,p))
        else:
            res=(1.0/48.0)*(1.0/math.sqrt(10.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(4,p)*(Bfunction(0,p,t)-2.0*Bfunction(2,p,t))+Afunction(1,p)*(2.0*Bfunction(3,p,t)-Bfunction(5,p,t))+Bfunction(1,p,t)*(Afunction(5,p)-2.0*Afunction(3,p))+Bfunction(4,p,t)*(2.0*Afunction(2,p)-Afunction(0,p)))
    if (type1==6 and (type2>=7 and type2<=9)) or (type2==6 and (type1>=7 and type1<=9)): #3s-3psigma
        if equal(t,0.0)==1:
            res=(1.0/12600.0)*(1.0/math.sqrt(3.0))*(p*p*p*p*p*p*p)*(35.0*Afunction(5,p)-14.0*Afunction(3,p)+3.0*Afunction(1,p))
        else:
            res=(1.0/480.0)*(1.0/math.sqrt(3.0))*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(-1.0*Afunction(6,p)*Bfunction(1,p,t)+Afunction(5,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(4,p)*(Bfunction(1,p,t)+2.0*Bfunction(3,p,t))+2.0*Afunction(3,p)*(Bfunction(4,p,t)-Bfunction(2,p,t))-Afunction(2,p)*(2.0*Bfunction(3,p,t)+Bfunction(5,p,t))+Afunction(1,p)*(Bfunction(4,p,t)-Bfunction(6,p,t))+Afunction(0,p)*Bfunction(5,p,t))
    return res;
def Smulliken_psigma_psigma(p,t,type1,type2):
    res=0;
    if(type1>=3 and type1<=5) and (type2>=3 and type2<=5): #2psigma-2psigma
        if equal(t,0.0)==1:
            res=(1.0/120.0)*(p*p*p*p*p)*(5.0*Afunction(4,p)-18.0*Afunction(2,p)+5.0*Afunction(0,p))
        else:
            res=(1.0/16.0)*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Bfunction(2,p,t)*(Afunction(0,p)+Afunction(4,p))-Afunction(2,p)*(Bfunction(0,p,t)+Bfunction(4,p,t)))
    if (type1>=3 and type1<=5) and (type2>=7 and type2<=9): #2psigma-3psigma
        if equal(t,0.0)==1:
            res=(1.0/120.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p*p)*(5.0*Afunction(5,p)-18.0*Afunction(3,p)+5.0*Afunction(1,p))
        else:
            res=(1.0/16.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(2,p)*(Bfunction(1,p,t)+Bfunction(5,p,t))-Afunction(3,p)*(Bfunction(0,p,t)+Bfunction(4,p,t))-Bfunction(3,p,t)*(Afunction(0,p)+Afunction(4,p))+Bfunction(2,p,t)*(Afunction(1,p)+Afunction(5,p)))
    if (type1>=7 and type1<=9) and (type2>=7 and type2<=9): #3psigma-3psigma
        if equal(t,0.0)==1:
            res=(1.0/25200.0)*(p*p*p*p*p*p*p)*(35.0*Afunction(6,p)-147.0*Afunction(4,p)+85.0*Afunction(2,p)-21.0*Afunction(0,p))
        else:
            res=(1.0/480.0)*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(Afunction(6,p)*Bfunction(2,p,t)-Afunction(4,p)*(Bfunction(0,p,t)+2.0*Bfunction(4,p,t))+Afunction(2,p)*(Bfunction(6,p,t)+2.0*Bfunction(2,p,t))-Afunction(0,p)*Bfunction(4,p,t))
    return res
def Smulliken_ppi_ppi(p,t,type1,type2):
    res=0
    if (type1>=3 and type1<=5) and (type2>=3 and type2<=5): #2ppi-2ppi
        if equal(t,0.0)==1:
            res=(1.0/120.0)*(p*p*p*p*p)*(5.0*Afunction(4,p)-6.0*Afunction(2,p)+Afunction(0,p))
        else:
            res=(1.0/32.0)*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Afunction(4,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(2,p)*(Bfunction(4,p,t)-Bfunction(0,p,t))+Afunction(0,p)*(Bfunction(2,p,t)-Bfunction(4,p,t)))
    if (type1>=3 and type1<=5) and (type2>=7 and type2<=9): #2ppi-3ppi
        if equal(t,0.0)==1:
            res=(1.0/120.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p*p)*(5.0*Afunction(5,p)-6.0*Afunction(3,p)+Afunction(1,p))
        else:
            res=(1.0/32.0)*(1.0/math.sqrt(30.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(5,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(4,p)*(Bfunction(3,p,t)-Bfunction(1,p,t))+Afunction(3,p)*(Bfunction(4,p,t)-Bfunction(0,p,t))+Afunction(2,p)*(Bfunction(1,p,t)-Bfunction(5,p,t))+Afunction(1,p)*(Bfunction(2,p,t)-Bfunction(4,p,t))+Afunction(0,p)*(Bfunction(5,p,t)-Bfunction(3,p,t)))
    if (type1>=7 and type1<=9) and (type2>=7 and type2<=9): #3ppi-3ppi
        if equal(t,0.0)==1:
            res=(1.0/25200.0)*(p*p*p*p*p*p*p)*(35.0*Afunction(6,p)-49.0*Afunction(4,p)+17.0*Afunction(2,p)-3.0*Afunction(0,p))
        else:
            res=(1.0/960.0)*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(Afunction(6,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(4,p)*(2.0*Bfunction(4,p,t)-Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(2,p)*(2.0*Bfunction(2,p,t)-Bfunction(4,p,t)-Bfunction(6,p,t))+Afunction(0,p)*(Bfunction(6,p,t)-Bfunction(4,p,t)))
    return res
def overlap(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2,type1,type2):
    res=0
    vector=[[1,0,0],[0,1,0],[0,0,1]]
    # s-s
    if type1==1 and type2==1:
        res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
    if type1==2 and type2==2:
        res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
    if type1==6 and type2==6:
        res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
    if type1==1 and type2==2:
        res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
    if type1==2 and type2==1:
        res=Smulliken_s_s(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
    if type1==1 and type2==6:
        res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
    if type1==6 and type2==1:
        res=Smulliken_s_s(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
    if type1==2 and type2==6:
        res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
    if type1==6 and type2==2:
        res=Smulliken_s_s(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
    # s-p
    if type1==1 and (type2>=3 and type2<=5):

        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        res=vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2]
        res=res*S
    
    if type2==1 and (type1>=3 and type1<=5):

        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        res=vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2]
        res=res*S
    
    if type1==2 and (type2>=3 and type2<=5):

        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        res=vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2]
        res=res*S
    
    if type2==2 and (type1>=3 and type1<=5):

        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        res=vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2]
        res=res*S
    
    if type1==1 and (type2>=7 and type2<=9):

        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        res=vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2]
        res=res*S
    
    if type2==1 and (type1>=7 and type1<=9):

        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        res=vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2]
        res=res*S
    
    if type1==2 and (type2>=7 and type2<=9):

        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        res=vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2]
        res=res*S
    
    if type2==2 and (type1>=7 and type1<=9):

        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        res=vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2]
        res=res*S
    
    if (type1>=3 and type1<=5) and type2==6:

        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        res=vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2]
        res=res*S
    
    if (type2>=3 and type2<=5) and type1==6:

        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        res=vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2]
        res=res*S
    
    if type1==6 and (type2>=7 and type2<=9):

        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        res=vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2]
        res=res*S
    
    if type2==6 and (type1>=7 and type1<=9):

        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        res=vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2]
        res=res*S
    
    if (type1>=3 and type1<=5) and (type2>=3 and type2<=5):

        S=Smulliken_psigma_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        Spi=Smulliken_ppi_ppi(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        r1=(vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2])
        r1=r1*(vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2])
        r2=(vector[type1-2-1][0]*eyb[0]+vector[type1-2-1][1]*eyb[1]+vector[type1-2-1][2]*eyb[2])
        r2=r2*(vector[type2-2-1][0]*eyb[0]+vector[type2-2-1][1]*eyb[1]+vector[type2-2-1][2]*eyb[2])
        r3=(vector[type1-2-1][0]*exb[0]+vector[type1-2-1][1]*exb[1]+vector[type1-2-1][2]*exb[2])
        r3=r3*(vector[type2-2-1][0]*exb[0]+vector[type2-2-1][1]*exb[1]+vector[type2-2-1][2]*exb[2])
        res=r1*S+(r2+r3)*Spi
    
    if (type1>=3 and type1<=5) and (type2>=7 and type2<=9):

        S=Smulliken_psigma_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        Spi=Smulliken_ppi_ppi(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        r1=   (vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2])
        r1=r1*(vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2])
        r2=   (vector[type1-2-1][0]*eyb[0]+vector[type1-2-1][1]*eyb[1]+vector[type1-2-1][2]*eyb[2])
        r2=r2*(vector[type2-6-1][0]*eyb[0]+vector[type2-6-1][1]*eyb[1]+vector[type2-6-1][2]*eyb[2])
        r3=   (vector[type1-2-1][0]*exb[0]+vector[type1-2-1][1]*exb[1]+vector[type1-2-1][2]*exb[2])
        r3=r3*(vector[type2-6-1][0]*exb[0]+vector[type2-6-1][1]*exb[1]+vector[type2-6-1][2]*exb[2])
        res=r1*S+(r2+r3)*Spi
    
    if (type2>=3 and type2<=5) and (type1>=7 and type1<=9):

        S=Smulliken_psigma_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        Spi=Smulliken_ppi_ppi(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1)
        exb,eyb,ezb=MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1)
        r1=   (vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2])
        r1=r1*(vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2])
        r2=   (vector[type2-2-1][0]*eyb[0]+vector[type2-2-1][1]*eyb[1]+vector[type2-2-1][2]*eyb[2])
        r2=r2*(vector[type1-6-1][0]*eyb[0]+vector[type1-6-1][1]*eyb[1]+vector[type1-6-1][2]*eyb[2])
        r3=   (vector[type2-2-1][0]*exb[0]+vector[type2-2-1][1]*exb[1]+vector[type2-2-1][2]*exb[2])
        r3=r3*(vector[type1-6-1][0]*exb[0]+vector[type1-6-1][1]*exb[1]+vector[type1-6-1][2]*exb[2])
        res=r1*S+(r2+r3)*Spi
    
    if (type1>=7 and type1<=9) and (type2>=7 and type2<=9):

        S=Smulliken_psigma_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        Spi=Smulliken_ppi_ppi(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2)
        exb,eyb,ezb=MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2)
        r1=   (vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2])
        r1=r1*(vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2])
        r2=   (vector[type1-6-1][0]*eyb[0]+vector[type1-6-1][1]*eyb[1]+vector[type1-6-1][2]*eyb[2])
        r2=r2*(vector[type2-6-1][0]*eyb[0]+vector[type2-6-1][1]*eyb[1]+vector[type2-6-1][2]*eyb[2])
        r3=   (vector[type1-6-1][0]*exb[0]+vector[type1-6-1][1]*exb[1]+vector[type1-6-1][2]*exb[2])
        r3=r3*(vector[type2-6-1][0]*exb[0]+vector[type2-6-1][1]*exb[1]+vector[type2-6-1][2]*exb[2])
        res=r1*S+(r2+r3)*Spi
    
    return res

def MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2):
    align=0
    ez=np.array([X2-X1,Y2-Y1,Z2-Z1])
    ez=ez/np.linalg.norm(ez)
    if equal(ez[1],0.0)==1 and equal(ez[2],0.0)==1:
        ey=np.array([0,0,1])
        align=1
    if equal(ez[0],0.0)==1 and equal(ez[2],0.0)==1:
        ey=np.array([1,0,0])
        align=1
    if equal(ez[0],0.0)==1 and equal(ez[1],0.0)==1:
        ey=np.array([0,1,0])
        align=1
    if align==0:
        t=[0,0,0]
        maxpos=0
        absmax=abs(ez[0])
        for i in [1,2]:
            if abs(ez[i])>absmax:
                maxpos=i
                absmax=abs(ez[i])
        for i in [0,1,2]:
            if i!=maxpos:
                t[i]=1.0
            else:
                t[i]=-1.0*(ez[(i+1)%3+1-1]+ez[(i+1+1)%3+1-1])/ez[i]
        ey=np.array(t)
        ey=ey/np.linalg.norm(ey)
    ex=np.cross(ey,ez)
    return ex,ey,-ez
