import math
res=0
def anIres(X1,Y1,Z1,a1,lx1,ly1,lz1,X2,Y2,Z2,a2,lx2,ly2,lz2):
    global res
    try:
        if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
            res=(2*math.sqrt(2)*pow(a1,0.75)*pow(a2,0.75))/(pow(a1+a2,1.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
            res=(4*math.sqrt(2)*pow(a1,1.75)*pow(a2,1.25)*(Z1-Z2))/(pow(a1+a2,2.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
            res=(4*math.sqrt(0.6666666666666666)*pow(a1,0.75)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(Z1-Z2,2)))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
            res=(4*math.sqrt(2)*pow(a1,1.75)*pow(a2,1.25)*(Y1-Y2))/(pow(a1+a2,2.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,2.75)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
            res=(4*math.sqrt(0.6666666666666666)*pow(a1,0.75)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(Y1-Y2,2)))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
            res=(4*math.sqrt(2)*pow(a1,1.75)*pow(a2,1.25)*(X1-X2))/(pow(a1+a2,2.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,2.75)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,2.75)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Y1-Y2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
            res=(4*math.sqrt(0.6666666666666666)*pow(a1,0.75)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(X1-X2,2)))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==0:
            res=(4*math.sqrt(2)*pow(a1,1.25)*pow(a2,1.75)*(-Z1+Z2))/(pow(a1+a2,2.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==1:
            res=(4*math.sqrt(2)*pow(a1,1.25)*pow(a2,1.25)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==2:
            res=(-8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(Z1-Z2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,2.25)*pow(a2,2.25)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,2.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(Y1-Y2)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==2 and lz2==0:
            res=(-8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(Y1-Y2,2))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,2.25)*pow(a2,2.25)*(X1-X2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,2.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==1 and ly2==1 and lz2==0:
            res=(-16*math.sqrt(2)*pow(a1,3.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==1 and lx2==2 and ly2==0 and lz2==0:
            res=(-8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(X1-X2,2))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==0 and lz2==0:
            res=(4*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*(a1+a2*(1+2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==0 and lz2==1:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(Z1-Z2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==0 and lz2==2:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-6*a1*a2*(-1+a2*pow(Z1-Z2,2))+pow(a2,2)*(3+2*a2*pow(Z1-Z2,2))+pow(a1,2)*(3-6*a2*pow(Z1-Z2,2)+4*pow(a2,2)*pow(Z1-Z2,4))+2*pow(a1,3)*pow(Z1-Z2,2)))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==1 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(Y1-Y2)*(a1+a2*(1+2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==1 and lz2==1:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(Y1-Y2)*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(Z1-Z2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==2 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(Y1-Y2,2))*(a1+a2*(1+2*a2*pow(Z1-Z2,2))))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==1 and ly2==0 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(X1-X2)*(a1+a2*(1+2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==1 and ly2==0 and lz2==1:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(Z1-Z2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==1 and ly2==1 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,2.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Y1-Y2)*(a1+a2*(1+2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==0 and lz1==2 and lx2==2 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(X1-X2,2))*(a1+a2*(1+2*a2*pow(Z1-Z2,2))))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
            res=(4*math.sqrt(2)*pow(a1,1.25)*pow(a2,1.75)*(-Y1+Y2))/(pow(a1+a2,2.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
            res=(-8*math.sqrt(2)*pow(a1,2.25)*pow(a2,2.25)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(-Y1+Y2)*(a1+a2+2*pow(a1,2)*pow(Z1-Z2,2)))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
            res=(4*math.sqrt(2)*pow(a1,1.25)*pow(a2,1.25)*(a2+a1*(1-2*a2*pow(Y1-Y2,2))))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,2.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(Y1-Y2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
            res=(-8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(Y1-Y2,2)))*(Y1-Y2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,2.25)*pow(a2,2.25)*(X1-X2)*(Y1-Y2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
            res=(16*math.sqrt(2)*pow(a1,3.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(-Y1+Y2)*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,2.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a2+a1*(1-2*a2*pow(Y1-Y2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(X1-X2,2))*(-Y1+Y2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.75)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==0 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.25)*(-Y1+Y2)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==0 and lz2==2:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(Y1-Y2)*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(Z1-Z2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==1 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.25)*(a2+a1*(1-2*a2*pow(Y1-Y2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==1 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(Y1-Y2,2)))*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==2 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(Y1-Y2,2)))*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==1 and ly2==0 and lz2==0:
            res=(16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,3.25)*(X1-X2)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==1 and ly2==0 and lz2==1:
            res=(16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(-Y1+Y2)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==1 and ly2==1 and lz2==0:
            res=(-16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a2+a1*(1-2*a2*pow(Y1-Y2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==1 and lz1==1 and lx2==2 and ly2==0 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(X1-X2,2))*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
            res=(4*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*(a1+a2*(1+2*a2*pow(Y1-Y2,2))))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(a1+a2*(1+2*a2*pow(Y1-Y2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2*(1+2*a2*pow(Y1-Y2,2)))*(a1+a2+2*pow(a1,2)*pow(Z1-Z2,2)))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(Y1-Y2,2)))*(Y1-Y2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(Y1-Y2,2)))*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-6*a1*a2*(-1+a2*pow(Y1-Y2,2))+pow(a2,2)*(3+2*a2*pow(Y1-Y2,2))+pow(a1,2)*(3-6*a2*pow(Y1-Y2,2)+4*pow(a2,2)*pow(Y1-Y2,4))+2*pow(a1,3)*pow(Y1-Y2,2)))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(X1-X2)*(a1+a2*(1+2*a2*pow(Y1-Y2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,2.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a1+a2*(1+2*a2*pow(Y1-Y2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(X1-X2)*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(Y1-Y2,2)))*(Y1-Y2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==0 and ly1==2 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2+2*pow(a1,2)*pow(X1-X2,2))*(a1+a2*(1+2*a2*pow(Y1-Y2,2))))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
            res=(-4*math.sqrt(2)*pow(a1,1.25)*pow(a2,1.75)*(X1-X2))/(pow(a1+a2,2.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
            res=(-8*math.sqrt(2)*pow(a1,2.25)*pow(a2,2.25)*(X1-X2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(-X1+X2)*(a1+a2+2*pow(a1,2)*pow(Z1-Z2,2)))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,2.25)*pow(a2,2.25)*(-X1+X2)*(Y1-Y2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
            res=(-16*math.sqrt(2)*pow(a1,3.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
            res=(-8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a1+a2+2*pow(a1,2)*pow(Y1-Y2,2)))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
            res=(4*math.sqrt(2)*pow(a1,1.25)*pow(a2,1.25)*(a2+a1*(1-2*a2*pow(X1-X2,2))))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,2.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,2.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-a2+a1*(-1+2*a2*pow(X1-X2,2)))*(Y1-Y2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
            res=(-8*math.sqrt(0.6666666666666666)*pow(a1,1.25)*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(X1-X2,2)))*(X1-X2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.75)*(X1-X2)*(Z1-Z2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==1:
            res=(-8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.25)*(X1-X2)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==2:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(Z1-Z2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==0:
            res=(16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,3.25)*(X1-X2)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==1:
            res=(-16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Y1-Y2)*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==2 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a1+a2+2*pow(a1,2)*pow(Y1-Y2,2))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.25)*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==1:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(a2+a1*(1-2*a2*pow(Z1-Z2,2))))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==1 and ly2==1 and lz2==0:
            res=(-16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==0 and lz1==1 and lx2==2 and ly2==0 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(X1-X2,2)))*(X1-X2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.75)*(X1-X2)*(Y1-Y2))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
            res=(16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,3.25)*(X1-X2)*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.75)*math.sqrt(pow(a2,2))*(X1-X2)*(Y1-Y2)*(a1+a2+2*pow(a1,2)*pow(Z1-Z2,2)))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
            res=(-8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.25)*(X1-X2)*(a2+a1*(1-2*a2*pow(Y1-Y2,2))))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
            res=(-16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(a2+a1*(1-2*a2*pow(Y1-Y2,2)))*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(X1-X2)*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(Y1-Y2,2)))*(Y1-Y2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,2.25)*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(-Y1+Y2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
            res=(16*math.sqrt(2)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(-Y1+Y2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a2+a1*(1-2*a2*pow(X1-X2,2)))*(a2+a1*(1-2*a2*pow(Y1-Y2,2))))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==1 and ly1==1 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.75)*math.sqrt(pow(a2,2))*(-(a1*a2)+pow(a2,2)+2*pow(a1,2)*(-1+a2*pow(X1-X2,2)))*(X1-X2)*(Y1-Y2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
            res=(4*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*(a1+a2*(1+2*a2*pow(X1-X2,2))))/(pow(a1+a2,3.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(a1+a2*(1+2*a2*pow(X1-X2,2)))*(Z1-Z2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2*(1+2*a2*pow(X1-X2,2)))*(a1+a2+2*pow(a1,2)*pow(Z1-Z2,2)))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(a1+a2*(1+2*a2*pow(X1-X2,2)))*(Y1-Y2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,2.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2*(1+2*a2*pow(X1-X2,2)))*(Y1-Y2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(a1+a2*(1+2*a2*pow(X1-X2,2)))*(a1+a2+2*pow(a1,2)*pow(Y1-Y2,2)))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
            res=(8*math.sqrt(0.6666666666666666)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,1.25)*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(X1-X2,2)))*(X1-X2))/(pow(a1+a2,4.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(X1-X2,2)))*(X1-X2)*(Z1-Z2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
            res=(16*math.sqrt(0.6666666666666666)*pow(a1,1.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(pow(a1,2)-2*pow(a2,2)+a1*a2*(-1+2*a2*pow(X1-X2,2)))*(X1-X2)*(Y1-Y2))/(pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
        if lx1==2 and ly1==0 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
            res=(8*math.sqrt(2)*pow(a1,0.75)*math.sqrt(pow(a1,2))*pow(a2,0.75)*math.sqrt(pow(a2,2))*(-6*a1*a2*(-1+a2*pow(X1-X2,2))+pow(a2,2)*(3+2*a2*pow(X1-X2,2))+pow(a1,2)*(3-6*a2*pow(X1-X2,2)+4*pow(a2,2)*pow(X1-X2,4))+2*pow(a1,3)*pow(X1-X2,2)))/(3.*pow(a1+a2,5.5)*math.exp((a1*a2*(pow(X1,2)-2*X1*X2+pow(X2,2)+pow(Y1,2)-2*Y1*Y2+pow(Y2,2)+pow(Z1,2)-2*Z1*Z2+pow(Z2,2)))/(a1+a2)))
    except OverflowError as err:
        pass
    return res
