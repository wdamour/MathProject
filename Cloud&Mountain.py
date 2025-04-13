import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

def sqrt(n):
    return n**(1/2)

# Permet de transposer une matrice dans le but de tracer une forme 
def transpose(A):  
    C=[] 
    for j in range (len(A[0])): 
        B=[] 
        for i in range (len(A)): 
            B.append(A[i][j]) 
        C.append(B) 
    return C
 
# Génère des plans
def plan(n,xmin,xmax,ymin,ymax, Xtime, Ytime, Ztime, AffineTime) : 
    r=int(n**(1/2)) 
    dx=(xmax-xmin)/(r+1) 
    dy=(ymax-ymin)/(r+1) 
    W=[] 
    for i in range(r+1): 
        for j in range(r+1): 
            x=xmin+dx*i 
            y=ymin+dy*j 
            W.append([x,y,-(Xtime/Ztime)*x - (Ytime/Ztime)*y - AffineTime/Ztime]) 
    return(transpose(W)) 


def projection(A,C, S, Xtime, Ytime, Ztime, B):
    x_a, y_a, z_a = A
    x_b, y_b, z_b = B
    x_s, y_s, z_s = S
    x_c, y_c, z_c = C
    
    d_x = x_c - x_s 
    d_y = y_c - y_s 
    d_z = z_c - z_s 
    
    c = Xtime * x_b + Ytime * y_b + Ztime * z_b
    
    k = (-Xtime * x_a - Ytime * y_a - Ztime * z_a + c)/(Xtime * d_x + Ytime * d_y + Ztime * d_z)
    
    return (d_x * k + x_a, d_y * k + y_a, d_z * k + z_a)

def homothetieRatio(A,B,Ap,Bp):
    x_a, y_a, z_a = A
    x_b, y_b, z_b = B
        
    x_ap, y_ap, z_ap = Ap
    x_bp, y_bp, z_bp = Bp
    
    l = sqrt((x_b-x_a)**2 + (y_b-y_a)**2 +(z_b-z_a)**2)
    lprime = sqrt((x_bp-x_ap)**2 + (y_bp-y_ap)**2 +(z_bp-z_ap)**2)
    
    return lprime/l
    
def homothetie(A, ratio, B):
    x_a, y_a, z_a = A
    x_b, y_b, z_b = B
    return (x_a * ratio - x_b*(ratio -1), y_a * ratio - y_b*(ratio -1), z_a * ratio - z_b*(ratio -1))

# N nb of point of each plain = 1000,
# Soleil coordiantes (4, 4, 7)
# PlaneTime = [XtimePlane, YTimePlane, ZtimePlane, AffinePlane] 2 3 5 -6
# NuageTime = [XtimeNuage, YTimeNuage, ZtimeNuage, AffineNuage] 2 3 5 -25
#LimitPlane = [xmin,xmax,ymin,ymax] -5 5 -5 5
#LimitNuage = [xmin,xmax,ymin,ymax] -1 1 -1 1
def ombre(N, Soleil, PlaneTime, NuageTime,LimitPlane, LimitNuage):
    (X,Y,Z)=plan(N, LimitPlane[0],LimitPlane[1],LimitPlane[2], LimitPlane[3], PlaneTime[0], PlaneTime[1], PlaneTime[2], PlaneTime[3])  #Mountain
    (X1,Y1,Z1)=plan(N, LimitNuage[0],LimitNuage[1],LimitNuage[2], LimitNuage[3], NuageTime[0], NuageTime[1], NuageTime[2], NuageTime[3])  #Nuage
    
    center_x = (LimitNuage[1] + LimitNuage[0])/2
    center_y = (LimitNuage[3] + LimitNuage[2])/2
    CenterNuage = [center_x, center_y,-( NuageTime[0]/NuageTime[2])*center_x - (NuageTime[1]/ NuageTime[2])*center_y -  NuageTime[3]/ NuageTime[2]]
    BordNuage = [LimitNuage[0], LimitNuage[2],-( NuageTime[0]/NuageTime[2])*LimitNuage[0] - (NuageTime[1]/ NuageTime[2])*LimitNuage[2] -  NuageTime[3]/ NuageTime[2]]
    
    
    XCenterPrime, YCenterPrime, ZCenterPrime = projection(CenterNuage, CenterNuage, Soleil, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
    XBordPrime, YBordPrime, ZBordPrime = projection(BordNuage, BordNuage, Soleil, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
    
    XOmbreTemp, YOmbreTemp, ZOmbreTemp = [],[],[]
    for i in range (len(X1)):
        (Xtemp, Ytemp, Ztemp) = projection([X1[i], Y1[i], Z1[i]], CenterNuage, Soleil, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
        XOmbreTemp.append(Xtemp)
        YOmbreTemp.append(Ytemp)
        ZOmbreTemp.append(Ztemp)
        
    ratio = homothetieRatio(CenterNuage, BordNuage, [XCenterPrime, YCenterPrime, ZCenterPrime], [XBordPrime, YBordPrime, ZBordPrime])
    XOmbre, YOmbre, ZOmbre = [], [], []
    for i in range (len(X1)):
        (Xtemp, Ytemp, Ztemp) = homothetie([XOmbreTemp[i], YOmbreTemp[i], ZOmbreTemp[i]], ratio, [XCenterPrime, YCenterPrime, ZCenterPrime])
        XOmbre.append(Xtemp)
        YOmbre.append(Ytemp)
        ZOmbre.append(Ztemp)
        
    fig = plt.figure("Ombre d’un nuage")
    axes = fig.add_subplot(111, projection="3d")
    print(axes, type(axes)) 

    # Trace les listes en 3D 
    axes.plot(X,Y,Z, label="Mountain") 
    axes.plot(X1,Y1,Z1, label="Cloud") 
    axes.plot(XOmbre, YOmbre, ZOmbre, label="Projection & Homotethie") 
    axes.plot(XOmbreTemp, YOmbreTemp, ZOmbreTemp, label="Projection") 
    
    axes.scatter(Soleil[0], Soleil[1], Soleil[2], color='yellow', s=50, label="Sun")
    axes.scatter(CenterNuage[0], CenterNuage[1], CenterNuage[2], color='red', s=50, label="MiddleCloud")
    axes.scatter(BordNuage[0], BordNuage[1], BordNuage[2], color='black', s=50, label="CornerCloud")
    axes.scatter(XCenterPrime, YCenterPrime, ZCenterPrime, color='red', s=50, label="MiddleCloudProjection")
    axes.scatter(XBordPrime, YBordPrime, ZBordPrime, color='black', s=50, label="CornerCloudProjection")

    # Ajoute des étiquettes pour les axes axes.set_xlabel(« X ») 
    axes.set_ylabel("Y") 
    axes.set_zlabel("Z") 
    axes.legend(bbox_to_anchor=(1, 1),bbox_transform=fig.transFigure)
    plt.show() 

ombre(1000, [4, 4, 7], [2, 3, 5, -6], [2, 3, 5, -25], [-5, 5, -5, 5], [-1, 1, -1, 1] )