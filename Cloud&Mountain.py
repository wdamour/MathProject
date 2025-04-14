import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider

#square root function
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
 
# Génère des plans (Montagne et Nuage)
# n -> number of point to generate
# xmin -> limite min on x axis
# xmax -> limite max on x axis
# ymin -> limite min on y axis
# ymax -> limite max on y axis
# (Xtime,Ytime,Ztime,AffineTime) -> (a,b,c,d) => ax + by + cz + d = 0
def plan(n:int, xmin:float, xmax:float, ymin:float, ymax:float, Xtime:float, Ytime:float, Ztime:float, AffineTime:float) : 
    r=int(n**(1/2)) 
    dx=(xmax-xmin)/(r+1) 
    dy=(ymax-ymin)/(r+1) 
    W=[] 
    for i in range(r+1): 
        for j in range(r+1): 
            x=xmin+dx*i 
            y=ymin+dy*j 
            if (Ztime == 0): #prevent division by 0
                W.append([x,y,0]) 
            else:  
                W.append([x,y,-(Xtime/Ztime)*x - (Ytime/Ztime)*y - AffineTime/Ztime]) 
    return(transpose(W)) 

# Projects a point on a plane
# A(x,y,z) -> point to project
# C(x,y,z) -> point on the parallel line
# S(x,y,z) -> other point on the parallel line, in our example it's the sun
# (Xtime,Ytime,Ztime) -> (a,b,c) => ax + by + cz + d = 0 of the arrival plane
# B(x,y,z) -> point on the arrival plane
def projection(A:tuple[float, float, float],C:tuple[float, float, float], S:tuple[float, float, float], Xtime: float, Ytime: float, Ztime:float, B: tuple[float, float, float]):
    x_a, y_a, z_a = A
    x_b, y_b, z_b = B
    x_s, y_s, z_s = S
    x_c, y_c, z_c = C
    
    #calculate the director vector
    d_x = x_c - x_s 
    d_y = y_c - y_s 
    d_z = z_c - z_s 
    
    #Constant (because the plane is affine and not linear)
    c = Xtime * x_b + Ytime * y_b + Ztime * z_b
    
    #isolation of k in the x',y' and z' equations
    k = (-Xtime * x_a - Ytime * y_a - Ztime * z_a + c)/(Xtime * d_x + Ytime * d_y + Ztime * d_z)
    
    # director vector * k + (x,y,z)
    return (d_x * k + x_a, d_y * k + y_a, d_z * k + z_a)

# Compute the homothety ratio
# A(x,y,z) -> center of plane 1
# B(x,y,z) -> point of plane 1
# Ap(x,y,z) -> projection of point A on plane 2
# Bp(x,y,z) -> projection of point B on plane 2
def homothetieRatio(A:tuple[float, float, float],B:tuple[float, float, float],Ap:tuple[float, float, float],Bp:tuple[float, float, float]):
    x_a, y_a, z_a = A
    x_b, y_b, z_b = B
        
    x_ap, y_ap, z_ap = Ap
    x_bp, y_bp, z_bp = Bp
    
    l = sqrt((x_b-x_a)**2 + (y_b-y_a)**2 +(z_b-z_a)**2) #norm of AB
    lprime = sqrt((x_bp-x_ap)**2 + (y_bp-y_ap)**2 +(z_bp-z_ap)**2) #norm of ApBp
    
    return lprime/l

#Compute the homothety
# A(x,y,z) -> point to which homothety is applied
# ratio -> homothety ratio
# B(x,y,z) -> center of homothety
def homothetie(A:tuple[float, float, float], ratio:float, B:tuple[float, float, float]):
    x_a, y_a, z_a = A
    x_b, y_b, z_b = B
    return (x_a * ratio - x_b*(ratio -1), y_a * ratio - y_b*(ratio -1), z_a * ratio - z_b*(ratio -1))

#Computes the projection of a cloud on a parallel plane
# N nb of point of each plain,
# Sun coordiantes
# PlaneTime = [XtimePlane, YTimePlane, ZtimePlane, AffinePlane] -> (a,b,c,d) => ax + by + cz + d = 0
# CloudTime = [XtimeCloud, YTimeCloud, ZtimeCloud, AffineCloud] -> (a,b,c,d) => ax + by + cz + d = 0
#LimitPlane = [xmin,xmax,ymin,ymax] 
#LimitCloud = [xmin,xmax,ymin,ymax] 
def ombreParallele(N:int, Sun:tuple[float, float, float], PlaneTime:tuple[float, float, float, float], CloudTime:tuple[float, float, float, float], LimitPlane:tuple[float, float, float, float], LimitCloud:tuple[float, float, float, float]):
    #Generate the two planes
    (X,Y,Z)=plan(N, LimitPlane[0],LimitPlane[1],LimitPlane[2], LimitPlane[3], PlaneTime[0], PlaneTime[1], PlaneTime[2], PlaneTime[3])  #Mountain
    (X1,Y1,Z1)=plan(N, LimitCloud[0],LimitCloud[1],LimitCloud[2], LimitCloud[3], CloudTime[0], CloudTime[1], CloudTime[2], CloudTime[3])  #Cloud
    
    #Compute the center of the cloud and a corner of the cloud
    center_x = (LimitCloud[1] + LimitCloud[0])/2
    center_y = (LimitCloud[3] + LimitCloud[2])/2
    CenterCloud = [center_x, center_y,-( CloudTime[0]/CloudTime[2])*center_x - (CloudTime[1]/ CloudTime[2])*center_y -  CloudTime[3]/ CloudTime[2]]
    BordCloud = [LimitCloud[0], LimitCloud[2],-( CloudTime[0]/CloudTime[2])*LimitCloud[0] - (CloudTime[1]/ CloudTime[2])*LimitCloud[2] -  CloudTime[3]/ CloudTime[2]]
    
    #Project the center and a corner of the cloud based on their own line with the sun
    XCenterPrime, YCenterPrime, ZCenterPrime = projection(CenterCloud, CenterCloud, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
    XBordPrime, YBordPrime, ZBordPrime = projection(BordCloud, BordCloud, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
    
    #project every point of the cloud based on the line between the sun and the center of the cloud
    XShadowTemp, YShadowTemp, ZShadowTemp = [],[],[]
    for i in range (len(X1)):
        (Xtemp, Ytemp, Ztemp) = projection([X1[i], Y1[i], Z1[i]], CenterCloud, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
        XShadowTemp.append(Xtemp)
        YShadowTemp.append(Ytemp)
        ZShadowTemp.append(Ztemp)
        
    #compute the homothety ratio
    ratio = homothetieRatio(CenterCloud, BordCloud, [XCenterPrime, YCenterPrime, ZCenterPrime], [XBordPrime, YBordPrime, ZBordPrime])
    
    # Compute the homothety of each projected points
    XShadow, YShadow, ZShadow = [], [], []
    for i in range (len(X1)):
        (Xtemp, Ytemp, Ztemp) = homothetie([XShadowTemp[i], YShadowTemp[i], ZShadowTemp[i]], ratio, [XCenterPrime, YCenterPrime, ZCenterPrime])
        XShadow.append(Xtemp)
        YShadow.append(Ytemp)
        ZShadow.append(Ztemp)
    
    #create the window and the 3D plane
    fig = plt.figure("Shadow of a Cloud")
    axes = fig.add_subplot(111, projection="3d")

    # plot the different planes 
    axes.plot(X,Y,Z, color='green', label="Mountain") 
    axes.plot(X1,Y1,Z1, color='gray', label="Cloud") 
    axes.plot(XShadow, YShadow, ZShadow,color='black', label="Projection & Homotethie") 
    axes.plot(XShadowTemp, YShadowTemp, ZShadowTemp,color='darkblue', label="Projection") 
    
    axes.scatter(Sun[0], Sun[1], Sun[2], color='yellow', s=50, label="Sun") #coordinates of the sun
    axes.scatter(CenterCloud[0], CenterCloud[1], CenterCloud[2], color='red', s=50, label="Middle of Cloud") # coordinates of the center of the cloud
    axes.scatter(BordCloud[0], BordCloud[1], BordCloud[2], color='black', s=50, label="Corner of Cloud") # coordinates of a corner of a cloud
    axes.scatter(XCenterPrime, YCenterPrime, ZCenterPrime, color='red', s=50, label="Middle Cloud Projection") # coordinates of the projection of the center of the cloud
    axes.scatter(XBordPrime, YBordPrime, ZBordPrime, color='black', s=50, label="Corner Cloud Projection") # coordinates of the projection of a corner of the cloud

    # Label the axis
    axes.set_xlabel("X") 
    axes.set_ylabel("Y") 
    axes.set_zlabel("Z") 
    axes.legend(bbox_to_anchor=(1, 1),bbox_transform=fig.transFigure)
    plt.show() 

# Computes the projection of a cloud to a plane with which it is not parallel
# N nb of point of each plain,
# Sun coordiantes
# PlaneTime = [XtimePlane, YTimePlane, ZtimePlane, AffinePlane] -> (a,b,c,d) => ax + by + cz + d = 0
# CloudTime = [XtimeCloud, YTimeCloud, ZtimeCloud, AffineCloud] -> (a,b,c,d) => ax + by + cz + d = 0
#LimitPlane = [xmin,xmax,ymin,ymax] 
#LimitCloud = [xmin,xmax,ymin,ymax] 
def ombreNonParallele(N:int, Sun:tuple[float, float, float], PlaneTime:tuple[float, float, float, float], CloudTime:tuple[float, float, float, float],LimitPlane:tuple[float, float, float, float], LimitCloud:tuple[float, float, float, float]):
        #Generate the two planes
    (X,Y,Z)=plan(N, LimitPlane[0],LimitPlane[1],LimitPlane[2], LimitPlane[3], PlaneTime[0], PlaneTime[1], PlaneTime[2], PlaneTime[3])  #Mountain
    (X1,Y1,Z1)=plan(N, LimitCloud[0],LimitCloud[1],LimitCloud[2], LimitCloud[3], CloudTime[0], CloudTime[1], CloudTime[2], CloudTime[3])  #Cloud
    
    #Compute the center of the cloud and a corner of the cloud
    center_x = (LimitCloud[1] + LimitCloud[0])/2
    center_y = (LimitCloud[3] + LimitCloud[2])/2
    CenterCloud = [center_x, center_y,-( CloudTime[0]/CloudTime[2])*center_x - (CloudTime[1]/ CloudTime[2])*center_y -  CloudTime[3]/ CloudTime[2]]
    BordCloud = [LimitCloud[0], LimitCloud[2],-( CloudTime[0]/CloudTime[2])*LimitCloud[0] - (CloudTime[1]/ CloudTime[2])*LimitCloud[2] -  CloudTime[3]/ CloudTime[2]]
    
    #Project the center and a corner of the cloud based on their own line with the sun
    XCenterPrime, YCenterPrime, ZCenterPrime = projection(CenterCloud, CenterCloud, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
    XBordPrime, YBordPrime, ZBordPrime = projection(BordCloud, BordCloud, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
    
    #project every point of the cloud based on the line between the sun and the point
    XShadow, YShadow, ZShadow = [],[],[]
    for i in range (len(X1)):
        (Xtemp, Ytemp, Ztemp) = projection([X1[i], Y1[i], Z1[i]], [X1[i], Y1[i], Z1[i]], Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
        XShadow.append(Xtemp)
        YShadow.append(Ytemp)
        ZShadow.append(Ztemp)
        
    #create the window and the 3D plane   
    fig = plt.figure("Shadow of a Cloud")
    axes = fig.add_subplot(111, projection="3d")
    print(axes, type(axes)) 

    # plot the different planes 
    axes.plot(X,Y,Z, color='green', label="Mountain") 
    axes.plot(X1,Y1,Z1, color='gray', label="Cloud") 
    axes.plot(XShadow, YShadow, ZShadow,color='black', label="Projection/Shadow")  
    
    axes.scatter(Sun[0], Sun[1], Sun[2], color='yellow', s=50, label="Sun") #coordinates of the sun
    axes.scatter(CenterCloud[0], CenterCloud[1], CenterCloud[2], color='red', s=50, label="Middle of Cloud") # coordinates of the center of the cloud
    axes.scatter(BordCloud[0], BordCloud[1], BordCloud[2], color='black', s=50, label="Corner of Cloud")  # coordinates of a corner of a cloud
    axes.scatter(XCenterPrime, YCenterPrime, ZCenterPrime, color='red', s=50, label="Middle Cloud Projection") # coordinates of the projection of the center of the cloud
    axes.scatter(XBordPrime, YBordPrime, ZBordPrime, color='black', s=50, label="Corner Cloud Projection") # coordinates of the projection of a corner of the cloud

    # Label the axis
    axes.set_xlabel("X") 
    axes.set_ylabel("Y") 
    axes.set_zlabel("Z") 
    axes.legend(bbox_to_anchor=(1, 1),bbox_transform=fig.transFigure)
    plt.show()
    
    
# Computes the projection of a cloud to a plane with which it is not parallel and with sliders to change 
# the coordinates of the Sun and the equation of the cloud
# N nb of point of each plain,
# Sun coordiantes
# PlaneTime = [XtimePlane, YTimePlane, ZtimePlane, AffinePlane] -> (a,b,c,d) => ax + by + cz + d = 0
# CloudTime = [XtimeCloud, YTimeCloud, ZtimeCloud, AffineCloud] -> (a,b,c,d) => ax + by + cz + d = 0
#LimitPlane = [xmin,xmax,ymin,ymax] 
#LimitCloud = [xmin,xmax,ymin,ymax] 
def ombreMovingSun(N:int, Sun:tuple[float, float, float], PlaneTime:tuple[float, float, float, float], CloudTime:tuple[float, float, float, float],LimitPlane:tuple[float, float, float, float], LimitCloud:tuple[float, float, float, float]):
    #create a window
    fig =plt.figure("Shadow of a Cloud")
    axes = fig.add_subplot(111, projection="3d")
    plt.subplots_adjust(bottom=0.25)
    
    #slider Sun
    ax_slider_x = plt.axes([0.1, 0.1, 0.33, 0.03])
    ax_slider_y = plt.axes([0.1, 0.05, 0.33, 0.03])
    ax_slider_z = plt.axes([0.1, 0.01, 0.33, 0.03])
    SunSliderX = Slider(ax_slider_x, 'Sun X value :', -5.0, 5.0, Sun[0])
    SunSliderY = Slider(ax_slider_y, 'Sun Y value :', -5.0, 5.0, Sun[1])
    SunSliderZ = Slider(ax_slider_z, 'Sun Z value :', -4.0, 7.0, Sun[2])
    
    #slider cloud plane
    ax_slider_a = plt.axes([0.6, 0.15, 0.33, 0.03])
    ax_slider_b = plt.axes([0.6, 0.1, 0.33, 0.03])
    ax_slider_c = plt.axes([0.6, 0.05, 0.33, 0.03])
    ax_slider_d = plt.axes([0.6, 0.01, 0.33, 0.03])
    CloudSliderA = Slider(ax_slider_a, 'Cloud A value :', -5.0, 5.0, CloudTime[0])
    CloudSliderB = Slider(ax_slider_b, 'Cloud b value :', -5.0, 5.0, CloudTime[1])
    CloudSliderC = Slider(ax_slider_c, 'Cloud C value :', -5.0, 5.0, CloudTime[2])
    CloudSliderD = Slider(ax_slider_d, 'Cloud d value :', -25.0, 25.0, CloudTime[3])
    
    # Mountain and Cloud planes
    X, Y, Z = plan(N, *LimitPlane, *PlaneTime)
    X1, Y1, Z1 = plan(N, *LimitCloud, *CloudTime)
    axes.plot(X, Y, Z, color='green', label="Mountain")
    axes.set_title(formatPlaneEquation(*CloudTime))
    cloud_plot, = axes.plot(X1, Y1, Z1, color='gray', label="Cloud")
    
    # Cloud
    cloud_points = [[X1[i], Y1[i], Z1[i]] for i in range(len(X1))]

    # Sun initial
    sun_scatter = axes.scatter(*Sun, color='yellow', s=50, label="Sun")

    # Computes shadow's initial position
    XShadow, YShadow, ZShadow = [], [], []
    for pt in cloud_points:
        x, y, z = projection(pt, pt, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
        XShadow.append(x)
        YShadow.append(y)
        ZShadow.append(z)

    shadow_plot, = axes.plot(XShadow, YShadow, ZShadow, color='black', label="Shadow")
    
    #changes the value everytime a slider is activated
    def update(val):
        #change the value of the sun
        Sun[0] = SunSliderX.val
        Sun[1] = SunSliderY.val
        Sun[2] = SunSliderZ.val
        
        #change the value of the cloud
        CloudTime[0] = CloudSliderA.val
        CloudTime[1] = CloudSliderB.val
        CloudTime[2] = CloudSliderC.val
        CloudTime[3] = CloudSliderD.val
        
        #the title of the graph is the equation of the cloud
        axes.set_title(formatPlaneEquation(*CloudTime))
        
        # Update the position of the sun
        sun_scatter._offsets3d = ([Sun[0]], [Sun[1]], [Sun[2]])
        
        #Compute the  cloud's new plane
        X1, Y1, Z1 = plan(N, *LimitCloud, *CloudTime)
        cloud_points = [[X1[i], Y1[i], Z1[i]] for i in range(len(X1))]
        #update the position of the cloud on the window
        cloud_plot.set_data(X1, Y1)
        cloud_plot.set_3d_properties(Z1)
        
        # Compute the new cloud projection
        XShadow, YShadow, ZShadow = [], [], []
        for pt in cloud_points:
            x, y, z = projection(pt, pt, Sun, PlaneTime[0], PlaneTime[1], PlaneTime[2], [X[0], Y[0], Z[0]])
            XShadow.append(x)
            YShadow.append(y)
            ZShadow.append(z)
        # update the position of the shadow on the new plane
        shadow_plot.set_data(XShadow, YShadow)
        shadow_plot.set_3d_properties(ZShadow)
        
        fig.canvas.draw_idle()

    
    # Link the sliders with the update function
    SunSliderX.on_changed(update)
    SunSliderY.on_changed(update)
    SunSliderZ.on_changed(update)
    
    CloudSliderA.on_changed(update)
    CloudSliderB.on_changed(update)
    CloudSliderC.on_changed(update)
    CloudSliderD.on_changed(update)

    axes.set_xlabel("X")
    axes.set_ylabel("Y")
    axes.set_zlabel("Z")
    axes.legend(loc='upper left')
    plt.show()

#hepls display the equation of the plane   
def formatPlaneEquation(A:float, B:float, C:float, D:float):
    terms = []
    terms.append("PlanNuage = ")
    if A != 0: terms.append(f"{A:+.2f}·x")
    if B != 0: terms.append(f"{B:+.2f}·y")
    if C != 0: terms.append(f"{C:+.2f}·z")
    if D != 0: terms.append(f"{D:+.2f}")
    return ' '.join(terms) + " = 0"
    
#Parallel plane
# 10 000 points / plane
# Sun(4,4,7)
# MountainPlane ={(x,y,z) e R^3| 2x + 3y + 5z -6 = 0} 
# CloudPlane ={(x,y,z) e R^3| 2x + 3y + 5z -25 = 0} 
# Mountain : -5 <= x <= 5; -5 <= y <= 5
# Cloud : -1 <= x <= 1; -1 <= y <= 1
ombreParallele(10000, [4, 4, 7], [2, 3, 5, -6], [2, 3, 5, -25], [-5, 5, -5, 5], [-1, 1, -1, 1] )

#Non-Parallel plane
# 10 000 points / plane
# Sun(4,4,7)
# MountainPlane ={(x,y,z) e R^3| 2x + 3y + 5z -6 = 0} 
# CloudPlane ={(x,y,z) e R^3| x + y + 5z -25 = 0} 
# Mountain : -5 <= x <= 5; -5 <= y <= 5
# Cloud : -1 <= x <= 1; -1 <= y <= 1
ombreNonParallele(10000, [4, 4, 7], [2, 3, 5, -6], [1,1, 5, -25], [-5, 5, -5, 5], [-1, 1, -1, 1] )


# Moving Sun & moving cloud
# 10 000 points / plane
# Sun(4,4,7)
# MountainPlane ={(x,y,z) e R^3| 2x + 3y + 5z -6 = 0} 
# CloudPlane ={(x,y,z) e R^3| x + y + 5z -25 = 0} 
# Mountain : -5 <= x <= 5; -5 <= y <= 5
# Cloud : -1 <= x <= 1; -1 <= y <= 1
ombreMovingSun(10000, [4, 4, 7], [2, 3, 5, -6], [1,1, 5, -25], [-5, 5, -5, 5], [-1, 1, -1, 1] )