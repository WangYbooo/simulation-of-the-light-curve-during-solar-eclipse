from math import *
import numpy as np
def GetKey(item): 
    # order by reverse clockwise order (trigonometric convention)
    r=fmod(atan2(item[1],item[0]),2.*pi)
    return r
def GetKey2(item):  
    
    r=fmod(atan2(item[0][1],item[0][0]),2.*pi)  
    return r
def Hull(ens):  
    lens=len(ens)
    hu=[]
    for i in range(lens):
        if len(ens[i])==1:
            a=ens[i][0]
            hu.append(a)
        if len(ens[i])==2:
            a=ens[i][0]
            b=ens[i][1]
            hu.append(a)
            hu.append(b)
    return tuple(hu)

def rang(intersect):

    lu=len(intersect)
    intersect_r=[]
    for i in range(lu):
        seg=intersect[i]
        seg=OrderedPolygon(seg)
        intersect_r.append(seg)
    return tuple(intersect_r)
def rang2(intersect):

    lu=len(intersect)
    intersect_r=[]
    for i in range(lu):
        seg=intersect[i]
        seg=OrderedPolygon(seg)    
        intersect_r.append(seg)

    intersec_r=sorted(intersect_r,key=GetKey2)  
    moons=[]
    for i in range(lu):
        u=(intersect_r[i][0],intersect_r[i-1][1])  
        moons.append(u)
    return tuple(moons)
def sgn(x):   
    if x < 0. :
        sg=-1.
    else:
        sg=1.
    return sg
def BaryCenter(polygon): 
    xmean=0.
    ymean=0.
    l=len(polygon)
    for i in range(l):
        xmean=xmean+polygon[i][0]
        ymean=ymean+polygon[i][1]
    xmean=xmean/l
    ymean=ymean/l
    return (xmean,ymean)
def OrderedPolygon(polygon):    
    # put the vertices of a polygon in the order specified by getKey
    # polygon vertices = ((xA,yA),(xB,yB),(xC,yC),....)
    sor=sorted(polygon,key=GetKey) 
    return tuple(sor)
def lenline(line): 
    x1=line[0][0]
    y1=line[0][1]
    x2=line[1][0]
    y2=line[1][1]
    d=sqrt((x2-x1)**2+(y2-y1)**2)
    return d
def IntersectLineCircle(line,circle):   
    # line=((x1,y1),(x2,y2)), circle=((x0,y0),r)
    x0=circle[0][0]
    y0=circle[0][1]
    r=circle[1] 
    x1=line[0][0]-x0
    y1=line[0][1]-y0
    x2=line[1][0]-x0
    y2=line[1][1]-y0
    dx=x2-x1
    dy=y2-y1
    dr=sqrt(dx**2+dy**2)
    D=x1*y2-x2*y1
    delta=r**2*dr**2-D**2 
    sol=[]
    if delta >= 0. :

       xA=(D*dy+sgn(dy)*dx*sqrt(delta))/dr**2
       xB=(D*dy-sgn(dy)*dx*sqrt(delta))/dr**2
       yA=(-D*dx+abs(dy)*sqrt(delta))/dr**2
       yB=(-D*dx-abs(dy)*sqrt(delta))/dr**2
       sol.append((xA+x0,yA+y0))
       sol.append((xB+x0,yB+y0))    
    return tuple(sol)

def IntersectInsideSegment(seg,line): 
    sol=[]
    if len(seg) == 0: return tuple(sol)
    xA=seg[0][0]
    yA=seg[0][1]
    xB=seg[1][0]
    yB=seg[1][1]
    x1=line[0][0]
    y1=line[0][1]
    x2=line[1][0]
    y2=line[1][1]
    if x1 != x2 :   
        if (xA >= x1 and xA <= x2):
            sol.append((xA,yA))
        if xA >= x2 and xA <= x1:
            sol.append((xA,yA))
        if (xB >= x1 and xB <= x2):
            sol.append((xB,yB))
        if xB >= x2 and xB <= x1:
            sol.append((xB,yB))
    if x1 == x2:
        if yA >= y1 and yA <= y2:
            sol.append((xA,yA))
        if yA >= y2 and yA <= y1:
            sol.append((xA,yA))
        if yB >= y1 and yB <= y2:
            sol.append((xB,yB))
        if yB >= y2 and yB <= y1:
            sol.append((xB,yB))
    return tuple(sol)


def IntersectInsideSegment_2(seg,line): 
    sol = []
    if len(seg) == 0:
        return tuple(sol)
    xA=seg[0][0]
    yA=seg[0][1]
    xB=seg[1][0]
    yB=seg[1][1]
    x1=line[0][0]
    y1=line[0][1]
    x2=line[1][0]
    y2=line[1][1]
    vector1 = ((xA-x1),(yA-y1))
    vector2 = ((xA-x2),(yA-y2))
    vector3 = ((xB-x1),(yB-y1))
    vector4 = ((xB-x2),(yB-y2))
    if np.dot(vector1,vector2) <= 0.:   
        sol.append((xA,yA))
    if np.dot(vector3,vector4) <= 0.:
        sol.append((xB,yB))
    return tuple(sol)

def PointInCircle(point,circle):  
    #circle=((x0,y0),r)
    # pts on the circle are inside the circle
    x=point[0]
    y=point[1]
    x0=circle[0][0]
    y0=circle[0][1]
    r=circle[1]
    if sqrt((x-x0)**2+(y-y0)**2) < r : return True
    return False

def HalfLens(a,r):
    theta=2.*asin(a/2./r) 
    R=r*cos(theta/2.) 
    area=0.5*(r**2*theta-a*R)
    return area
def TriangleSurface(triangle): 
    # triangle vertices = ((xA,yA),(xB,yB),(xC,yC))
    if len(triangle) != 3 :
        raise BaseException
    xA=triangle[0][0]
    yA=triangle[0][1]
    xB=triangle[1][0]
    yB=triangle[1][1]
    xC=triangle[2][0]
    yC=triangle[2][1]
    a=sqrt((xA-xB)**2+(yA-yB)**2)
    b=sqrt((xB-xC)**2+(yB-yC)**2)
    c=sqrt((xC-xA)**2+(yC-yA)**2)
    s=0.5*(a+b+c)
    area=sqrt(s*(s-a)*(s-b)*(s-c))
    return area
def PolygonSurface(polygon): 
    polygon=list(polygon)
    n=len(polygon)
    # print("1=",polygon)
    # print("np=",n)
    if  n == 0 : 
        raise BaseException
    # print("0=",polygon[0])
    polygon.append(polygon[0]) 
    # print("2=",polygon)
    s=0.
    for i in range(n):
        xA=polygon[i][0]
        yA=polygon[i][1]
        xB=polygon[i+1][0]
        yB=polygon[i+1][1]
        #print(xA,xB,yA,yB)
        s=s+(xA*yB-xB*yA)
    return 0.5*s
def PointInTriangle(point,triangle):
    # triangle vertices = ((xA,yA),(xB,yB),(xC,yC))
    xA=triangle[0][0]
    yA=triangle[0][1]
    xB=triangle[1][0]-xA
    yB=triangle[1][1]-yA
    xC=triangle[2][0]-xA
    yC=triangle[2][1]-yA
    x=point[0]-xA
    y=point[1]-yA
    a=(x*yC-xC*y)/(xB*yC-xC*yB)
    b=-(x*yB-xB*y)/(xB*yC-xC*yB)
    if a < 0 : return False
    if b < 0 : return False
    if a+b > 1 : return False
    return True

def VerticesInsideCircle(triangle,circle): #在圆内的点形成列表

    inside=[]
    if PointInCircle(triangle[0],circle): inside.append(triangle[0])
    if PointInCircle(triangle[1],circle): inside.append(triangle[1])
    if PointInCircle(triangle[2],circle): inside.append(triangle[2])
    return inside
def VerticesOutsideCircle(triangle,circle):  #不在圆内的点形成列表
    outside=[]
    if not PointInCircle(triangle[0],circle): outside.append(triangle[0])
    if not PointInCircle(triangle[1],circle): outside.append(triangle[1])
    if not PointInCircle(triangle[2],circle): outside.append(triangle[2])
    return outside
def EdgesInterceptCircle(triangle,circle):  #三角形三边与圆的交点，且交点在三角形三边上；用点的位置判断完全在圆内的边；返回的结果将每条边上的点放在一起e[0],e[1],e[2]
    ens=[]
    line1=(triangle[0],triangle[1])
    intersect1=IntersectLineCircle(line1,circle)  #直线与圆的交点

    intersect1=IntersectInsideSegment_2(intersect1,line1) #三角形一边所在的直线与圆的交点如果在边内则将其组成一个元组

    if len(intersect1) == 2 : ens.append(intersect1)
    line2=(triangle[0],triangle[2])
    intersect2=IntersectLineCircle(line2,circle)
    intersect2=IntersectInsideSegment_2(intersect2,line2)

    if len(intersect2) == 2 : ens.append(intersect2)
    line3=(triangle[1],triangle[2])
    intersect3=IntersectLineCircle(line3,circle)
    intersect3=IntersectInsideSegment_2(intersect3,line3)

    if len(intersect3) == 2 : ens.append(intersect3)
    return tuple(ens)
def EdgesInterceptCircle_single(triangle,circle): #三角形三边与圆的交点 用点的位置判断完全在圆内的边
    ens=[]
    line1=(triangle[0],triangle[1])

    intersect1=IntersectLineCircle(line1,circle)  #直线与圆的交点

    intersect1=IntersectInsideSegment_2(intersect1,line1) #三角形一边所在的直线与圆的交点如果在边内则将其组成一个元组

    if len(intersect1) >= 1 : ens.append(intersect1)
    line2=(triangle[0],triangle[2])

    intersect2=IntersectLineCircle(line2,circle)

    intersect2=IntersectInsideSegment_2(intersect2,line2)

    if len(intersect2) >= 1 : ens.append(intersect2)
    line3=(triangle[1],triangle[2])

    intersect3=IntersectLineCircle(line3,circle)

    intersect3=IntersectInsideSegment_2(intersect3,line3)

    if len(intersect3) >= 1 : ens.append(intersect3)
    return tuple(ens)
def Closepoint(a,b,ens):  
        if len(ens)==2:
            if len(ens[0])==len(ens[1])==2:   #m=0,n=2 
                x11=ens[0][0][0]
                y11=ens[0][0][1]
                x12=ens[0][1][0]
                y12=ens[0][1][1]
                x21=ens[1][0][0]
                y21=ens[1][0][1]
                x22=ens[1][1][0]
                y22=ens[1][1][1]
                l1=(x11-x21)**2+(y11-y21)**2
                l2=(x11-x22)**2+(y11-y22)**2
                if l1 < l2:
                    return ((x11,y11),(x21,y21)),((x12,y12),(x22,y22))
                else:
                    return ((x11,y11),(x22,y22)),((x12,y12),(x21,y21))
            if len(ens[0])==1:  
                x11=ens[0][0][0]
                y11=ens[0][0][1]
                x21=ens[1][0][0]
                y21=ens[1][0][1]
                x22=ens[1][1][0]
                y22=ens[1][1][1]
                l1=(x11-x21)**2+(y11-y21)**2
                l2=(x11-x22)**2+(y11-y22)**2
                if l1<l2:
                    return ((x11,y11),(x21,y21)),((0,0),(0,0))  
                else:
                    return ((x11,y11),(x22,y22)),((0,0),(0,0))
            if len(ens[1])==1:  
                x11=ens[0][0][0]
                y11=ens[0][0][1]
                x12=ens[0][1][0]
                y12=ens[0][1][1]
                x21=ens[1][0][0]
                y21=ens[1][0][1]
                l1=(x11-x21)**2+(y11-y21)**2
                l2=(x21-x12)**2+(y21-y12)**2
                if l1<l2:
                    return ((x11,y11),(x21,y21)),((0,0),(0,0))
                else:
                    return ((x21,y21),(x21,y21)),((0,0),(0,0))
            
        if len(ens)==3:
            if len(ens[0])==1 or len(ens[1])==1 or len(ens[2])==1:  #m=1,n=1
                if len(ens[0])==2:  
                    x11=ens[0][0][0]
                    y11=ens[0][0][1]
                    x12=ens[0][1][0]
                    y12=ens[0][1][1]
                    x2=ens[1][0][0]
                    y2=ens[1][0][1]
                    x3=ens[2][0][0]
                    y3=ens[2][0][1]
                if len(ens[1])==2:
                    x11=ens[1][0][0]
                    y11=ens[1][0][1]
                    x12=ens[1][1][0]
                    y12=ens[1][1][1]
                    x2=ens[0][0][0]
                    y2=ens[0][0][1]
                    x3=ens[2][0][0]
                    y3=ens[2][0][1]
                if len(ens[2])==2:
                    x11=ens[2][0][0]
                    y11=ens[2][0][1]
                    x12=ens[2][1][0]
                    y12=ens[2][1][1]
                    x2=ens[1][0][0]
                    y2=ens[1][0][1]
                    x3=ens[0][0][0]
                    y3=ens[0][0][1]
                l1=(x2-x11)**2+(y2-y11)**2
                l2=(x2-x12)**2+(y2-y12)**2
                if l1<l2:
                    return ((x11,y11),(x2,y2)),((x12,y12),(x3,y3))
                else:
                    return ((x11,y11),(x3,y3)),((x12,y12),(x2,y2))
            if len(ens[0])==len(ens[1])==len(ens[2])==2:
                x11=ens[a][0][0]
                y11=ens[a][0][1]
                x12=ens[a][1][0]
                y12=ens[a][1][1]
                x21=ens[b][0][0]
                y21=ens[b][0][1]
                x22=ens[b][1][0]
                y22=ens[b][1][1]
                l1=(x11-x21)**2+(y11-y21)**2
                l2=(x11-x22)**2+(y11-y22)**2
                if l1 < l2:
                    return ((x11,y11),(x21,y21))
                else:
                    return ((x11,y11),(x22,y22))
            


def MidSegment(ens): 
    if len(ens) != 2 : raise BaseException
    xA=ens[0][0]
    yA=ens[0][1]
    xB=ens[1][0]
    yB=ens[1][1]
    # coordinates of the mid-point of the chord
    xM=(xB+xA)/2.
    yM=(yB+yA)/2.
    return(xM,yM)
def CrossChord(ens,circle):  
    # coordinates of the mid-point of the chord
    (xM,yM)=MidSegment(ens)
    r=circle[1]
    xA=ens[0][0]
    yA=ens[0][1]
    xB=ens[1][0]
    yB=ens[1][1]
    xD=xB-xA
    yD=yB-yA
    norm=sqrt(xD**2+yD**2)
    xD=xD/norm
    yD=yD/norm
    # print(xD,yD)
    # normal to the chord, pointing right
    xNright=yD
    yNright=-xD

    # compute two points perpendicular to the chord 
    x1=xM+r*xNright
    y1=yM+r*yNright
    x2=xM-r*xNright
    y2=yM-r*yNright
    line=((x1,y1),(x2,y2))  

    # compute the intersection of the line with the circle
    intersect=IntersectLineCircle(line,circle)
    d=lenline(intersect)
    return intersect




def AreaInterceptCircle(ens,r):    
    areas=[]
    for intersect in ens:
        area=0.
        if len(intersect) == 2 :
           pt1=intersect[0]
           pt2=intersect[1]
           a=sqrt((pt1[0]-pt2[0])**2+(pt1[1]-pt2[1])**2)
           area=HalfLens(a,r)
           areas.append(area)
    return tuple(areas)  
def InterceptTriangleCircle(triangle,circle): 

    # common area between a triangle and a circle (any case)
    (x0,y0)=BaryCenter(triangle)
    # all is centered on triangle barycenter

    
    triangle=((triangle[0][0]-x0,triangle[0][1]-y0),  
              (triangle[1][0]-x0,triangle[1][1]-y0),
              (triangle[2][0]-x0,triangle[2][1]-y0))
    circle=((circle[0][0]-x0,circle[0][1]-y0),circle[1])
    
    # print('triangle:',triangle)
    # print('circle:',circle)

    centercircle=circle[0]
    r=circle[1]

    verticesinside=VerticesInsideCircle(triangle,circle)
    VerticesOutside=VerticesOutsideCircle(triangle,circle)

    # m number of vertices inside or on the circle
    m=len(verticesinside)

    edgesintercept=EdgesInterceptCircle(triangle,circle)  #
    # n number of edges inside or on the circle 
    n=len(edgesintercept)
    if m == 3 :                            
        area=TriangleSurface(triangle)
        return area                                  
    # m = 0 all the vertices are not in the circle
    if m == 0:  
       # no intersection with the edges 
       if n == 0 :                         
          if PointInTriangle(centercircle,triangle): 
             # circle completely in triangle
             area=pi*r**2
             return area
          else :                                     
             # circle completely outside triangle
             area=0.
             return area
       elif n == 1:                               
             edgesnew=Hull(edgesintercept)
             line=[]
             # only one edge of the triangle intercepts the circle
            #  print(edgesintercept)
             cross=CrossChord(edgesintercept[0],circle) 
             a=lenline(edgesintercept[0])
             midchord=MidSegment(edgesintercept[0])
             point=cross[0]
             line.append(point)
             line.append(midchord)
             line=tuple(line) 
             d=lenline(line) 
             if PointInTriangle(point,triangle):  
                if d > r :       
                 area=pi*r**2-HalfLens(a,r)  
                 return area
                else:
                 area=HalfLens(a,r)    
                 return area
             line=[]
             point=cross[1]
             line.append(point)
             line.append(midchord)
             line=tuple(line)
             d=lenline(line)
             if PointInTriangle(point,triangle):    
                if d > r :
                 area=pi*r**2-HalfLens(a,r)
                 return area
                else:
                 area=HalfLens(a,r)
                 return area
       elif n == 2 :                       
             # two and only two edges of the triangle intercepts the circle
            #  print(edgesintercept)
             secor1,sector2=Closepoint(0,1,edgesintercept)  #两个截出弧形面积的点
             a1=lenline(secor1)
             a2=lenline(sector2)
             area_half1=HalfLens(a1,r)
             area_half2=HalfLens(a2,r)
             edgesnew=Hull(edgesintercept)
             convexhull=OrderedPolygon(edgesnew)
             area=PolygonSurface(convexhull)+area_half1+area_half2   #同样只用了PS这个函数，没考虑弧线
             return area   
       elif n == 3:                                 #6. 0顶点，3相交边(5)    有问题
             # the three vertices of the triangle are outside the circle
             # but the edges of the triangle are intercepting the circle
             # 3 possible cases, each one for each edge
             areas = AreaInterceptCircle(edgesintercept,r)
             area=pi*r**2-areas[0]-areas[1]-areas[2]
             return area
       else : raise BaseException 
    if m == 2 :                                  #7.  2顶点，3边（8） 实现？
       # two vertices of the triangle are inside the circle
       edgesintercept_single=EdgesInterceptCircle_single(triangle,circle)
       edgesintercept_single=Hull(edgesintercept_single)
       if len(edgesintercept_single) == 1:
           edgesintercept_single=list(edgesintercept_single)
           x1 = edgesintercept_single[0][0]
           y1 = edgesintercept_single[0][1]
           vector1 = ((x1-triangle[0][0]),(y1-triangle[0][1]))
           vector2 = ((x1-triangle[1][0]),(y1-triangle[1][1]))
           vector3 = ((x1-triangle[2][0]),(y1-triangle[2][1]))
           if abs(1-abs(np.dot(vector1,vector2)/(np.linalg.norm(vector1)*np.linalg.norm(vector2))))<1e-10:
               edgesintercept_single.append((triangle[2][0],triangle[2][1]))
           if abs(1-abs(np.dot(vector1,vector3)/(np.linalg.norm(vector1)*np.linalg.norm(vector3))))<1e-10:
               edgesintercept_single.append((triangle[1][0],triangle[1][1]))
           if abs(1-abs(np.dot(vector2,vector3)/(np.linalg.norm(vector2)*np.linalg.norm(vector3))))<1e-10:
               edgesintercept_single.append((triangle[0][0],triangle[0][1]))
       a=lenline(edgesintercept_single)
       area_half=HalfLens(a,r)       
       convexhull=list(verticesinside)+list(edgesintercept_single)
       convexhull=OrderedPolygon(convexhull)
       area=PolygonSurface(convexhull) + area_half 

       return area
    if m == 1 :                                  
       # one vertice of the triangle is inside the circle
       if n == 0 :                             
             edgesintercept_single=EdgesInterceptCircle_single(triangle,circle)

             edgesintercept_single=Hull(edgesintercept_single)
             a=lenline(edgesintercept_single)
             area_half=HalfLens(a,r)
             convexhull=list(verticesinside)+list(edgesintercept_single)
             convexhull=OrderedPolygon(convexhull)
             area=PolygonSurface(convexhull)+area_half
             return area
       if n == 1 :                               #9.   1顶点，3边（7）  实现？
             edgesintercept_single=EdgesInterceptCircle_single(triangle,circle)
             sector=Closepoint(0,1,edgesintercept_single)
             a1=lenline(sector[0])
             a2=lenline(sector[1])
             area_half1=HalfLens(a1,r)
             area_half2=HalfLens(a2,r)
             edgesintercept_single=Hull(edgesintercept_single)
             convexhull=list(verticesinside)+list(edgesintercept_single)
             convexhull=OrderedPolygon(convexhull)
             area = PolygonSurface(convexhull)+area_half1+area_half2
             return area
