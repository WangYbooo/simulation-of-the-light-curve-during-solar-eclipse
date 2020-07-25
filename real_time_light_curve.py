import math
import scipy as sci
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
from math import *
from matplotlib.pyplot import MultipleLocator
from astropy import units as u
from astropy.time import Time
import time as tm
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_moon
from astropy import constants
from datetime import datetime,tzinfo,timedelta
from matplotlib.patches import Ellipse, Circle, Polygon


cellsize = 1/256     
N_limb = int(360/cellsize)  
theta = np.arange(0,360,cellsize)

time_duration = Time(datetime(2012,11,13,20,36,00,0),scale='utc') + (np.linspace(0,3,4))*u.second

spice.furnsh("path/lakeland/eclipse_ground.tm")

file = open('path/light_limb_real_5min.txt','w')
file_perfect = open('path/light_perfect_real_5min.txt','w')

def momentlight(circle):  
    xA=0.
    yA=0.
    area_sum=0.
    for i in range(len(r_limb)-1): 
        xB=r_limb[i]*cos(deg2rad(theta[i]))
        yB=r_limb[i]*sin(deg2rad(theta[i]))
        xC=r_limb[i+1]*cos(deg2rad(theta[i+1]))
        yC=r_limb[i+1]*sin(deg2rad(theta[i+1]))

        triangle=((xA,yA),(xB,yB),(xC,yC))
        
        area=its.InterceptTriangleCircle(triangle,circle)

        area_sum = area_sum+area
    return area_sum

def mask(r_s,R_m,d):
    # r : sun, R : moon, d : distance between centers
    if R_m >= r_s:
       if r_s == 0. : return 0.
       if d >= r_s+R_m:
           light=pi*r_s**2
       elif d <= R_m-r_s:
           light=0.
       else:
           a=2.*acos((r_s**2+d**2-R_m**2)/2./r_s/d)
           A=2.*acos((R_m**2+d**2-r_s**2)/2./R_m/d)
           light=pi*r_s**2-0.5*r_s**2*(a-sin(a))-0.5*R_m**2*(A-sin(A))
       return light
    else:
       if R_m == 0. : return 0.

def project_unit(v1,n):

    dis_v = np.dot(v1,n)
    v_n = n*dis_v
    v_project = v1 - v_n
    v_project = v_project/np.linalg.norm(v_project)
    return v_project


def selframe_axis(et):

    moon_obsever,lt = spice.spkpos('moon',et,'GROUNDSTATION_TOPO','LT+S','399123')
    sun_observer,lt = spice.spkpos('sun',et,'GROUNDSTATION_TOPO','LT+S','399123')
    dis_moon = np.linalg.norm(moon_obsever)
    dis_sun = np.linalg.norm(sun_observer)

    axis_x = moon_obsever/dis_moon

    ele_ref_obs=90.
    azi_ref_obs=0.

    ele_ref_obs=deg2rad(ele_ref_obs)
    azi_ref_obs=deg2rad(azi_ref_obs)
    ref_obs=pol2car((ele_ref_obs,azi_ref_obs))

    axis_y=cross_prod(axis_x,ref_obs) # i_obs is NOT normalized
    axis_y=normalize(axis_y)

    axis_z = np.cross(axis_x,axis_y)    #y与z叉乘才满足右手定则
    
    return axis_x,axis_y,axis_z


def rotation_matrix(et):

    axis_x,axis_y,axis_z = selframe_axis(et)
    coord_topo = np.array([axis_x,axis_y,axis_z]).T  #将三个坐标写成矩阵形式
    R_topo2self = np.linalg.inv(coord_topo)
    R_topo2j2000 = spice.pxform('GROUNDSTATION_TOPO','J2000',et)
    R_j20002me = spice.pxform('j2000','moon_me',et)
    
    return R_topo2self,R_topo2j2000,R_j20002me


def deg2rad(x):
    y=x/180.*pi
    return y
def rad2deg(x):
    y=x/pi*180.
    return y
def dot_prod(x,y):
    p=x[0]*y[0]+x[1]*y[1]+x[2]*y[2]
    return p
def norm(x):
    n=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
    return n
def normalize(x):
    n=norm(x)
    return (x[0]/n,x[1]/n,x[2]/n)
def cross_prod(x,y):
    p=[0.,0.,0.]
    p[0]=x[1]*y[2]-x[2]*y[1]
    p[1]=x[2]*y[0]-x[0]*y[2]
    p[2]=x[0]*y[1]-x[1]*y[0]
    return p
def triple_prod(x,y,z):
    p=dot_prod(x,cross_prod(y,z))
    return p
def pol2car(pol):
    lat=pol[0]
    lon=pol[1]
    car=[0.,0.,0.]
    car[0]=cos(lat)*cos(lon)
    car[1]=cos(lat)*sin(lon)
    car[2]=sin(lat)
    return car
def car2pol(car):
    lat=asin(car[2])
    lon=atan2(car[1],car[0])
    return (lat,lon)

def datetime_tostr(t):   
    return(t.strftime('%Y-%m-%d %H:%M:%S.%f'))

def datetime_2str_hms(t):
    return(t.strftime('%H:%M:%S'))

def datetime_tosecond(t):
    t = datetime_tostr(t)
    struct_t = tm.strptime(t,'%Y-%m-%d %H:%M:%S.%f')
    t_t=tm.mktime(struct_t)
    return(t_t)




value_light=[]
value_light_perfect=[]

for time in time_duration:

    t = datetime_tostr(time)
    t_pic = datetime_2str_hms(time)
    et = spice.utc2et(t)

    moon_obsever,lt = spice.spkpos('moon',et,'GROUNDSTATION_TOPO','LT+S','399123')
    sun_observer,lt = spice.spkpos('sun',et,'GROUNDSTATION_TOPO','LT+S','399123')
    dis_moon = np.linalg.norm(moon_obsever)
    dis_sun = np.linalg.norm(sun_observer)

    r_moon=1737.4

    D=dis_moon*(1-(r_moon/dis_moon)**2)
    r=r_moon*sqrt(1-(r_moon/dis_moon)**2)
    om=moon_obsever
    print(D)


    half_angle_view_moon=atan(r_moon/dis_moon)

    los_obs=normalize(moon_obsever)
    ele_ref_obs=90.
    azi_ref_obs=0.
    print("ele_ref_obs=",ele_ref_obs,", azi_ref_obs=",azi_ref_obs)
    ele_ref_obs=deg2rad(ele_ref_obs)
    azi_ref_obs=deg2rad(azi_ref_obs)
    ref_obs=pol2car((ele_ref_obs,azi_ref_obs))

    i_obs=cross_prod(los_obs,ref_obs) # i_obs is NOT normalized
    i_obs=normalize(i_obs)
    j_obs=cross_prod(los_obs,i_obs)

    vol=triple_prod(los_obs,i_obs,j_obs)

    R_me2topo=spice.pxform('Moon_ME','GROUNDSTATION_TOPO',et)

    North_Pole_moon=np.dot(R_me2topo,(0,0,1)) 
    M_moon=np.dot(R_me2topo,(1,0,0))
    K_moon=np.dot(R_me2topo,(0,1,0))

    delta_N_limb=360./N_limb
    plot_lat_limb=[]
    plot_lon_limb=[]
    x_p=[]
    y_p=[]
    z_p=[]
    lon_limb=[]
    lat_limb=[]
    for i in range(N_limb):
        w=[0.,0.,0.]
        theta_limb=float(i-1)*delta_N_limb
        theta_limb=deg2rad(theta_limb)
        cos_theta=cos(theta_limb)
        sin_theta=sin(theta_limb)  
        w[0]=(D*los_obs[0]+r*(cos_theta*i_obs[0]+sin_theta*j_obs[0]))-om[0]
        w[1]=(D*los_obs[1]+r*(cos_theta*i_obs[1]+sin_theta*j_obs[1]))-om[1]
        w[2]=(D*los_obs[2]+r*(cos_theta*i_obs[2]+sin_theta*j_obs[2]))-om[2]
        w=normalize(w)
        (lat_w,lon_w)=car2pol(w)
        x=dot_prod(w,M_moon) 
        y=dot_prod(w,K_moon) 
        z=dot_prod(w,North_Pole_moon) 

        x_p.append(x)
        y_p.append(y)
        z_p.append(z)
        (lat,lon)=car2pol((x,y,z)) # polar coordinates of limb on moon crust frame

        lat_limb.append(rad2deg(lat))
        lon_limb.append(rad2deg(lon))


    p_n = np.arange(N_limb) 


    dem_line_1=[]
    dem_line_2=[]
    dem_row_1=[]
    dem_row_2=[]

    for i in range(len(lat_limb)):

        line_1 = math.floor(((90.0-cellsize/2) - float(lat_limb[i]))/cellsize) + 6   
        line_2 = line_1 + 1
        
        row_1 = abs(math.floor((float(lon_limb[i])-(-180+cellsize/2))/cellsize))  
        row_2 = row_1 + 1

        dem_line_1.append(line_1)
        dem_row_1.append(row_2)
        dem_line_2.append(line_2)
        dem_row_2.append(row_2)





    def list_2d(list1,list2,list3):
        first = []                      
        for i in range(len(list1)):
            first.append([int(list1[i])])    
            first[i].append(int(list2[i]))
            first[i].append(int(list3[i]))
        return(first)




    point_list = list_2d(dem_line_1,dem_row_1,p_n)  

    point_order_lat = sorted(point_list,key=(lambda x:x[0]))   


    lat_order=[]
    lon_order=[]
    h_order=[]

    for i in range(len(point_order_lat)):
        lat_order.append(point_order_lat[i][0])
        lon_order.append(point_order_lat[i][1])


    latdict_dem={}

    for i in range(len(lat_order)):   
        if lat_order[i] not in latdict_dem:
            latdict_dem[lat_order[i]]=[]
            latdict_dem[lat_order[i]].append(lon_order[i])
        else:
            latdict_dem[lat_order[i]].append(lon_order[i])

    n = lat_order[0]
    j = 0

    with open('path_DEM','r') as file:
        for i in range(46086):
            rows = file.readline()

            if i == n:   
                for l in range(len(latdict_dem[lat_order[0+j]])):
                    m = latdict_dem[lat_order[0+j]][l]
                    h_order.append(rows.split(' ')[m])
                j = j + len(latdict_dem[n])  
                if j < len(lat_order):        
                    n = lat_order[j]



    point_reorder = []

    for i in range(len(h_order)):
        point_reorder.append([point_order_lat[i][2]])
        point_reorder[i].append(h_order[i])


    point_limb = sorted(point_reorder,key=(lambda x:x[0]),reverse=False)  #使用sorted排序


    h_limb=[]

    for i in range(len(point_limb)):
        h_limb.append(point_limb[i][1])



    r_limb=[]

    for j in range(len(h_limb)):
        r_limb.append(float(1737400.0+(0.5 * float(h_limb[j]))))


    r_limb.append(r_limb[0])  
    theta=list(theta)
    theta.append(theta[0])

    R_topo2self,R_topo2j2000,R_j20002me = rotation_matrix(et)
    


    angle_separation = np.arccos(np.dot(moon_obsever,sun_observer)/(dis_moon*dis_sun))
    a = np.dot(moon_obsever,sun_observer)/(dis_moon*dis_sun)
    sun_apparent_topo = sun_observer/np.linalg.norm(sun_observer)*(dis_moon/a)  #将太阳位置归算到与observe——moon垂直的平面
    
    d=dis_moon*np.tan(angle_separation) *1000     #太阳在XY平面上


    R_sun = constants.R_sun.to(u.km)
    R_sun = R_sun.value
    R_moon = 1737.4
    R_moon_adjuse=1736.85

    sunsize = R_sun * (dis_moon/a)/dis_sun*1000
    moonsize = R_moon*1000
    moonadjusesize=R_moon_adjuse*1000

    d_size = moonsize-sunsize

    sunpos_eclipse = np.dot(R_topo2self,sun_apparent_topo-moon_obsever)    #坐标转换时，需要先平移到M点的位置再乘以旋转矩阵

    sunpos_moon_XY = (project(sunpos_eclipse,np.array([1,0,0]))[1]*1000,project(sunpos_eclipse,np.array([1,0,0]))[2]*1000)

    circle = (sunpos_moon_XY,sunsize)

    print(t_pic)

    light = np.pi*sunsize*sunsize-momentlight(circle)
    light_perfect = mask(sunsize,moonsize,d)

    t_light.append(t_pic)
    value_light.append(light)
    value_light_perfect.append(light_perfect)

    file.write(str(light))
    file.write('\n')
    file_perfect.write(str(light_perfect))
    file_perfect.write('\n')