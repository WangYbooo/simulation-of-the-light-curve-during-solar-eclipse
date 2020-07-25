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

d2r = np.pi / 180.0

spice.furnsh("path/lakeland/eclipse_ground.tm")


time = Time(datetime(2012,11,13,20,56,00,0),scale='utc')


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

    axis_z = np.cross(axis_x,axis_y)    
    
    return axis_x,axis_y,axis_z


def rotation_matrix(et):
    axis_x,axis_y,axis_z = selframe_axis(et)
    coord_topo = np.array([axis_x,axis_y,axis_z]).T  
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
print(i_obs)
print(j_obs)
vol=triple_prod(los_obs,i_obs,j_obs)
# volume is positive is the 3-vectors define a direct frame,

print("volume=",vol)

R_me2topo=spice.pxform('Moon_ME','GROUNDSTATION_TOPO',et)

North_Pole_moon=np.dot(R_me2topo,(0,0,1))  #ME系的北极方向
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

with open('path/LRO_dem_256.txt','r') as file:
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


point_limb = sorted(point_reorder,key=(lambda x:x[0]),reverse=False) 


h_limb=[]

for i in range(len(point_limb)):
    h_limb.append(point_limb[i][1])



r_limb=[]

for j in range(len(h_limb)):
    r_limb.append(float(1737400.0+(0.5 * float(h_limb[j]))))
    

r_limb.append(r_limb[0])   #首尾相连
theta=list(theta)
theta.append(theta[0])

print('r_limb',len(r_limb))

x_limb=[]
x_mean=[]
y_limb=[]
y_mean=[]
for i in range(len(r_limb)):
    x_limb.append(r_limb[i]*np.cos(deg2rad(theta[i])))
    x_mean.append(1737400.0*np.sin(deg2rad(theta[i])))
    y_limb.append(r_limb[i]*np.sin(deg2rad(theta[i])))
    y_mean.append(1737400.0*np.cos(deg2rad(theta[i])))

x_limb=tuple(x_limb)
x_mean=tuple(x_mean)
y_limb=tuple(y_limb)
y_mean=tuple(y_mean)


plt.figure(figsize=(10, 10))
ax = plt.subplot(111)
ax.plot(x_mean,y_mean,'--',linewidth=2,color='black',alpha=0.5)
ax.plot(x_limb,y_limb,linewidth=1,color='black',alpha=2)

plt.axis('scaled')

ax=plt.gca()


x_major_locator=MultipleLocator(500000)   
y_major_locator=MultipleLocator(500000)
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)

plt.xlim(-2000000,2000000)
plt.ylim(-2000000,2000000)



ax.spines['top'].set_color('none')   
ax.spines['right'].set_color('none')

ax.xaxis.set_ticks_position('bottom')   
ax.spines['bottom'].set_position(('data', 0))    
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))



plt.xlabel('$meters$',fontsize=10,position='1900000,0')
plt.ylabel('$meters$',rotation=90,horizontalalignment='center',position='-1,1900000')




plt.show()
