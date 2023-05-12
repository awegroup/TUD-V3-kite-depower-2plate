'''
-*- coding: utf-8 -*-
Written by J.A.W. Poland on 12/05/2023
This script computes the width of a triangular two-plate representation of the V3 kite
It plots the following results:
    - Tetrahedron algorithm for the CAD shape and the presimulated shape
    - Trilateration algorithm for the CAD shape and the presimulated shape
    - Photogrammetry measurements at the powered (u_p = 1) and depowered (u_p = 0) position
'''
#%%

### Import modules
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

### Define design shape (CAD)
P0 = np.array([    0.        ,     0.        ,     0.        ])
P1 = np.array([ 1545.03267768, -4130.03979687,  7261.36495027])
P2 = np.array([    0.        ,     0.        , 11000.        ])
P3 = np.array([ 1545.03267768,  4130.03979687,  7261.36495027])
P4 = np.array([ 2200.        ,     0.        , 11000.        ])
points_CAD = [P0,P1,P2,P3,P4]

# Defining line lengths
a       = np.linalg.norm(P2-P3)
b       = np.linalg.norm(P0-P3)
cref    = np.linalg.norm(P2-P4)
d       = np.linalg.norm(P0-P2)
e       = np.linalg.norm(P4-P3)
lines_CAD = [a,b,cref,d,e]

print('--- design CAD shape, line lengths ---')
print('a    : ',round(a   ,1), 'mm')  
print('b    : ',round(b   ,1), 'mm')  
print('cref : ',round(cref,1), 'mm') 
print('d    : ',round(d   ,1), 'mm')    
print('e    : ',round(e   ,1), 'mm') 
print(' ')


### Define presimulated shape
P0 = np.array([    0.        ,     0.        ,     0.        ])
P1 = np.array([  896.63860741, -4152.35709084,  7427.97340237])
P2 = np.array([    0.        ,     0.        , 11000.        ])
P3 = np.array([  896.63860741,  4152.35709084,  7427.97340237])
P4 =  np.array([ 2199.87097737,     0.        , 11021.07129607])
points_presimulated = [P0,P1,P2,P3,P4]

# Defining line lengths
a       = np.linalg.norm(P2-P3)
b       = np.linalg.norm(P0-P3)
cref    = np.linalg.norm(P2-P4)
d       = np.linalg.norm(P0-P2)
e       = np.linalg.norm(P4-P3)
lines_presimulated = [a,b,cref,d,e]

print('--- presimulated shape, line lengths ---')
print('a    : ',round(a   ,1), 'mm')  
print('b    : ',round(b   ,1), 'mm')  
print('cref : ',round(cref,1), 'mm') 
print('d    : ',round(d   ,1), 'mm')    
print('e    : ',round(e   ,1), 'mm') 
print(' ')

### Define depower-tape extenstion relation to the power-setting u_p
def up_to_ld(u_p,delta_ld):
    '''
    input:  u_p         = [0-1]   (power-setting)
            delta_ld    = 8 or 13 [%] (depower-tape extension, in percentage)
    output: depower-tape extension [mm]
    '''
    if delta_ld == 8:
        depower_tape_max = 1.482
    elif delta_ld == 13:
        depower_tape_max = 1.722
    else:
        raise Exception('delta_ld wrong value, should be: 8 or 13')

    depower_tape_0 = 1.098
    return 1e3*0.5*(1-u_p)*(depower_tape_max-depower_tape_0)


### Define Tetrahedron algorithm
def tetrahedron_algorithm(u_p,delta_ld,points,variables):
    '''
    --- Tetrahedron algorithm ---
    input:  u_p         = [0-1]   (power-setting)
            delta_ld    = 8 or 13 [%] (depower-tape extension)
            points      = [P0,P1,P2,P3,P4] (P0 is the origin)
            variables   = [a,b,cref,d,e] (a,b,cref,d,e are the line lengths)
    output: w           = width [mm]    (width of the kite)
    '''
    ### defining the line lengths
    a,b,cref,d,e = variables[0],variables[1],variables[2],variables[3],variables[4]
    
    ### finding the new length by adding the depower-tape extension
    lim  = np.linalg.norm(points[0]-points[4]) +up_to_ld(u_p,delta_ld)
    
    ### calculating the width
    Q1 = d**2 + lim**2 -cref**2
    Q2 = b**2 + lim**2 -e**2
    Q3 = b**2 + d**2  -a**2

    # Volume_t is the volume of the tetrahedron spanned between P0,P2,P3,P4
    Volume_t = (1/12)*np.sqrt(4*(b**2)*(d**2)*(lim**2) -(b**2)*(Q1**2) -(d**2)*(Q2**2) -(lim**2)*(Q3**2) +Q1*Q2*Q3)
    
    # Area_t is the area of the triangle between P0,P2,P4
    Area_t = (1/4) *np.sqrt((lim+d+cref)*(-lim+d+cref)*(lim-d+cref)*(lim+d-cref))
    
    # w is the width of the kite
    w = 6*(Volume_t/Area_t)

    return w /1e3 # getting it back to meters from mm

### Define Tetrahedron plot function
def plotting_tetrahedron_algorithm(delta_ld,points,variables,label,line_font_size,color):
    ## define a list of u_p values from 0 to 1
    u_p_list = np.arange(0, 1.1, 0.1)

    ## loop through each u_p value and calculate the width
    w_list_tetrahedron = []
    for i in u_p_list:
        w_list_tetrahedron.append(tetrahedron_algorithm(i,delta_ld,points,variables))
    
    ## plot the results
    plt.plot(u_p_list, w_list_tetrahedron, label=label, color=color, linewidth=line_font_size)
    return


### Defining trilaterate, used in trilateration_algorithm
def trilaterate(P1, r1, P2, r2, P3, r3, side):
    # Find the intersection of three spheres
    # P1,P2,P3 are the centers, r1,r2,r3 are the radii
    # Implementaton based on Wikipedia Trilateration article.

    P1, P2, P3 = np.array(P1), np.array(P2), np.array(P3)
    temp1 = P2 - P1
    e_x = temp1 / norm(temp1)  # unit vector from P1 to P2
    temp2 = P3 - P1
    # By multiplying the unit vector lying on the base axis with the vector from to 1 to 3
    # One does a linear transformation to project temp2 onto this base line, thereby finding the
    # x-coordinate of point 3 in the new coordinate system, which is called i
    i = np.dot(e_x, temp2)

    temp3 = temp2 - i * e_x  # Removing the x-coordinate component of the vector from 1 to 3
    if norm(temp3) == 0:  # Getting rid of the division by zeo
        e_y = temp3
        print('error, norm(temp3) = 0')
    else:
        e_y = temp3 / norm(temp3)  # Normalizing that vector

    e_z = np.cross(e_x, e_y)  # Using the cross product to find the last unit vector, perpendicular to both
    j = np.dot(e_y, temp2)  # y-coordinate of P3, calculated in the same fashion

    d = norm(P2 - P1)  # distance between point 1 & 2

    # Actual solution, note that x,y,z are defined in the new coordinate system!
    x = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    if j == 0:  # Getting rid of the division by zeo
        y = (r1 * r1 - r3 * r3 + i * i + j * j - 2 * i * x)
        print('error, j = 0')
    else:
        y = (r1 * r1 - r3 * r3 + i * i + j * j - 2 * i * x) / (2 * j)

    temp4 = r1 * r1 - x * x - y * y
    if temp4 < 0:
        raise Exception("The three spheres do not intersect!")
    z = sqrt(temp4)

    # Factor determines sign of z
    if side == 'right':
        factor = -1
    elif side == 'left':
        factor = 1

    closest_point = P1 + x * e_x + y * e_y + factor * z * e_z

    return np.array([closest_point[0], closest_point[1], closest_point[2]])


### Define Trilateration algorithm
def trilateration_algorithm(u_p,delta_ld,points,variables):
    '''
    --- Trilateration algorithm --- 
    input:  u_p         = [0-1]   (power-setting)
            l_d_percentage = 8 or 13 [%] (depower-tape extension)
            points      = [P0,P1,P2,P3,P4] (P0 is the origin)
            variables   = [a,b,cref,d,e] (a,b,cref,d,e are the line lengths)
    output: w           = width [mm]    (width of the kite)
    '''
    ### defining the points and line lengths
    P0,P1,P2,P3,P4 = points[0],points[1],points[2],points[3],points[4]
    a,b,cref,d,e = variables[0],variables[1],variables[2],variables[3],variables[4]
    
    ### finding the new length by adding the depower-tape extension
    lim  = np.linalg.norm(P0-P4) +up_to_ld(u_p,delta_ld)

    ## Using trigonometry to find the new location of point P4
    q1       = (1/(2*d))*(lim**2 -cref**2 - d**2)
    q2       = np.sqrt(cref**2 - q1**2)
    P4_new  = [q2,0,d+q1]

    P1_new = trilaterate(P0,b, P2, a, P4_new,e,'left')
    P3_new = [P1_new[0],-P1_new[1],P1_new[2]]
    
    ### calculating the width
    w = np.linalg.norm(P1_new-P3_new)
    return w/1e3 # getting it back to meters from mm


### Define Trilateration plot function
def plotting_trilateration_algorithm(delta_ld,points,variables,label,line_font_size,line_style,color):
    ## define a list of u_p values from 0 to 1
    u_p_list = np.arange(0, 1.1, 0.1)
    
    ## loop through each u_p value and calculate the width
    w_list_2plate = []

    for i in u_p_list:
        w_list_2plate.append(trilateration_algorithm(i,delta_ld,points,variables))

    ## plot the results
    plt.plot(u_p_list, w_list_2plate,label=label,color=color,linestyle =line_style,linewidth =line_font_size)
    return


### Photogrammetry measurements
data_photogrammetry_width = np.array([8.2377444 , 8.2114677 , 8.25745192, 8.25088275, 8.2377444 , 7.86407394,
                                7.84413124, 8.0435583 , 7.96378748, 8.00367289, 8.34285119, 8.39540458,
                                8.2377444 , 8.17205266, 8.25745192, 7.7922407 , 7.91707424, 7.68054754,
                                7.68054754, 7.67397736,  ])
data_photogrammetry_width_presimulated = [8.2822583 , 8.25583961, 8.30207231, 8.29546764, 8.2822583,  7.90656865,
                                    7.88651818, 8.08702288, 8.006821  , 8.04692194, 8.38793304, 8.44077042,
                                    8.2822583 , 8.21621158, 8.30207231, 7.83434725, 7.95985535, 7.72205054,
                                    7.72205054, 7.71544485,]
data_photogrammetry_up = [1.,1., 1., 1., 1., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.]

### PLOTTING THE RESULTS

# defining settings
label_font_size  = 12
legend_font_size = 11
dot_size = 2
line_font_size = 4.2
factor_white_dashes = 0.45 

# plotting design CAD shape
plotting_tetrahedron_algorithm(13,points_CAD,lines_CAD,r'Tetrahedron ($\delta l_d$ = 13%), from design geometry',line_font_size,'blue')
plotting_trilateration_algorithm(13,points_CAD,lines_CAD,r'_Trilateration ($\delta l_d$ 13%), from design geometry',factor_white_dashes*line_font_size,'dashed','white')

# plotting presimulated shape 8% 
plotting_tetrahedron_algorithm(8,points_presimulated,lines_presimulated,r'Tetrahedron ($\delta l_d$ = 8%), from presimulation',line_font_size,'red')
plotting_trilateration_algorithm(8,points_presimulated,lines_presimulated,r'_Trilateration ($\delta l_d$ 8%), from presimulation',factor_white_dashes*line_font_size,'dashed','white')

# plotting presimulated shape 8%
plotting_tetrahedron_algorithm(13,points_presimulated,lines_presimulated,r'Tetrahedron ($\delta l_d$ = 13%), from presimulation',line_font_size,'green')
plotting_trilateration_algorithm(13,points_presimulated,lines_presimulated,r'_Trilateration ($\delta l_d$ 13%), from presimulation',factor_white_dashes*line_font_size,'dashed','white')


# plotting photogrammetry measurements
plt.plot(data_photogrammetry_up, data_photogrammetry_width, 'o', color='white',linewidth = dot_size,mec = 'black',label = 'Photogrammetry measurements')

# configuring the plot and saving it
plt.grid(True)
plt.xlabel(r'$u_p\,[-]$',fontsize=label_font_size)
plt.ylabel(r'$w\,[m]$',fontsize=label_font_size)
plt.ylim(7.4,8.6)
plt.legend(fontsize=legend_font_size,loc='lower right', bbox_to_anchor=(1.0, 0.0))
plt.rcParams['pdf.fonttype'] = 42 # Output Type 3 (Type3) or Type 42 (TrueType)
plt.savefig('results_2plate_graph.svg')
plt.savefig('results_2plate_graph.pdf')
plt.show()


