import matplotlib.pyplot as plt
import math as m
import numpy as np
from math import sqrt

#input variables for calculating gsd
effective_imaging_element_size = 1.7
focal_length = 40
horizontal_pixel = 2560
vertical_pixel = 2560
altitude = .0255

#convert pixel to mm 
Xsensor = effective_imaging_element_size*0.001*horizontal_pixel
Ysensor = effective_imaging_element_size*0.001*vertical_pixel
pixel_vector = ([Xsensor],[Ysensor],[0])  
print(pixel_vector)

#solve for horizontal and vertical FOV
HFOV = 2*m.atan2(Xsensor,2*focal_length)
VFOV = 2*m.atan2(Ysensor,2*focal_length)
print("your HFOV is",HFOV*(180/m.pi),"degrees")
print("your VFOV is",VFOV*(180/m.pi),"degrees")
   
#solve for horizontal and vertical Dimension at NADIR
Horizontal_Dim = 2*altitude*m.tan(HFOV/2)
Vertical_Dim = 2*altitude*m.tan(VFOV/2)
vector_NADIR = ([Horizontal_Dim],[Vertical_Dim],[-altitude])
vector_intial_middle = ([0],[0],[-altitude])
vector_intial_first = ([Horizontal_Dim/2],[Vertical_Dim/2],[-altitude])
vector_intial_second = ([-Horizontal_Dim/2],[-Vertical_Dim/2],[-altitude])
vector_intial_third = ([-Horizontal_Dim/2],[Vertical_Dim/2],[-altitude])
vector_intial_fourth = ([Horizontal_Dim/2],[-Vertical_Dim/2],[-altitude])
print("your Horizontal Dimension is",Horizontal_Dim,"meters")
print("your Vertical Dimension is",Vertical_Dim,"meters")
print("your Vector at NADIR is", vector_NADIR)

#Euler Rotation Matrices
def Rx(theta):
  return np.matrix([[ 1, 0           , 0           ],
                   [ 0, m.cos(theta),-m.sin(theta)],
                   [ 0, m.sin(theta), m.cos(theta)]])
  
def Ry(theta):
  return np.matrix([[ m.cos(theta), 0, m.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-m.sin(theta), 0, m.cos(theta)]])
  
def Rz(theta):
  return np.matrix([[ m.cos(theta), -m.sin(theta), 0 ],
                   [ m.sin(theta), m.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

#Euler Rotation Input Values for Inertial Sensor
phi = .10*m.pi/180 
theta = 10*m.pi/180
psi = .10*m.pi/180 
print("phi =", phi)
print("theta  =", theta)
print("psi =", psi)
  
#Obtaining rotation matrix based on order rotated, change order depending on rotation
R = Rz(psi) * Ry(theta) * Rx(phi)
print("Your Euler Rotation Matrix =", np.round(R, decimals=2)) 

#New projected vector after Euler Rotation
vector_rotation_middle = R * vector_intial_middle
print("Your vector middle after Euler rotation =", np.round(vector_rotation_middle, decimals=2))
vector_rotation_first = R * vector_intial_first
print("Your vector first after Euler rotation =", np.round(vector_rotation_first, decimals=2))
vector_rotation_second = R * vector_intial_second
print("Your vector second after Euler rotation =", np.round(vector_rotation_second, decimals=2))
vector_rotation_third = R * vector_intial_third
print("Your vector third after Euler rotation =", np.round(vector_rotation_third, decimals=2))
vector_rotation_fourth = R * vector_intial_fourth
print("Your vector fourth after Euler rotation =", np.round(vector_rotation_fourth, decimals=2))


#Assigning direction values using Euler matrix
#middle vector direction
if vector_rotation_middle[0] < 0:
    direction_middle_x = -1
else:
    direction_middle_x = 1
print("your middle x direction is",direction_middle_x)

if vector_rotation_middle[1] < 0:
    direction_middle_y = -1
else:
    direction_middle_y = 1
print("your middle y direction is",direction_middle_y)                

#first vector direction
if vector_rotation_first[0] < 0:
    direction_first_x = -1
else:
    direction_first_x = 1
print("your first x direction is",direction_first_x)

if vector_rotation_first[1] < 0:
    direction_first_y = -1
else:
    direction_first_y = 1
print("your first y direction is",direction_first_y)   

#second vector direction
if vector_rotation_second[0] < 0:
    direction_second_x = -1
else:
    direction_second_x = 1
print("your second x direction is",direction_second_x)

if vector_rotation_second[1] < 0:
    direction_second_y = -1
else:
    direction_second_y = 1
print("your second y direction is",direction_second_y)  

#third vector direction
if vector_rotation_third[0] < 0:
    direction_third_x = -1
else:
    direction_third_x = 1
print("your third x direction is",direction_third_x)

if vector_rotation_third[1] < 0:
    direction_third_y = -1
else:
    direction_third_y = 1
print("your third y direction is",direction_third_y)  

#fourth vector direction
if vector_rotation_fourth[0] < 0:
    direction_fourth_x = -1
else:
    direction_fourth_x = 1
print("your fourth x direction is",direction_fourth_x)

if vector_rotation_fourth[1] < 0:
    direction_fourth_y = -1
else:
    direction_fourth_y = 1
print("your fourth y direction is",direction_fourth_y)  

#Solving for magnitude of vectors
#If just a theta or psi rotation
if ((phi == 0) and (-m.pi/2 <= theta <= m.pi/2) and (-2*m.pi <= psi <= 2*m.pi)):
    x_middle = abs(m.tan(theta)*altitude)*direction_middle_x
    x_hyp_middle = sqrt(x_middle**2 + altitude**2)
    y_middle = abs(x_hyp_middle*m.tan(VFOV/2))*direction_middle_y
    vector_middle = ([x_middle],[y_middle],[-altitude])
    print("your new middle vector =", vector_middle)
    
    x_1 = abs(m.tan(theta+(HFOV/2))*altitude)*direction_first_x
    x_hyp1 = sqrt(x_1**2 + altitude**2)
    y_1 = abs(x_hyp1*m.tan(VFOV/2))*direction_first_y
    vector_first = ([x_1],[y_1],[-altitude])
    print("your new first vector =", vector_first)
    
    x_2 = abs(m.tan(theta-(HFOV/2))*altitude)*direction_second_x
    x_hyp2 = sqrt(x_2**2 + altitude**2)
    y_2 = abs(x_hyp2*m.tan(VFOV/2))*direction_second_y
    vector_second = ([x_2],[y_2],[-altitude])
    print("your new second vector =", vector_second)
    
    x_3 = abs(m.tan(theta-(HFOV/2))*altitude)*direction_third_x
    x_hyp3 = sqrt(x_3**2 + altitude**2)
    y_3 = abs(x_hyp3*m.tan(VFOV/2))*direction_third_y
    vector_third = ([x_3],[y_3],[-altitude])
    print("your new third vector =", vector_third)
   
    x_4 = abs(m.tan(theta+(HFOV/2))*altitude)*direction_fourth_x
    x_hyp4 = sqrt(x_4**2 + altitude**2)
    y_4 = abs(x_hyp4*m.tan(VFOV/2))*direction_fourth_y
    vector_fourth = ([x_4],[y_4],[-altitude])
    print("your new fourth vector =", vector_fourth)
    
    #3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #axis
    ax.quiver(0, 0, 0, 0, 0, 2, color='grey', linestyle="dashed")
    ax.quiver(0, 0, 0, 0, 4, 0, color='grey', linestyle="dashed")
    ax.quiver(0, 0, 0, 4, 0, 0, color='grey', linestyle="dashed")
    #Initial vectors
    ax.quiver(0, 0, 0, vector_NADIR[0], vector_NADIR[1], vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, Horizontal_Dim/2, Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, -Horizontal_Dim/2, -Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, -Horizontal_Dim/2, Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, Horizontal_Dim/2, -Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    #showing 3 points just for theta rotation
    ax.quiver(0, 0, 0, vector_middle[0], vector_middle[1], vector_middle[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_first[0], vector_first[1], vector_first[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_second[0], vector_second[1], vector_second[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_third[0], vector_third[1], vector_third[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_fourth[0], vector_fourth[1], vector_fourth[2], arrow_length_ratio=0, color='blue')
    
    #axis limits
    ax.set_xlim([-4, 4])
    ax.set_ylim([-4, 4])
    ax.set_zlim([-4, 0])
    plt.savefig("fov.png")
    plt.show()
    
    #2D plot
    x_coordinates = [vector_middle[0],vector_first[0],vector_second[0],vector_third[0],vector_fourth[0]]
    y_coordinates = [vector_middle[1],vector_first[1],vector_second[1],vector_third[1],vector_fourth[1]]
    plt.scatter(x_coordinates, y_coordinates)
    
#If just a phi or psi rotation
elif ((-m.pi/2 <= phi <= m.pi/2) and (theta == 0) and (-2*m.pi <= psi <= 2*m.pi)):
    y_middle = abs(m.tan(phi)*altitude)*direction_middle_y
    y_hyp_middle = sqrt(y_middle**2 + altitude**2)
    x_middle = abs(y_hyp_middle*m.tan(HFOV/2))*direction_middle_x
    vector_middle = ([x_middle],[y_middle],[-altitude])
    print("your new middle vector =", vector_middle)
    
    y_1 = abs(m.tan(phi+(VFOV/2))*altitude)*direction_first_y
    y_hyp1 = sqrt(y_1**2 + altitude**2)
    x_1 = abs(y_hyp1*m.tan(HFOV/2))*direction_first_x
    vector_first = ([x_1],[y_1],[-altitude])
    print("your new first vector =", vector_first)
    
    y_2 = abs(m.tan(phi-(VFOV/2))*altitude)*direction_second_y
    y_hyp2 = sqrt(y_2**2 + altitude**2)
    x_2 = abs(y_hyp2*m.tan(HFOV/2))*direction_second_x
    vector_second = ([x_2],[y_2],[-altitude])
    print("your new second vector =", vector_second)    
    
    y_3 = abs(m.tan(phi-(VFOV/2))*altitude)*direction_third_y
    y_hyp3 = sqrt(y_3**2 + altitude**2)
    x_3 = abs(y_hyp3*m.tan(HFOV/2))*direction_third_x
    vector_third = ([x_3],[y_3],[-altitude])
    print("your new third vector =", vector_third) 

    y_4 = abs(m.tan(phi+(VFOV/2))*altitude)*direction_fourth_y
    y_hyp4 = sqrt(y_4**2 + altitude**2)
    x_4 = abs(y_hyp4*m.tan(HFOV/2))*direction_fourth_x
    vector_fourth = ([x_4],[y_4],[-altitude])
    print("your new fourth vector =", vector_fourth) 

    #3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #axis
    ax.quiver(0, 0, 0, 0, 0, 2, color='grey', linestyle="dashed")
    ax.quiver(0, 0, 0, 0, 4, 0, color='grey', linestyle="dashed")
    ax.quiver(0, 0, 0, 4, 0, 0, color='grey', linestyle="dashed")
    #Initial vectors
    ax.quiver(0, 0, 0, vector_NADIR[0], vector_NADIR[1], vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, Horizontal_Dim/2, Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, -Horizontal_Dim/2, -Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, -Horizontal_Dim/2, Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, Horizontal_Dim/2, -Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    #showing 3 points just for theta rotation
    ax.quiver(0, 0, 0, -1*x_middle, vector_middle[1], vector_middle[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, -1*x_1, vector_first[1], vector_first[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, -1*x_2, vector_second[1], vector_second[2], arrow_length_ratio=0, color='blue')
    #changing y to negative
    ax.quiver(0, 0, 0, vector_middle[0], vector_middle[1], vector_middle[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_first[0], vector_first[1], vector_first[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_second[0], vector_second[1], vector_second[2], arrow_length_ratio=0, color='blue')
    #axis limits
    ax.set_xlim([-4, 4])
    ax.set_ylim([-4, 4])
    ax.set_zlim([-4, 0])
    plt.savefig("fov.png")
    plt.show()

    #2D plot
    x_coordinates = [vector_middle[0],vector_first[0],vector_second[0],vector_third[0],vector_fourth[0]]
    y_coordinates = [vector_middle[1],vector_first[1],vector_second[1],vector_third[1],vector_fourth[1]]
    plt.scatter(x_coordinates, y_coordinates) 


#Compound rotation of theta,phi,psi
elif ((-m.pi/2 <= phi <= m.pi/2) and (-m.pi/2 <= theta <= m.pi/2) and (-2*m.pi <= psi <= 2*m.pi)):
    x_middle = abs(m.tan(theta)*altitude)*direction_middle_x
    x_1 = abs(m.tan(theta+(HFOV/2))*altitude)*direction_first_x
    x_2 = abs(m.tan(theta-(HFOV/2))*altitude)*direction_second_x
    x_3 = abs(m.tan(theta-(HFOV/2))*altitude)*direction_third_x
    x_4 = abs(m.tan(theta+(HFOV/2))*altitude)*direction_fourth_x
    
    y_middle = abs(m.tan(phi)*altitude)*direction_middle_y
    y_1 = abs(m.tan(phi+(VFOV/2))*altitude)*direction_first_y
    y_2 = abs(m.tan(phi-(VFOV/2))*altitude)*direction_second_y
    y_3 = abs(m.tan(phi-(VFOV/2))*altitude)*direction_third_y
    y_4 = abs(m.tan(phi+(VFOV/2))*altitude)*direction_fourth_y
    
    vector_middle = ([x_middle],[y_middle],[-altitude])
    vector_first = ([x_1],[y_1],[-altitude])
    vector_second = ([x_2],[y_2],[-altitude])
    vector_third = ([x_3],[y_3],[-altitude])
    vector_fourth = ([x_4],[y_4],[-altitude])
    
    print("your new middle vector =", vector_middle)        
    print("your new first vector =", vector_first)
    print("your new second vector =", vector_second)
    print("your new third vector =", vector_third) 
    print("your new fourth vector =", vector_fourth) 

    x_hyp1 = sqrt(x_1**2 + altitude**2)
    y_1 = abs(x_hyp1*m.tan(VFOV/2))*direction_first_y
    vector_first = ([x_1],[y_1],[-altitude])
    print("your new first vector =", vector_first)

   
    #3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #axis
    # ax.quiver(0, 0, 0, 0, 0, .5, color='orange', linestyle="dashed")
    # ax.quiver(0, 0, 0, 0, .5, 0, color='green', linestyle="dashed")
    # ax.quiver(0, 0, 0, .5, 0, 0, color='purple', linestyle="dashed")
    #Initial vectors
    ax.quiver(0, 0, 0, Horizontal_Dim/2, Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, -Horizontal_Dim/2, -Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, -Horizontal_Dim/2, Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    ax.quiver(0, 0, 0, Horizontal_Dim/2, -Vertical_Dim/2, vector_NADIR[2], arrow_length_ratio=0, color='red')
    
    #ax.quiver(0, 0, 0, vector_middle[0], vector_middle[1], vector_middle[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_first[0], vector_first[1], vector_first[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_second[0], vector_second[1], vector_second[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_third[0], vector_third[1], vector_third[2], arrow_length_ratio=0, color='blue')
    ax.quiver(0, 0, 0, vector_fourth[0], vector_fourth[1], vector_fourth[2], arrow_length_ratio=0, color='blue')
    #axis limits
    ax.set_xlim([-.006, .006])
    ax.set_ylim([-.006, .006])
    ax.set_zlim([-.0255,0])
    plt.savefig("fov.png")
    plt.show()

    # #2D plot
    # plt.plot([vector_second[0],vector_middle[0],vector_first[0]],[vector_second[1],vector_middle[1],vector_first[1]])
    # plt.plot([x_middle,x_1,x_2],[y_middle*-1,y_1*-1,y_2*-1])
    # plt.plot([x_2,x_2],[y_2*-1,y_2])
    # plt.plot([x_1,x_1],[y_1*-1,y_1]) 
    
    #2D plot
    # x_coordinates = [vector_middle[0],vector_first[0],vector_second[0],vector_third[0],vector_fourth[0]]
    # y_coordinates = [vector_middle[1],vector_first[1],vector_second[1],vector_third[1],vector_fourth[1]]
    # plt.xlim(-8000,8000)
    # plt.ylim(-8000,8000)
    # plt.scatter(x_coordinates, y_coordinates)
    
else:
    print("something went wrong not sure what")
    






