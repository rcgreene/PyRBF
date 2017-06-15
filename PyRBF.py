import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

res_x = 15
res_y = 15
test_func = lambda x,y: np.sin(3*x)*np.sin(y) #function to be approximated
epsilon = .5 #shape parameter
f = lambda r: np.exp(-((r)**2)/epsilon) #kernel function
dist = lambda x_1, y_1, x_2, y_2: np.sqrt((x_1 - x_2)**2 + (y_1 - y_2)**2)


x = np.linspace(0,np.pi,res_x)      #Points used for the fit
y = np.linspace(0,np.pi,res_y)
xx, yy = np.meshgrid(x, y)
z_course = test_func(xx,yy)

x_fine = np.linspace(0,np.pi,100)   #Points used to test the fit
y_fine = np.linspace(0,np.pi,100)
xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)
z_true = test_func(xx_fine,yy_fine)#np.sin(x_fine)

interp = np.empty((res_x*res_y,res_x*res_y))
points = [np.reshape(xx,res_x*res_y),np.reshape(yy,res_x*res_y)]    #Put points in
points_fine = [np.reshape(xx_fine,10000),np.reshape(yy_fine,10000)] #easy-to-use list

for i in range(res_x*res_y):
    for j in range(res_x*res_y):
        #interp i,j = KERNEL(|P[j] - P[i]|)
        interp[i,j] = f(dist(points[0][i],points[1][i],points[0][j],points[1][j]))


inv_interp = lin.inv(interp) #This is used to determine the model
#(NOTE: DIRECT INVERSION CAN HAVE A HIGH CONDITION NUMBER)
print z_course.shape
model = np.dot(inv_interp,np.reshape(z_course,res_x*res_y))
print model
z_interp = np.zeros(10000)

for i in range(len(z_interp)):
    for j in range(len(model)):
        z_interp[i] += model[j]*f(dist(points[0][j],points[1][j],points_fine[0][i],points_fine[1][i]))

fig = plt.figure() #Plot interpolated function
ax = fig.add_subplot(111, projection='3d')
z_interp = np.reshape(z_interp, (100,100))
ax.plot_wireframe(xx_fine,yy_fine, z_interp)
plt.show()

fig = plt.figure() #Plot error
ax = fig.add_subplot(111, projection='3d')
z_interp = np.reshape(z_interp, (100,100))
ax.plot_wireframe(xx_fine,yy_fine, z_interp - z_true)
plt.show()
