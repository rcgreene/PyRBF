import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt


epsilon = .3 #shape parameter
x = np.linspace(0,np.pi,20)
y = np.linspace(0,np.pi,10)
xx, yy = np.meshgrid(x, y)

f = lambda r: np.exp(-((r)**2)/epsilon)
dist = lambda x_1, y_1, x_2, y_2: np.sqrt((x_1 - x_2)**2 + (y_1 - y_2)**2)

#def g(x_1, x_2):
#    x = np.abs(x_1 - x_2)
#    if (x) > 1:
#        return 0
#    else:
#        return np.exp(-1/(1-x**2) + 1)
z_course = np.sin(xx)*np.sin(yy)

x_fine = np.linspace(0,np.pi,100)
y_fine = np.linspace(0,np.pi,100)
xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)

z_true = [g(c, np.pi/2) for c in x_fine]#np.sin(x_fine)

interp = np.empty((20,20))
for i in range(20):
    for j in range(20):
        interp[i,j] = f(x[i],x[j])


inv_interp = lin.inv(interp)

model = np.dot(inv_interp,y_course)

y_interp = np.zeros(x_fine.shape)

for i in range(len(y_interp)):
    for j in range(len(model)):
        y_interp[i] += model[j]*f(x[j],x_fine[i])

print model

plt.plot(y_true)
plt.show()
