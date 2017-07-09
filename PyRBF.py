import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def gauss(r, epsilon):
    return np.exp(-((r)**2)*epsilon) #kernel function

def dist(x_1,y_1,x_2,y_2):
    return np.sqrt((x_1 - x_2)**2 + (y_1 - y_2)**2)

class rbf_interpolator:
    def __init__(self, epsilon = 1, rbf = gauss):
        self.epsilon = epsilon
        self.rbf = rbf

    def eval_rbf(self, r):
        return self.rbf(r, self.epsilon)

    def create_eval_matrix(self, points_1, points_2, func = None):
        #create RBF interpolation matrix from points
        if not func:
            func = self.eval_rbf
        res_1 = len(points_1[0])
        res_2 = len(points_2[0])
        eval_mat = np.empty((res_1,res_2))
        for i in range(res_1):
            for j in range(res_2):
                r = dist(points_1[0][i],points_1[1][i],points_2[0][j],points_2[1][j])
                eval_mat[i,j] = func(r)
        return eval_mat

    def determine_model(self, points, values):
        #determines rbf_model for a function which takes on 'values' at 'points'
        self.model_points = points
        interp = self.create_eval_matrix(points, points)
        self.model = lin.solve(interp, values)

    def evaluate_model(self, eval_points):
        #interpolate based on model_points
        if not hasattr(self, 'model'):
            print "must determine model first"
            return
        eval_mat = self.create_eval_matrix(eval_points, self.model_points)
        return np.dot(eval_mat, self.model)

    def plot_model(self, res, x_range = [0,np.pi], y_range = [0,np.pi]):
        x = np.linspace(x_range[0], x_range[1], res)
        y = np.linspace(y_range[0], y_range[1], res)
        XX,YY = np.meshgrid(x,y)
        plot_points = [np.reshape(XX,res**2), np.reshape(YY,res**2)]
        self.z_calc = self.evaluate_model(plot_points)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        self.z_calc = np.reshape(self.z_calc, (res,res))
        ax.plot_wireframe(XX, YY, self.z_calc)
        plt.show()

def test(epsilon, res):
    #interpolates a simple test function using a gaussian with shape parameter epsilon over a grid with res by res points
    res_x = res
    res_y = res
    test_func = lambda x,y: np.sin(np.pi*(np.cos(x)+1)/2)*np.sin(np.pi*(np.cos(y)+1)/2) #function to be approximated


    x = np.linspace(0,np.pi,res_x)      #Points used for the fit
    y = np.linspace(0,np.pi,res_y)
    xx, yy = np.meshgrid(x, y)
    z_course = test_func(xx,yy)

    x_fine = np.linspace(0,np.pi,100)   #Points used to test the fit
    y_fine = np.linspace(0,np.pi,100)
    xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)
    z_true = test_func(xx_fine,yy_fine)#np.sin(x_fine)

    points = [np.reshape(xx,res_x*res_y),np.reshape(yy,res_x*res_y)]    #Put points in
    points_fine = [np.reshape(xx_fine,10000),np.reshape(yy_fine,10000)] #easy-to-use list

    interpolator = rbf_interpolator(epsilon = epsilon)

    interpolator.determine_model(points, np.reshape(z_course, res_x*res_y))
    interpolator.plot_model(100)

    fig = plt.figure() #Plot error
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(xx_fine,yy_fine, interpolator.z_calc - z_true)
    plt.show()
