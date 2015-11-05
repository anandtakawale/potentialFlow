import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt

def polarR(x, y):
    """
    (numpy array, numpy array) -> numpy array
    Returns the radius for number of points in the numpy array
    x: X position
    y: Y position
    """
    return np.sqrt(x**2 + y**2)

def polarTheta(x, y):
    """
    (float , float) - > float
    Returns the angle of position vector of x, y w.r.t origin
    """
    return np.arctan2(y , x)

class potential(object):
    def __init__(self, strength):
        """
        Defines the basic characteristics for the potential
        strength: strength of the object
        Strength represents differnt things for differnt potentials
        It is 
        For uniformFlow - velocity [m/s]
        For source/sink, vortex - strength [m^2 / s]
        """
        self.m = strength 
                
class source(potential):
    """
    A typical potential represnting source/sink
    """
    def __init__(self, strength, position):
        """
        position: Origin of the object in the XY plane (a tuple)
        """
        potential.__init__(self, strength)
        self.pos = position
        
    def fi(self, X, Y):
        """
        (self, array, array) -> array
        Returns the potential function values for given X and Y co-ordinates
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return self.m / (2 * pi) * np.where(polarR(X, Y) > 0.1 ,np.log(polarR(X, Y)), 0)
    
    def si(self, X, Y):
        """
        (self, array, array) -> array
        Returns the stream function values for the given X and Y co-ordinates"
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return self.m / (2 * pi) * polarTheta(X, Y)

    def getu(self, X, Y):
        """
        (self, array, array) -> array
        Returns the u velocity for the given field
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return self.m / (2 * pi) * np.where(polarR(X, Y) > 0.2, np.cos(polarTheta(X, Y)) / polarR(X, Y), 0)

    def getv(self, X, Y):
        """
        (self, array, array) -> array
        Returns the verticle component of velocity
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return self.m / (2 * pi) * np.where(polarR(X, Y) > 0.2, np.sin(polarTheta(X, Y)) / polarR(X, Y), 0)
    
    def  __str__(self):
        return "Source with m =" + str(self.m) + " @ " + str(self.pos)

class uniform(potential):
    """
    Represnts uniform flow in particular direction specified
    """
    def __init__(self, velocity, angle):
        """
        (float, float) -> obj
        Initializes uniform flow with following parameters
        velocity: Magnitude of velocity [m/s]
        angle: angle made by velocity wrt X axis [radians]
        """
        potential.__init__(self, velocity)
        self.alpha = angle * pi / 180

    def fi(self, X, Y):
        """
        Returns the value of potential function at a position X and Y
        """
        return self.m * X * np.cos(self.alpha) + self.m * Y * np.sin(self.alpha)
    
    def si(self, X, Y):
        """
        Returns the value of stream function at a position X and Y
        """
        return self.m * -X * np.sin(self.alpha) + self.m * Y * np.cos(self.alpha)

    def getu(self, X, Y):
        """
        Returns the horizontal component of velocity at position (X, Y)
        """
        return self.m * np.cos(self.alpha)

    def getv(self, X, Y):
        """
        Returns the verticle component of velocity at position (X, Y)
        """
        return  self.m * np.sin(self.alpha)

    def __str__(self):
        return "Uniform flow of "+ str(self.m) + "m/s @ " + str(self.alpha * 180 / pi) + " degrees"


class vortex(potential):
    """
    A typical potential represnting vortex
    """
    def __init__(self, strength, position):
        """
        position: Origin of the object in the XY plane (a tuple)
        """
        potential.__init__(self, strength)
        self.pos = position
        
    def si(self, X, Y):
        """
        (self, array, array) -> array
        Returns the stream function values for given X and Y co-ordinates
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return - self.m / (2 * pi) * np.where(polarR(X, Y) > 0.1, np.log(polarR(X, Y)), 0)
    
    def fi(self, X, Y):
        """
        (self, array, array) -> array
        Returns the potential function values for the given X and Y co-ordinates"
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return self.m / (2 * pi) * np.where(polarR(X, Y) > 0.1, polarTheta(X, Y), 0)


    def getu(self, X, Y):
        """
        Returns the horizontal component of velocity at position (X, Y)
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return - self.m / (2 * pi) * np.where(polarR(X, Y) > 0.5, np.sin(polarTheta(X, Y)) / polarR(X, Y), 0)

    def getv(self, X, Y):
        """
        Returns the vertical component of velocity at position (X, Y)
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return self.m / (2 * pi) * np.where(polarR(X, Y) > 0.5, np.cos(polarTheta(X, Y)) / polarR(X, Y), 0)
  
    def  __str__(self):
        return "Vortex with m =" + str(self.m) + " @ " + str(self.pos)

class doublet(potential):
    """
    A potential reprsenting doublet
    """
    def __init__(self, strength, position):
        potential.__init__(self, strength)
        self.pos = position

    def si(self, X, Y):
        """
        Retruns the stream function value at X and Y
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return np.where(polarR(X, Y) > 0.1, -self.m / polarR(X, Y) * np.sin(polarTheta(X, Y)), 0)

    def fi(self, X, Y):
        """
        Returns the potential function value at X and Y
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return np.where(polarR(X, Y) > 0.1, self.m / polarR(X, Y) * np.cos(polarTheta(X, Y)), 0)
    
    def getu(self, X, Y):
        """
        Returns the horizontal component of velocity
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return np.where(polarR(X, Y) > 0.5, -self.m / (polarR(X, Y))**2 * np.cos(2 * polarTheta(X, Y)), 0.0)
    
    def getv(self, X, Y):
        """
        Returns the horizontal component of velocity
        """
        X = X - self.pos[0]
        Y = Y - self.pos[1]
        return np.where(polarR(X, Y) > 0.5, -self.m / (polarR(X, Y))**2 * np.sin(2 * polarTheta(X, Y)), 0.0)
    
    def __str__(self):
        return "Doublet with strength " + str(self.m) + " @ " + str(self.pos)


def contourplot(X, Y, Z, title):
    """
    Plots the contour graph of Z in the XY plane with title
    """
    Zmax = np.max(Z)
    Zmin = np.min(Z)
    V = np.linspace(Zmin, Zmax, 101)
    CP = plt.contour(X, Y, Z, V, linewidth = 2)
    plt.axes().set_aspect('equal', 'datalim')
    plt.colorbar(CP, shrink = 0.8, extend = "both")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(title)

def quiverplot(X, Y, U, V, title):
    """
    Plots the quiverplot 
    """
    plt.quiver(X, Y, U, V, linewidth = 1, color = 'b', units = 'x', width = 0.022)
    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")


def plotsifi(objs, p1=(-5, -5), p2 = (5, 5), fi = False):
    """
    Plots the stream function(by default) and potentialfunction(optional) contours for the given objects in the flow field
    objs: list of objects in the flow field
    p1: bottom left corner of the flow field (x0, y0)
    p2: top right corner of the flow field (x, y)
    """
    x = np.linspace(p1[0], p2[0], 801)
    y = np.linspace(p1[1], p2[1], 801)
    X, Y = np.meshgrid(x, y)
    Zsi = np.zeros(X.shape)
    name = ""
    i = 0
    for obj in objs:
        Zsi = Zsi + obj.si(X, Y)
        if i%2 == 0:
            name += "\n"
        else:
            name+= " "
        name+= obj.__str__()
        i += 1
    contourplot(X, Y, Zsi, "Streamlines in the flowfield with " + name)
    if fi:
        Zfi = np.zeros(X.shape)
        Zfi = Zfi + obj.fi(X, Y)
        contourplot(X, Y, Zfi, "Potential lines in the flowfield")
    plt.show()

def plotquiver(objs, p1=(-5, -5), p2 = (5, 5)):
    """
    Plots the velocity quiver for the given superposition of potential flow
    """
    x = np.linspace(p1[0], p2[0], 50)
    y = np.linspace(p1[1], p2[1], 50)
    X, Y = np.meshgrid(x, y)
    u = np.zeros(X.shape)
    v = np.zeros(X.shape)
    for obj in objs:
        u += obj.getu(X, Y)
        v += obj.getv(X, Y)
    quiverplot(X, Y, u, v, "Quiver plot for velocities")
    plt.show()
