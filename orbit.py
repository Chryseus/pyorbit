import math
import numpy
import os
import re
import csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# ============
# === Data ===
# ============

G: float = 6.67430e-11 # gravitational constant m3 kg-1 s-2 
epoch = 2451545.0 # J2000 epoch - 1 January 2000 - 12:00:00 UTC

Sun = {"μ": 1.32712440018e20,
       "mass": 1.9885e30,
       "radius": 695700}

# Planets
Earth = {"μ": 3.986004418e14,
          "mass": 5.97237e24,
          "radius": 6378.1366}
Mercury = {"μ": 2.2032e13,
          "mass": 0.330114e24,
          "radius": 2440.53}
Venus = {"μ": 3.24859e14,
          "mass": 4.86747e24,
          "radius": 6051.8}
Mars = {"μ": 4.282837e13,
          "mass": 0.641712e24,
          "radius": 3396.19}
Ceres = {"μ": 6.26325e10,
          "mass": 0.00093835e24,
          "radius": 469.73}
Jupiter = {"μ": 1.26686534e17,
          "mass": 1898.187e24,
          "radius": 71492}
Saturn = {"μ": 3.7931187e16,
          "mass": 568.3174e24,
          "radius": 60268}
Uranus = {"μ": 5.793939e15,
          "mass": 86.8127e24,
          "radius": 25559}
Neptune = {"μ": 6.836529e15,
          "mass": 102.4126e24,
          "radius": 24764}
Pluto = {"μ": 8.71e11,
          "mass": 0.013030e24,
          "radius": 1188.3}

# Moons
Moon = {"μ": 4.9048695e12,
          "mass": 0.07342e24,
          "radius": 1738.1,
          "parent": "Earth" }

# Other objects

# =========================
# === Utility Functions ===
# =========================

class vec3:
    def __init__(self, x, y , z):
        self.x = x
        self.y = y
        self.z = z

def calcEccentricAnomalyNewton(M, e, epsilon):
    "Calculate eccentric anomaly using Newton's method"
    # M = Mean anomaly in degrees
    # e = Eccentricity
    # epsilon = Target accuracy
    i = 0
    M = math.radians(M)
    if e < 0.8:
        E = M
    else:
        E = math.pi
    
    while E != 0:
        F = E - e * math.sin(E) - M # function
        Fp = 1 - e * math.cos(E)    # derivative
        d = F / Fp
        d2 = E - d
        E = d2
        i = i + 1
        
        if i >= 40:
            print(f"Failed to converge after {i} iterations.")
            break
        if math.fabs(F) < epsilon: # Target accuracy has been reached
            E = math.degrees(E)
            return E
    return E

def calcJDN_GC(Y, M, D):
    "Calculate Julian day number from a Gregorian calendar date"
    a = (14 - M)//12
    y = Y + 4800 - a
    m = M + 12*a - 3
    JDN =  D + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045
    return JDN
def calcJDN_JC(Y, M, D):
    "Calculate Julian day number from a Julian calendar date"
    a = (14 - M)//12
    y = Y + 4800 - a
    m = M + 12*a - 3
    JDN = D + (153*m+2)//5 + y*365 + y//4 - 32083
    return JDN
def calcJD(JDN, H, M, S):
    "Calculate Julian date from Julian day number and time"
    JD = JDN + ((H-12)/24) + (M/1440) + (S/86400)
    return JD
def JDDeltaTimeSeconds(JD1, JD2):
    "Calculate time difference in seconds between two Julian dates"
    delta = JD1 - JD2
    return delta*86400

def calcSMA(ApA, PeA, r):
    a = ((ApA+PeA) / 2) + r
    return a
def calcSMA_c(ApC, PeC):
     a = ((ApA+PeA) / 2)
     return a


# ==================================================================================================    

class Body:
    def __init__(self, name, mass, radius, μ = Sun["μ"]):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.μ = μ

    def getKeplerianDebug(self):
        print(f"a = {self.a/1000} km")
        print(f"e = {self.e}")
        print(f"i = {self.i} degrees")
        print(f"Ω = {self.Ω} degrees")
        print(f"ω = {self.ω} degrees")
        print(f"M0 = {self.M0} degrees")
        print(f"E = {self.E} degrees")
        print(f"v = {self.v} degrees")

    def loadKeplerian(self,a,e,i,Ω,ω,M0=0,t=epoch,M=0,v=0,E=0):
        "Load Keplerian elements a, e, i, Ω, ω and JD time t since epoch"
        self.a = a*1000 # Semi-major axis - km (convert to meters)
        self.e = e # Eccentricity
        self.i = i # Inclination - degrees
        self.Ω = Ω # Right Ascension of Ascending Node (RAAN/LAN) - degrees
        self.ω = ω # Argument of Periapsis - degrees
        self.t = t # time from epoch - JD
        self.M = M # Mean Anomaly - degrees
        self.M0 = M0 # Mean Anomaly at epoch - degrees
        self.v = v # True Anomaly - degrees
        self.E = E # Eccentric Anomaly - degrees
        return 0

    def loadKeplerianHorizons(self, file):
        print(f"Loading Keplerian data for {self.name} from {file}")
        with open(file,'r') as handle:
            data = handle.read()
        m = re.search('(?<=\$\$SOE)(.*?)(?=\$\$EOE)',data,re.DOTALL | re.MULTILINE)
        string = m.group()
        string = string.strip()
        string = string.split(',')
        e = float(string[2])
        i = float(string[4])
        Ω = float(string[5])
        ω = float(string[6])
        M0 = float(string[9])
        a = float(string[11])
        self.loadKeplerian(a,e,i,Ω,ω,M0)

    def loadStateVectors(self, position, velocity):
        # TODO
        return 0

    def calcPericenter(self):
        "Calculate the pericentre distance PeC in meters"
        self.PeC = (1 - self.e) * (self.a)
        return self.PeC
    def calcApocenter(self):
        "Calculate the apocentre ApC distance in meters"
        self.ApC = (1 + self.e) * (self.a)
        return self.ApC
    def calcApoapsis(self):
        "Calculate the apoapsis ApA in meters"
        self.ApA = self.calcApocenter() - (Sun["radius"]*1000)
        return self.ApA
    def calcPeriapsis(self):
        "Calculate periapsis PeA in meters"
        self.PeA = self.calcPericenter() - (Sun["radius"]*1000)
        return self.PeA


    def getMeanMotion(self):
        "Calculate mean motion n in degrees per second"
        self.n = math.degrees(math.sqrt(self.μ/self.a**3))
        return self.n

    def calcMeanAnomalyFromEpoch(self):
        "Calculate mean anomaly M since the epoch and JD time t"
        self.M = self.M0 + self.getMeanMotion() * JDDeltaTimeSeconds(self.t,epoch)
        if self.M < 0:
            self.M = self.M % 360
        if self.M > 360:
            self.M = self.M % 360
        return self.M

    def calcEccentricAnomaly(self):
        "Calculate eccentric anomaly from mean anomaly M"
        # TODO add a tidy way for the user to set tolerance
        self.E = calcEccentricAnomalyNewton(self.M, self.e, 1e-15)
        return self.E
    
    def calcTrueAnomaly(self):
        "Calculate true anomaly from eccentric anomaly E"
        E = math.radians(self.E)
        e = self.e
        #v = 2 * math.atan(math.sqrt((1+e)/(1-e)) * math.tan(E/2))
        v = 2 * math.atan2(math.sqrt(1+e) * math.sin(E/2), math.sqrt(1-e) * math.cos(E/2))
        self.v = math.degrees(v)
        return v
    def calcAngularMomentum(self):
        "Calculates the angular momentum h"
        # TODO Check equation is correct
        self.h = math.sqrt(self.μ * self.a * (1 - self.e**2))
        #self.h = self.mass * self.calcVelocityPeriapsis() * self.calcPericenter()
        return self.h
    def calcVelocityPeriapsis(self):
        "Calculate velocity at periapsis"
        self.vp = math.sqrt(((1+self.e) * self.μ) / ((1-self.e) * self.a))
        return self.vp
    def calcVelocityApoapsis(self):
        "Calculate velocity at apoapsis va"
        self.va = math.sqrt(((1-self.e) * self.μ) / ((1+self.e) * self.a))
        return self.va
    def calcVelocityMean(self):
        "Calculates the mean velocity vm"
        self.vm = math.sqrt(self.μ / self.a) # m/s
        return self.vm
    def calcPeriod(self):
        "Calculates the orbital period p in seconds"
        self.p = 2 * math.pi * math.sqrt(self.a**3 / self.μ)
        return self.p
    def keplerianToStateVectors(self, parent=None):
        "Converts Keplerian elements to state vectors"
        self.position = vec3(0,0,0)
        self.velocity = vec3(0,0,0)
        
        # Convert to radians
        E = math.radians(self.E)
        Ω = math.radians(self.Ω)
        ω = math.radians(self.ω)
        inc = math.radians(self.i)
        v = math.radians(self.v)

        r = self.a * ((1-self.e**2) / (1+self.e * math.cos(v)))

        # Convert to rectangular coordinates
        x = r * math.cos(v)
        y = r * math.sin(v)
        z = 0
        # Local velocity vector
        rp = math.sqrt(self.μ*self.a) / r
        xd = rp * -(math.sin(E))
        yd = rp * math.sqrt(1-self.e**2) * math.cos(E)
        zd = 0

        o = numpy.array([[x],[y],[z]])
        odot = numpy.array([[xd],[yd],[zd]])
        
        # Matrices
        R3 = numpy.array([[ math.cos(Ω), -math.sin(Ω), 0 ],
                        [ math.sin(Ω), math.cos(Ω), 0 ],
                        [ 0, 0, 1 ]])
        R2 = numpy.array([[ 1, 0, 0 ],
                        [ 0, math.cos(inc), -math.sin(inc) ],
                        [ 0, math.sin(inc), math.cos(inc) ]])                
        R1 = numpy.array([[ math.cos(ω), -math.sin(ω), 0 ],
                        [ math.sin(ω), math.cos(ω), 0 ],
                        [ 0, 0, 1 ]])
        # Matrix multiplication
        Rsum = numpy.matmul(R1, o)
        Rsum2 = numpy.matmul(R2, Rsum)
        Rsum3 = numpy.matmul(R3, Rsum2)

        if(parent != None): # Move body centre to parent position
            self.position.x = Rsum3[0,0] + parent.position.x
            self.position.y = Rsum3[1,0] + parent.position.y
            self.position.z = Rsum3[2,0] + parent.position.z
        else:
            self.position.x = Rsum3[0,0]
            self.position.y = Rsum3[1,0]
            self.position.z = Rsum3[2,0]

        Rsum = numpy.matmul(R1, odot)
        Rsum2 = numpy.matmul(R2, Rsum)
        Rsum3 = numpy.matmul(R3, Rsum2)

        self.velocity.x = Rsum3[0,0]
        self.velocity.y = Rsum3[1,0]
        self.velocity.z = Rsum3[2,0]

        return 0
    def StateVectorsToKeplerian(self):
        return 0

def getOrbitData(body, points, parent=None):
    "Helper function to get position and velocity data for a body over one complete orbit"
    # TODO better data packing
    positions = []
    velocities = []
    body.M = 0
    step = 360/(points)
    for i in range(0,points+1):
        body.M = i*step
        body.calcEccentricAnomaly()
        body.calcTrueAnomaly()
        body.keplerianToStateVectors(parent)
        p = [body.position.x, body.position.y, body.position.z]
        v = [body.velocity.x, body.velocity.y, body.velocity.z]
        positions.append(p)
        velocities.append(v)
    return positions, velocities

def getPVatDate(body, JD, parent):
    body.t = JD
    body.calcMeanAnomalyFromEpoch()
    body.calcEccentricAnomaly()
    body.calcTrueAnomaly()
    body.keplerianToStateVectors(parent)
    return body.position.x, body.position.y, body.position.z, body.velocity.x, body.velocity.y, body.velocity.z

def unpackData(data):
    X = []
    Y = []
    Z = []
    for i in range(0,len(data)):
        X.append(data[i][0])
        Y.append(data[i][1])
        Z.append(data[i][2])
    return X, Y, Z

def plotBody(body, points, JD, parent=None):
    p, v = getOrbitData(body,points, parent)
    X, Y, Z = unpackData(p)
    Xd, Yd, Zd = unpackData(v)
    plt.plot(X,Y,Z)
    Xi, Yi, Zi, Xiv, Yiv, Ziv = getPVatDate(body,JD, parent)
    plt.plot(Xi,Yi,Zi,'bo')
    plt.quiver(Xi,Yi,Zi,Xiv,Yiv,Ziv, length=2e6, color='red')


os.system("cls")
earth = Body("Earth", Earth["mass"], Earth["radius"])
mercury = Body("Mercury", Mercury["mass"], Mercury["radius"])
venus = Body("Venus", Venus["mass"], Venus["radius"])
mars = Body("Mars", Mars["mass"], Mars["radius"])
halley1p = Body("Halley's Comet", 5e6, 5.5)
moon = Body("Moon",Moon["mass"],Moon["radius"])

earth.loadKeplerianHorizons("data/earth.txt")
mercury.loadKeplerianHorizons("data/mercury.txt")
venus.loadKeplerianHorizons("data/venus.txt")
mars.loadKeplerianHorizons("data/mars.txt")
halley1p.loadKeplerianHorizons("data/halley1p.txt")
moon.loadKeplerianHorizons("data/moon.txt")
moon.μ = Earth["μ"]

time = calcJD(calcJDN_GC(2021,3,27),12,0,0)


fig = plt.figure()
ax = plt.axes(projection = '3d')
plt.axis([-2e11,2e11,-2e11,2e11])
ax.set_zlim(-1.5e11,1.5e11)
plotBody(earth, 100, time)
plotBody(mercury, 50, time)
plotBody(venus, 100, time)
plotBody(mars, 100, time)
plotBody(halley1p, 2000, time)
plotBody(moon,100,time, earth)
moon.getKeplerianDebug()
plt.plot(0,0,0,'yo') # sun
plt.show()

