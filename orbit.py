import math
import numpy
import os
from engineering_notation import EngNumber
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

G = 6.67430e-11 # m3 kg-1 s-2 (gravitational constant)
Sun_Mass = 1.9885e30 # kg
μ = G*Sun_Mass

class vec3:
    def __init__(self, x, y , z):
        self.x = x
        self.y = y
        self.z = z

"""def calc_single_body_acceleration(bodies, body_index):
    acceleration = vec3(0,0,0)
    target_body = bodies[body_index]
    for index, external_body in enumerate(bodies):
        if index != body_index:
            r = (target_body.position.x - external_body.position.x)**2 + (target_body.position.y - external_body.position.y)**2 + (target_body.position.z - external_body.position.z)**2
            r = math.sqrt(r)
            tmp = G * external_body.mass / r**3
            acceleration.x += tmp * (external_body.position.x - target_body.position.x)
            acceleration.y += tmp * (external_body.position.y - target_body.position.y)
            acceleration.z += tmp * (external_body.position.z - target_body.position.z)
    return acceleration """


def calcEccentricAnomalyBisection(M, e, a, b, i, tol):
    # Calculate eccentric anomaly via bisection method
    # M = Mean anomaly in degrees
    # e = Eccentricity
    # a,b = interval [a,b]
    # i = number of iterations
    # tol = Target accuracy
    M = math.radians(M)
    f = lambda E: E - e * math.sin(E) - M # function
    
    if (f(a)*f(b)) >= 0: # ensure that there is a zero crossing in the interval
        print("Bisection failed.")
       
    a_n = a
    b_n = b
    for n in range(1, i):
        m_n = (a_n + b_n)/2 # midpoint
        print("bisection = " + str(m_n))
        f_m_n = f(m_n) # calculate function on midpoint
        print("f(bisection)" + str(f_m_n))
        print("a = " + str(a_n))
        print("b = " + str(b_n))
        print("f(a) = " + str(f(a_n)))
        print("f(b) = " + str(f(b_n)) + "\n")
        if (b_n - a_n) /2 < tol: # check if result is within tolerance
            print("Tolerance met after " + str(n) + " iterations.")
            return math.degrees(m_n)
        elif f(a_n)*f_m_n < 0: # new interval
                a_n = a_n
                b_n = m_n
        elif f(b_n)*f_m_n < 0: # # new interval
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0.00: # found exact result
            print("Found exact solution")
            return m_n
        else:
            print("Bisection method failed.")
            return None
    result = (a_n + b_n)/2 # return the midpoint we don't have the exact result
    return math.degrees(result)


def calcEccentricAnomalyNewton(M, e, epsilon):
    # Calculate eccentric anomaly via Newton's method
    # M = Mean anomaly in degrees
    # e = Eccentricity
    # epsilon = Target accuracy
    i = 0
    M = math.radians(M)
    #e = math.radians(e)
    if e < 0.8:
        E = M
    else:
        E = math.pi
    
    #print("Solving eccentric anomaly.")
    #print("E = " + str(E))
    #print("e = " + str(e))
    #print("M = " + str(M))
    while E != 0:
        F = E - e * math.sin(E) - M # function
        Fp = 1 - e * math.cos(E)    # derivative
        d = F / Fp
        d2 = E - d
        E = d2
        i = i + 1
        
        #print("Iteration " + str(i))
        #print("f(E) = " + str(F))
        #print("f'(E) = " + str(Fp))
        #print("F(E) / F'(E) " + str(d))
        #print("E - (F(E) / F'(E)) = " + str(d2) + "\n")
        
        if i >= 40:
            print("Failed to converge.")
            break
        if math.fabs(F) < epsilon:
            E = math.degrees(E)
            #print("Eccentric anomaly " + str(E))
            return E
    return E

# ==================================================================================================    

#T = 2451545.0 # J2000 epoch
T = 0
M0 = 3.586172562406391e2 # Mean anomaly at epoch
Sun_Radius = 695508 # km

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

class Body:
    def __init__(self, name, mass, radius):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.μ = G * mass

    def loadKeplarian(self,a,e,i,Ω,ω,t=T,M=0,v=0,E=0):
        "Load available Keplerian elements"
        self.a = a*1000 # Semi-major axis - km (convert to meters)
        self.e = e # Eccentricity
        self.i = i # Inclination - degrees
        self.Ω = Ω # Right Ascension of Ascending Node (RAAN/LAN) - degrees
        self.ω = ω # Argument of Periapsis - degrees
        self.t = t # time from epoch - JD
        self.M = M # Mean Anomaly - degrees
        self.v = v # True Anomaly - degrees
        self.E = E # Eccentric Anomaly - degrees
        return 0
    def loadStateVectors(self, position, velocity):
        return 0
    def calcPericenter(self):
        "Calculate the pericentre distance in meters"
        self.PeC = (1 - self.e) * (self.a) # m
        return self.PeC
    def calcApocenter(self):
        "Calculate the apocentre distance in meters"
        self.ApC = (1 + self.e) * (self.a) # m
        return self.ApC
    def calcApoapsis(self):
        "Calculate the apoapsis in meters"
        self.ApA = self.calcApocenter() - (Sun_Radius/1000)
        return self.ApA
    def calcPeriapsis(self):
        "Calculate periapsis in meters"
        self.PeA = self.calcPericenter() - (Sun_Radius/1000)
        return self.PeA
    def getMeanMotion(self):
        "Calculate mean motion in degrees per second"
        #self.n = math.degrees(math.sqrt(G*(Sun_Mass+self.mass)/self.a**3))
        self.n = math.degrees(math.sqrt(G*(Sun_Mass+self.mass)/self.a**3))
        return self.n
    def calcMeanAnomalyFromEpoch(self):
        "Calculate mean anomaly at time t"
        self.M = M0 + self.getMeanMotion() * (self.t - T)
        if self.M > 360:
            self.M = self.M % 360
        return self.M
    def calcEccentricAnomaly(self):
        "Calculate eccentric anomaly from mean anomaly"
        self.E = calcEccentricAnomalyNewton(self.M, self.e, 1e-15)
        return self.E
    def calcTrueAnomaly(self):
        "Calculate true anomaly from eccentric anomaly"
        E = math.radians(self.E)
        e = self.e
        #v = 2 * math.atan(math.sqrt((1+e)/(1-e)) * math.tan(E/2))
        v = 2 * math.atan2(math.sqrt(1+e) * math.sin(E/2), math.sqrt(1-e) * math.cos(E/2))
        self.v = math.degrees(v)
        return v
    def calcAngularMomentum(self):
        self.h = math.sqrt(μ * self.a * (1 - self.e**2))
        #self.h = self.mass * self.calcVelocityPeriapsis() * self.calcPericenter()
        return self.h
    def calcVelocityPeriapsis(self):
        "Calculate velocity at periapsis"
        t_a = self.a # m
        self.vp = math.sqrt(((1+self.e) * μ) / ((1-self.e) * t_a))
        return self.vp
    def calcVelocityApoapsis(self):
        "Calculate velocity at apoapsis"
        t_a = self.a # m
        self.va = math.sqrt(((1-self.e) * μ) / ((1+self.e) * t_a))
        return self.va
    def calcVelocityMean(self):
        self.vm = math.sqrt((G * Sun_Mass) / self.a) # m/s
        return self.vm
    def calcPeriod(self):
        self.p = 2 * math.pi * math.sqrt(self.a**3 / (G * Sun_Mass))
        return self.p
    def keplerianToStateVectors(self):
        self.position = vec3(0,0,0)
        self.velocity = vec3(0,0,0)
        self.calcAngularMomentum()
        self.calcPeriod()
        
        # Convert to radians
        E = math.radians(self.E)
        Ω = math.radians(self.Ω)
        ω = math.radians(self.ω)
        inc = math.radians(self.i)
        v = math.radians(self.v)
        #r = self.a * (1-self.e) * math.cos(E)
        r = self.a * ((1-self.e**2) / (1+self.e * math.cos(v)))

        # Convert to rectangular
        rx = r * math.cos(v)
        ry = r * math.sin(v)
        rz = 0

        rp = math.sqrt((G*Sun_Mass)*self.a) / r
        rxd = rp * -(math.sin(E))
        ryd = rp * math.sqrt(1-self.e**2) * math.cos(E)
        rzd = 0

        #o = numpy.array([[x],[y],[z]])
        
        # Position
        x = rx * (math.cos(ω) * math.cos(Ω) - math.sin(ω) * math.cos(inc) * math.sin(Ω)) - ry * (math.sin(ω) * math.cos(Ω) + math.cos(ω) * math.cos(inc) * math.sin(Ω))
        y = rx * (math.cos(ω) * math.sin(Ω) + math.sin(ω) * math.cos(inc) * math.cos(Ω)) + ry * (math.cos(ω) * math.cos(inc) * math.cos(Ω) - math.sin(ω) * math.sin(Ω))
        z = rx * (math.sin(ω) * math.sin(inc)) + ry * (math.cos(ω) * math.sin(inc))

        # Velocity
        xdot = rxd * (math.cos(ω) * math.cos(Ω) - math.sin(ω) * math.cos(inc) * math.sin(Ω)) - ryd * (math.sin(ω) * math.cos(Ω) + math.cos(ω) * math.cos(inc) * math.sin(Ω))
        ydot = rxd * (math.cos(ω) * math.sin(Ω) + math.sin(ω) * math.cos(inc) * math.cos(Ω)) + ryd * (math.cos(ω) * math.cos(inc) * math.cos(Ω) - math.sin(ω) * math.sin(Ω))
        zdot = rxd * (math.sin(ω) * math.sin(inc)) + ryd * (math.cos(ω) * math.sin(inc))
        

        '''R3 = numpy.array([[ math.cos(Ω), -math.sin(Ω), 0 ],
                        [ math.sin(Ω), math.cos(Ω), 0 ],
                        [ 0, 0, 1 ]])
        old R2 = numpy.array([[ math.cos(inc), 0, -math.sin(inc) ],
                        [ 0, 1, 0 ],
                        [ math.sin(inc), 0, math.cos(inc) ]])   
        R2 = numpy.array([[ 1, 0, 0 ],
                        [ 0, math.cos(inc), -math.sin(inc) ],
                        [ 0, math.sin(inc), math.cos(inc) ]])                
        R1 = numpy.array([[ math.cos(ω), -math.sin(ω), 0 ],
                        [ math.sin(ω), math.cos(ω), 0 ],
                        [ 0, 0, 1 ]])
       
        Rsum = numpy.matmul(R1, o)
        Rsum2 = numpy.matmul(R2, Rsum)
        Rsum3 = numpy.matmul(R3, Rsum2)

        self.position.x = Rsum3[0,0]
        self.position.y = Rsum3[1,0]
        self.position.z = Rsum3[2,0]'''

        self.position.x = x
        self.position.y = y
        self.position.z = z
        self.velocity.x = xdot
        self.velocity.y = ydot
        self.velocity.z = zdot

        return 0
    def StateVectorsToKeplarian(self):
        return 0

        
earth = Body("Earth", 5.97237e24, 6371)
earth.loadKeplarian(1.496534962730336e8, # a
              1.711862905357640e-2, # e
              4.181344269688850e-4, # i
              1.350829426264774e2, # Ω
              3.267259945200456e2, # ω
              0) # t

os.system("cls")
print ("Epoch J2000 (JD " + str(T) + ")")
print ("Julian date: " + str(calcJD(calcJDN_GC(2000,2,1),12,0,0)))
print("Semi-major axis: " + format(earth.a,'1e') + " km")
print("Orbit: " + format(earth.calcApocenter()/1000,'3e') + " x " + format(earth.calcPericenter()/1000,'3e') + " (" + format(earth.calcApoapsis()/1000,'3e') + " x " + format(earth.calcPeriapsis()/1000,'3e') + ") km")
print("Eccentricity: " + str(earth.e))
print("Mean motion: " + str(earth.getMeanMotion()) + " degrees/s")
print ("Orbital period: " + str(earth.calcPeriod()/(60*60*24)) + " days/year")
print ("Mean Velocity: " + str(earth.calcVelocityMean()) + " m/s")
print ("Velocity at apoapsis: " + str(earth.calcVelocityApoapsis()))
print ("Velocity at periapsis: " + str(earth.calcVelocityPeriapsis()))
print ("Mean Anomaly: " + str(earth.calcMeanAnomalyFromEpoch()) + " degrees")
print ("Eccentric Anomaly: " + str(earth.calcEccentricAnomaly()) + " degrees")
print ("True Anomaly: " + str(earth.calcTrueAnomaly()) + " degrees")
earth.keplerianToStateVectors()
print ("State Vectors")
'''
print ("X: " + format(earth.position.x,'1e'))
print ("Y: " + format(earth.position.y,'1e'))
print ("Z: " + format(earth.position.z,'1e'))
print ("vX: " + format(earth.velocity.x,'1e'))
print ("vY: " + format(earth.velocity.y,'1e'))
print ("vZ: " + format(earth.velocity.z,'1e'))'''
print("Position")
print(f"X: {earth.position.x:e}")
print(f"Y: {earth.position.y:e}")
print(f"Z: {earth.position.z:e}")
print("Velocity")
print(f"X: {earth.velocity.x:e}")
print(f"Y: {earth.velocity.y:e}")
print(f"Z: {earth.velocity.z:e}")

print ("Angular momentum: " + format(earth.h,'3e'))


fig = plt.figure()
ax = plt.axes(projection = '3d')

X = []
Y = []
Z = []
Xint = [earth.position.x]
Yint = [earth.position.y]
Zint = [earth.position.z]
Xvint = [earth.velocity.x]
Yvint = [earth.velocity.y]
Zvint = [earth.velocity.z]
Xv = []
Yv = []
Zv = []
Xh = []
Yh = []
Zh = []


positions = [0,0,0]
for i in range(0,700,1):
    earth = Body("Earth", 5.97237e24, 6371)
    venus = Body("Venus", 48.685e23, 6051.84)
    halley1p = Body("Halley's Comet", 5000, 5.5)
    earth.loadKeplarian(1.496534962730336e8, # a
              1.711862905357640e-2, # e
              4.181344269688850e-4, # i
              1.350829426264774e2, # Ω
              3.267259945200456e2, # ω
              i*60*60*24) # t
    venus.loadKeplarian(1.082081681708332e8, # a
              6.755786250503024e-3, # e
              3.394589648659516, # i
              7.667837463924961e1, # Ω
              5.518596653686583e1, # ω
              i*60*60*24) # t
    halley1p.loadKeplarian(2.681019359492625e9, # a
              9.672701666828314e-1, # e
              1.621960405866950e2, # i
              5.950786974713448e1, # Ω
              1.124496743913674e2, # ω
              i*60*60*24) # t
    earth.calcMeanAnomalyFromEpoch()
    venus.calcMeanAnomalyFromEpoch()
    halley1p.calcMeanAnomalyFromEpoch()
    earth.calcEccentricAnomaly()
    venus.calcEccentricAnomaly()
    halley1p.calcEccentricAnomaly()
    earth.calcTrueAnomaly()
    venus.calcTrueAnomaly()
    halley1p.calcTrueAnomaly()
    earth.keplerianToStateVectors()
    venus.keplerianToStateVectors()
    halley1p.keplerianToStateVectors()
    X.append(earth.position.x)
    Y.append(earth.position.y)
    Z.append(earth.position.z)
    Xv.append(venus.position.x)
    Yv.append(venus.position.y)
    Zv.append(venus.position.z)
    Xh.append(halley1p.position.x)
    Yh.append(halley1p.position.y)
    Zh.append(halley1p.position.z)
plt.plot(X, Y, Z)
plt.plot(Xv, Yv, Zv)
plt.plot(Xh, Yh, Zh)
plt.plot([0],[0],[0],'yo')
plt.plot(Xint,Yint,Zint,'bo')
plt.quiver(Xint,Yint, Zint,Xvint,Yvint, Zvint, length=2e6, color='red')
plt.axis([-1.5e11,1.5e11,-1.5e11,1.5e11])
ax.set_zlim(-1.5e11,1.5e11)
plt.show()