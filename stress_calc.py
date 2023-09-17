import math
import numpy as np
import matplotlib.pyplot as plt

# Pipe constants
m_sch160 = 7.25
m_sch80s = 5.41
m_xs = 5.41
m_sch40s = 4.05

# Material class is defined as a class with attributes olet class, pipe schedule, and nipple schedule
class materialClass:
    def __init__(self, oletClass, pipeSch, nipSch, is9SA0E): # 9SA0E has a different valve mass
        # All length data is in mm
        # Olet
        if oletClass == 3000:
            self.oletLength = 34.925
            self.oletMass = 1
            # Olet dimensions in mm, refer to bonney forge sockolet specs
            self.oletA = (1 + 3/8)*25.4
            self.oletB = (2 + 7/8)*25.4
            self.oletC = 2*25.4
        elif oletClass == 6000:
            self.oletLength = 42.8625
            self.oletMass = 2
            self.oletA = (1 + 11/16)*25.4
            self.oletB = (3 + 1/4)*25.4
            self.oletC = (1 + 15/16)*25.4

        # Nipple
        self.nipLength = 148 - self.oletLength
        if nipSch == "SCH80S":
            self.nipMass = self.nipLength * m_sch80s/1000
            self.nipID = 38.14
            self.nipOD = 48.3
        elif nipSch == "SCH160":
            self.nipMass = self.nipLength * m_sch160/1000
            self.nipID = 34.02
            self.nipOD = 48.3

        # Valve
        self.valveLength = 90
        if is9SA0E:
            self.valveMass = 10.87
        else:
            self.valveMass = 6.2
        
        # Pipe
        self.pipeLength = 151
        self.pipeOD = 48.3
        if pipeSch == 'XS' or pipeSch == "SCH80S":
            self.pipeMass = self.pipeLength * m_sch80s/1000
            self.pipeID = 38.14
        elif pipeSch == "SCH160":
            self.pipeMass = self.pipeLength * m_sch160/1000
            self.pipeID = 34.02
        elif pipeSch == "SCH40S":
            self.pipeMass = self.pipeLength * m_sch40s/1000
            self.pipeID =  40.94

        # Total mass and length of SBC
        self.totalMass = self.oletMass + self.nipMass + self.valveMass + self.pipeMass
        self.totalLength = self.oletLength + self.nipLength + self.valveLength + self.pipeLength

        # Center of Gravity location, assuming the cog of each component is located right in the middle of its length, units in mm
        self.cogOlet = self.oletLength/2
        self.cogNipple = self.oletLength + self.nipLength/2
        self.cogValve = self.oletLength + self.nipLength + self.valveLength/2
        self.cogPipe = self.oletLength + self.nipLength + self.valveLength + self.pipeLength/2
        self.overallCog = (self.cogOlet*self.oletMass + self.cogNipple*self.nipMass + self.cogValve*self.valveMass + self.cogPipe*self.pipeMass)/(self.oletMass + self.nipMass + self.valveMass + self.pipeMass)

        # Acceleration, force, and moment initialization in m/s^2, N, and N mm respectively
        self.oletAcc = 0
        self.nipAcc = 0
        self.valveAcc = 0
        self.pipeAcc = 0
        self.cogAcc = 0

        self.oletForce = 0
        self.nipForce = 0
        self.valveForce = 0
        self.pipeForce = 0

        self.oletMoment = 0
        self.nipMoment = 0
        self.valveMoment = 0
        self.pipeMoment = 0

# Function for calculating acceleration based on vibration frequency declaration
def accelFromAmp(freq):
    freq = float(freq)
    exp = (math.log10(freq) + 1.871083)/2.084547
    amp = (pow(10, exp)*math.sqrt(2)/(2*math.pi*freq))/1000
    accel = amp * (2*math.pi*freq) ** 2
    return accel

# ======== Start of analysis ===========
# Material and constants declaration (Sockolet class, Pipe Schedule, Nipple Schedule, is9SA0E)
material = materialClass(3000, "SCH40S", "SCH80S", False)
vibFreq = input("Vibration frequency in Hz: ")
inspectAcc = accelFromAmp(vibFreq)
stressLimit = 35
weldPositions = [0, material.oletA, 148, 238] # Location of weldings in SBC, measured from the root of SBC
# Upper and lower position initialization for bisection method
x_lower = 1
x_upper = material.totalLength

# Function declarations
# Acceleration interpolator
def interpolate(x_want, x_known, y_known):
    # if x_known != 0:
    y_want = y_known * x_want / x_known
    return y_want

# Area moment of inertia function:
def areaInertia(od, id):
    I = math.pi*(od**4 - id**4)/64
    return I

# Shear stress calculation
def shearStress(V, M, r, I):
    h = 7 # h is assumed to be 7 mm
    A = 1.414 * math.pi * h * r
    shear1 = V/A
    # shear2 = M * r/I
    shear = 6*(shear1 ** 2)
    return shear

# Bending stress calculation
def bendStress(m):
    oletRoot = m[0] * material.oletB/(2*areaInertia(material.oletB, material.nipID))
    oletTip = m[1] * material.nipOD/(2*areaInertia(material.nipOD, material.nipID))
    nipValve = m[2] * material.nipOD/(2*areaInertia(material.nipOD, material.nipID))
    valvePipe = m[3] * material.pipeOD/(2*areaInertia(material.pipeOD, material.pipeID))
    return [oletRoot, oletTip, nipValve, valvePipe]

# Internal force plotting, requires the reaction forces to be calculated previously
def intForcePlot():
    section = np.linspace(0, material.totalLength, int(material.totalLength))
    intForceArr = [0 for i in range(len(section))]
    intMomentArr = [0 for i in range(len(section))]
    for n, i in enumerate(section):
        intForceArr[n], intMomentArr[n] = intForce(i, reactForce, reactMoment)

    plt.figure("Fig 1")
    plt.plot(section, intForceArr)
    plt.title("Internal Shear Force Diagram")
    plt.xlabel("Position")
    plt.ylabel("Shear Force (N)")
    plt.grid()
    plt.figure("Fig 2")
    plt.plot(section, intMomentArr)
    plt.title("Internal Moment Diagram")
    plt.xlabel("Position")
    plt.ylabel("Moment (N mm)")
    plt.grid()
    plt.show()

# Reaction force and moment calculation
def intForceProcedure(insPos):
    material.oletAcc = interpolate(material.cogOlet, insPos, inspectAcc)
    material.nipAcc = interpolate(material.cogNipple, insPos, inspectAcc)
    material.valveAcc = interpolate(material.cogValve, insPos, inspectAcc)
    material.pipeAcc = interpolate(material.cogPipe, insPos, inspectAcc)
    material.cogAcc = interpolate(material.overallCog, insPos, inspectAcc)

    # Reaction force calculation
    material.oletForce = material.oletMass * material.oletAcc
    material.nipForce = material.nipMass * material.nipAcc
    material.valveForce = material.valveMass * material.valveAcc
    material.pipeForce = material.pipeMass * material.pipeAcc
    reactForce = material.oletForce + material.nipForce + material.valveForce + material.pipeForce # in N

    # Reaction moment calculation
    material.oletMoment = material.oletMass * material.oletAcc * material.cogOlet
    material.nipMoment = material.nipMass * material.nipAcc * material.cogNipple
    material.valveMoment = material.valveMass * material.valveAcc * material.cogValve
    material.pipeMoment = material.pipeMass * material.pipeAcc * material.cogPipe
    reactMoment = material.oletMoment + material.nipMoment + material.valveMoment + material.pipeMoment # in N mm
    return reactForce, reactMoment

# Internal force calculation function
def intForce(pos, reactForce, reactMoment):
    intForce = 0
    intMoment = 0
    if pos < material.cogOlet:
        intForce = -1 * reactForce
        intMoment = reactMoment - reactForce*pos
    elif pos >= material.cogOlet and pos < material.cogNipple:
        intForce = -1 * (reactForce - material.oletForce)
        intMoment = reactMoment + material.oletForce*(pos - material.cogOlet) - reactForce*pos
    elif pos >= material.cogNipple and pos < material.cogValve:
        intForce = -1 * (reactForce - material.oletForce - material.nipForce)
        intMoment = reactMoment + material.oletForce*(pos-material.cogOlet) + material.nipForce*(pos-material.cogNipple) - reactForce*pos
    elif pos >= material.cogValve and pos < material.cogPipe:
        intForce = -1 * (reactForce - material.oletForce - material.nipForce - material.valveForce)
        intMoment = reactMoment + material.oletForce*(pos-material.cogOlet) + material.nipForce*(pos-material.cogNipple) + material.valveForce*(pos-material.cogValve) - reactForce*pos
    elif pos >= material.cogPipe and pos <= material.totalLength:
        intForce = 0
        intMoment = reactMoment + material.oletForce*(pos-material.cogOlet) + material.nipForce*(pos-material.cogNipple) + material.valveForce*(pos-material.cogValve) + material.pipeForce*(pos-material.cogPipe) - reactForce*pos
    else:
        print("Invalid length")
        return
    return intForce, intMoment

def vonmises(bend, shear):
    vonStress = 0.707 * math.sqrt(bend ** 2 + shear)
    return vonStress

# ====== Begin iteration =========
isFound = False
iternum = 1
while isFound == False: # Repeat until the optimal inspection point is found
    x_r = x_lower/2 + x_upper/2 # Midpoint of x_lower and x_upper

    # Reset array values at each iteration
    v_lower = []
    v_upper = []
    m_lower = []
    m_upper = []
    stressLower = []
    stressUpper = []

    # Calculate internal forces due to displacement at x_lower and x_r, respectively.
    for i in range(2): 
        if i == 0:
            reactForce, reactMoment = intForceProcedure(x_lower)
            for i in weldPositions:
                intV, intM = intForce(i, reactForce, reactMoment)
                v_lower.append(abs(intV))
                m_lower.append(intM)
        else:
            reactForce, reactMoment = intForceProcedure(x_r)
            for i in weldPositions:
                intV, intM = intForce(i, reactForce, reactMoment)
                v_upper.append(abs(intV))
                m_upper.append(intM)

    # Bending stress calculation
    stressLower = bendStress(m_lower)
    stressUpper = bendStress(m_upper)

    # Von Mises Stress calculation
    stressLower[0] = vonmises(stressLower[0], shearStress(v_lower[0], m_lower[0], material.oletB, areaInertia(material.oletB, material.nipID)))
    stressUpper[0] = vonmises(stressUpper[0], shearStress(v_upper[0], m_upper[0], material.oletB, areaInertia(material.oletB, material.nipID)))
    stressLower[1] = vonmises(stressLower[1], shearStress(v_lower[1], m_lower[1], material.nipOD, areaInertia(material.nipOD, material.nipID)))
    stressUpper[1] = vonmises(stressUpper[1], shearStress(v_upper[1], m_upper[1], material.nipOD, areaInertia(material.nipOD, material.nipID)))
    stressLower[2] = vonmises(stressLower[2], shearStress(v_lower[2], m_lower[2], material.nipOD, areaInertia(material.nipOD, material.nipID)))
    stressUpper[2] = vonmises(stressUpper[2], shearStress(v_upper[2], m_upper[2], material.nipOD, areaInertia(material.nipOD, material.nipID)))
    stressLower[3] = vonmises(stressLower[3], shearStress(v_lower[3], m_lower[3], material.pipeOD, areaInertia(material.pipeOD, material.pipeID)))
    stressUpper[3] = vonmises(stressUpper[3], shearStress(v_upper[3], m_upper[3], material.pipeOD, areaInertia(material.pipeOD, material.pipeID)))

    # Stress difference calculation for each upper and lower bound
    stressDiffLower = max(stressLower) - 35
    stressDiffUpper = max(stressUpper) - 35
    maxmaxStress = max(max(stressLower), max(stressUpper))

    # Uses bisection method to determine the optimal inspection location
    if stressDiffLower*stressDiffUpper < 0 and round(stressDiffLower*stressDiffUpper) != 0: # Optimal inspection location is between x_lower and x_r
        x_upper = x_r
    elif stressDiffLower*stressDiffUpper > 0 and round(stressDiffLower*stressDiffUpper) != 0: # Optimal inspection location is between x_r and x_upper
        x_lower = x_r
    elif round(stressDiffLower*stressDiffUpper) == 0 or abs(stressDiffLower*stressDiffUpper) < 0.0000001: # Optimal inspection location found
        isFound = True
        inspect_loc = x_r
    
    if maxmaxStress < 0.1 and x_r > 388:
        print("Inspection point not found because maximum stress is too small. SBC is safe in this frequency")
        break
    elif x_r == 389.0 and x_lower == 389.0:
        print("Optimal inspection point is beyond SBC length. SBC will most likely fail")
        break

    iternum = iternum + 1
    
if isFound:
    print(f"Optimal inspect location for vibration frequency {vibFreq} Hz: {inspect_loc} mm. Number of iterations: {iternum}")
    intForcePlot()
