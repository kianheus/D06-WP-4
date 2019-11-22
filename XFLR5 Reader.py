"""
XFLR5 data reader

This program reads text files from XFLR5
"""
import scipy as sp
from math import *
from scipy import interpolate
from scipy import integrate
from matplotlib import pyplot as plt

#--------------------Read files from XFLR5 for basic data---------------------
#Read the XFLR file for a=0
XFLRfile0 = open("MainWing_a=0_Useful.dat","r")
lines0=XFLRfile0.readlines()
XFLRfile0.close()

#Read the XFLR file for a=10
XFLRfile10 = open("MainWing_a=10_Useful.dat","r")
lines10=XFLRfile10.readlines()
XFLRfile10.close()

#Create lists for variables, for a=0 and a=10
yspanlst        = []
chordlst        = []
InducedAoA0lst  = []
Cl0lst          = []
InducedCd0lst   = []
CmQuarter0lst   = []

#Chord and yspan only defined once
InducedAoA10lst = []
Cl10lst         = []
InducedCd10lst  = []
CmQuarter10lst  = []

#----------------------------For alpha=0-------------------------------
#This loop reads the file and creates the tables
for line in lines0: 

    line = line.strip()

    #Save total CL of wing, and Induced CD
    if line.count("CL    =")==1:
        CL0  = float(line.strip("CL    ="))
    if line.count("ICd   =")==1:
        ICD0 = float(line.strip("ICd   ="))
    if line.count("Cm   =")==1:
        Cm0  = float(line.strip("Cm   ="))
        
    #Protect against empty line and use only columns with numerical (positive or negative) values
    if len(line)>0 and (line[0].isdigit()==True or line[0]=="-"):

        #Separate into columns
        columns = line.split()
        
        #Interpret columns
        yspanlst.append(float(columns[0].strip()))
        chordlst.append(float(columns[1].strip()))
        InducedAoA0lst.append(float(columns[2].strip()))
        Cl0lst.append(float(columns[3].strip()))
        InducedCd0lst.append(float(columns[5].strip()))
        CmQuarter0lst.append(float(columns[7].strip()))

#Generate function for each of the variables with SciPy, as a function of y
chord       = sp.interpolate.interp1d(yspanlst,chordlst,kind='cubic',fill_value="extrapolate")
InducedAoA0 = sp.interpolate.interp1d(yspanlst,InducedAoA0lst,kind='cubic',fill_value="extrapolate")
Cl0         = sp.interpolate.interp1d(yspanlst,Cl0lst,kind='cubic',fill_value="extrapolate")
InducedCd0  = sp.interpolate.interp1d(yspanlst,InducedCd0lst,kind='cubic',fill_value="extrapolate")
CmQuarter0  = sp.interpolate.interp1d(yspanlst,CmQuarter0lst,kind='cubic',fill_value="extrapolate")

#-------------------------For alpha=10---------------------
#This loop reads the file and creates the tables
for line in lines10: 

    line = line.strip()

    #Save total CL of wing and Induced CD
    if line.count("CL    =")==1:
        CL10  = float(line.strip("CL    ="))
    if line.count("ICd   =")==1:
        ICD10 = float(line.strip("ICd   ="))
    if line.count("Cm   =")==1:
        Cm10  = float(line.strip("Cm   ="))
        
        
    #Protect against empty line and use only columns with numerical (positive or negative) values
    if len(line)>0 and (line[0].isdigit()==True or line[0]=="-"):

        #Separate into columns
        columns = line.split()
        
        #Interpret columns
        InducedAoA10lst.append(float(columns[2].strip()))
        Cl10lst.append(float(columns[3].strip()))
        InducedCd10lst.append(float(columns[5].strip()))
        CmQuarter10lst.append(float(columns[7].strip()))

#Generate function for each of the variables with SciPy, as a function of y
InducedAoA10 = sp.interpolate.interp1d(yspanlst,InducedAoA10lst,kind='cubic',fill_value="extrapolate")
Cl10         = sp.interpolate.interp1d(yspanlst,Cl10lst,kind='cubic',fill_value="extrapolate")
InducedCd10  = sp.interpolate.interp1d(yspanlst,InducedCd10lst,kind='cubic',fill_value="extrapolate")
CmQuarter10  = sp.interpolate.interp1d(yspanlst,CmQuarter10lst,kind='cubic',fill_value="extrapolate")



#-------------------------Functions for forces-------------------------

#Define flow conditions and required total CLd
density = 0.3796    #[kg/m^3]
V       = 228.31    #[m/s]
CLd     = 0.4545    #[-]

#Required AoA for a certain total CLd
AoAd = asin((CLd-CL0)/(CL10-CL0)*sin(10*pi/180))

#--Function for CL(y) distribution for a required (known) total CLd--
def CL(CLd,y):
    return Cl0(y) + (CLd-CL0)/(CL10-CL0)*(Cl10(y)-Cl0(y))

#Induced angle of attack as a function of the AoA
def InducedAoA(AoAd,y):
    return InducedAoA0(y) + (AoAd*180/pi-0)/(10-0)*(InducedAoA10(y)-InducedAoA0(y))

#D'= L' * inducedAoA
def InducedCd(CLd,y,AoA):
    return abs(InducedAoA(AoAd,y))*pi/180* CL(CLd,y)

def CmQuarter(AoA,y):
    return CmQuarter0(y) + (AoA*180/pi-0)/(10-0)*(CmQuarter10(y)-CmQuarter0(y))

#-------------------------Distribution-------------------------------
#Lift per unit span: L'
def L(y):
    #L'(y) = C_L(y) * q * c(y)
    L = CL(CLd,y)*0.5*density*V**2*chord(y)
    return L

#Induced drag per unit span: D'_i
def D(y):
    #D'(y) = C_D(y) * q * c(y)
    D = InducedCd(CLd,y,AoAd)*0.5*density*V**2*chord(y)
    return D

#Pitching moment:
def M(y):
    #M = Cm * q * c^2
    M = CmQuarter(AoAd,y)*0.5*density*V**2*(chord(y))**2
    return M

#----------------------------Distribution of forces-----------------
ylst              = []
Liftlst           = []
Draglst           = []
PitchingMomentlst = []
InducedAoA0lst    = []
InducedAoA10lst   = []


#Plot the L, D and M at 125 distributed along the half span:
for i in range(0,1000):
    #0.096 metres per point, and 1200 points calculated per half wing
    y =i*12/1000
    #Create lists for lift, induced drag, pitching moment at each y point
    ylst.append(y)
    Liftlst.append(L(y))
    Draglst.append(D(y))
    PitchingMomentlst.append(M(y))
    InducedAoA0lst.append(InducedAoA0(y))
    InducedAoA10lst.append(InducedAoA10(y))



#-------------------------Total forces-----------------------

#Total lift, drag, and pitching moment
LiftTotal = sp.integrate.quad(L, -12.005, 12.005)[0]
DragTotal = sp.integrate.quad(D, 0, 12.005)[0]+ sp.integrate.quad(D, -12.005, 0)[0]
PitchingMomentTotal = sp.integrate.quad(M, -12.005, 12.005)[0]

print("Total lift",round(LiftTotal,4))
print("Total drag:",round(DragTotal,4))
print("Total pitching moment:",round(PitchingMomentTotal,4))



print("Ready")
