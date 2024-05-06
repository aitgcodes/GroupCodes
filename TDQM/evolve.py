import sys
import re
import numpy as np
from numpy import linalg as LA
import cmath
from fields import Laser
from qsystem import Hamiltonian
from propagators import Propagator

## Read in Hamiltonian and Perturbation matrix
fname = sys.argv[1]
fp=open(fname,"r")
lines=fp.readlines()
fp.close()

nst = int(lines[0].strip().split()[0])
print("No. of states :", nst)

diags0 = np.zeros(nst, dtype=float)
diags = np.zeros(nst, dtype=float)
hpert = np.zeros((nst,nst), dtype=float)
ziden = np.complex(1.0,0.0)

arr=lines[1].strip().split()
if len(arr) < nst:
        print("Not enough entries for H0!")
        quit()
for i in range(nst):
        diags0[i] = float(arr[i])

arr=lines[2].strip().split()
if len(arr) < nst:
        print("Not enough entries for diagonals of perturbation!")
        quit()
for i in range(nst):
        diags[i] = float(arr[i])

i = 1
for line in lines[3:nst+2]:
        arr = line.strip().split()
        if len(arr) < i:
                print("Not enough entries for lower off-diagonal of perturbation!")
                quit()
        for j in range(i):
	        hpert[i,j] = float(arr[j])
        i += 1

hpert += hpert.T
hpert += np.diag(diags)

print(hpert)

ham = Hamiltonian(nst,diags0,hpert) ## Initialize the Hamiltonian object

#Read in laser and simulation details

fname = sys.argv[2]
fp1 = open(fname,"r")
lines = fp1.readlines()
fp1.close()

## Set defaults
f0 = 0.0
w0 = 0.0
fpol = [0.0, 0.0, 0.0]
t0 = 0.0
tlas = False
envlp = "None"
tsig = 0.0
nt = 1
delt = 0.01
swexp = "power"
swprop = "etrs"

for line in lines:

        tdelt = re.search("DELT",line,re.I)
        tnt = re.search("NT",line,re.I)
        tw0 = re.search("W0",line, re.I) 
        tf0 = re.search("EAMP",line, re.I) 
        tt0 = re.search("T0",line, re.I) 
        ttsig = re.search("TSIG",line, re.I) 
        tlaser = re.search("TLAS",line, re.I)
        tfpol = re.search("FPOL",line,re.I)
        texp = re.search("EXP_METHOD",line,re.I)
        tprop = re.search("PROP",line,re.I)
        tenvlp = re.search("ENVLP",line,re.I)

        if tdelt:
            arr = line[tdelt.start():].strip().split()
            delt = float(arr[1])
        if tnt:
            arr = line[tnt.start():].strip().split()
            nt = int(arr[1])
        if tw0:
            arr = line[tw0.start():].strip().split()
            w0 = float(arr[1])
        if tf0:
            arr = line[tf0.start():].strip().split()
            f0 = float(arr[1])
        if tt0:
            arr = line[tt0.start():].strip().split()
            t0 = float(arr[1])
        if ttsig:
            arr = line[ttsig.start():].strip().split()
            tsig = float(arr[1])
        if tfpol:
            arr = line[tfpol.start():].strip().split()
            for k in range(3):
                fpol[k] = float(arr[k+1])
        if texp:
            arr = line[texp.start():].strip().split()
            swexp = arr[1]	
        if tprop:
            arr = line[tprop.start():].strip().split()
            swprop = arr[1]	
        if tlaser:
            arr = line[tlaser.start():].strip().split()
            tlas = arr[1]
        if tenvlp:
            arr = line[tenvlp.start():].strip().split()
            envlp = arr[1]

if tlas:
        laser = Laser(tlas,f0,w0,fpol)
        laser.InitEnvelope(envlp,[t0,tsig]) ### Initiate the laser

tin = 0.0
prop = Propagator(tin,nt,delt,swexp,swprop) ### Initiate the propagator

if tlas:
        print("Laser turned on with ...")
        print("	Electric Field Amplitude (a.u.): ",f0)
        print("	Central frequency of laser (a.u.): ",w0)
        print("	Laser polarization: ",fpol)
        flas = open("laser.dat","w")
        for it in range(nt):
                ft = laser.Amplitude(prop.time[it])
                flas.write("%10.6f %10.6f %10.6f %10.6f\n" %(prop.time[it],ft*fpol[0],ft*fpol[1],ft*fpol[2]))
        flas.close()
		
else:
        print("Not turned on!")
print("\n")
print("Propagator details ...")
print("	Intial time for simulation: ", tin)
print("	Total number of simulation steps: ", nt)
print("	Time step used for simulation (a.u.): ", delt)
print(" Exponential method used: ", swexp)
print(" Propagation method used: ", swprop)
print("\n\n")
print("Starting propagation ...")
#### Start propagation
cvec = np.zeros((nst,nt), dtype = 'complex')
amat = np.zeros((nst,nst), dtype = 'complex')
cvec[0,0] = ziden

fpt = open("evolve.dat", "w")

for it in range(nt-1):
        t = it*delt
        prop.time[it] = t
        fpt.write("%10.6f " %(t))
        for ist in range(nst):
             fpt.write("%10.6f " %(abs(cvec[ist,it])**2))
        fpt.write("\n")

        amat = prop.TDPropagate(ham,laser,it)
        cvec[:,it+1] = np.dot(amat,cvec[:,it])

it += 1
t = it*delt
prop.time[it] = t
fpt.write("%10.6f " %(t))
for ist in range(nst):
        fpt.write("%10.6f " %(abs(cvec[ist,it])**2))
fpt.write("\n")
