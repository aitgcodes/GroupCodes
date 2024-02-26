def get_energies(fdir,mult,nroots):

   if mult == "Singlets":
      srchstr1="E(RB3LYP"
   else:
      srchstr1="E(UB3LYP"

   srchstr2="Total energy after correction"

   fgs=Path(fdir/mult/"clrpcm.log")
   #print(fgs)

   eners=[]
   if fgs.exists():
      lines=fgs.read_text().splitlines()
      gslines=[line for line in lines if srchstr1 in line]
      eners.append(float(gslines[-1].split()[4]))
      #print(gsener)
   else:
      print(f"File {str(fgs)} not found!")

   for root in range(nroots):
       fex=Path(fdir/mult/f"ss_pcm_root{root+1}.log")
       if fex.exists():
          lines=fex.read_text().splitlines()
          exlines=[line for line in lines if srchstr2 in line]
          eners.append(float(exlines[-1].split()[5]))

   return eners

from pathlib import Path
import sys

com_dir=sys.argv[1]
calc_dir=sys.argv[2]

com_dir=Path(com_dir)
com_files = [f for f in com_dir.glob("*.com") if f.is_file()]

calc_dir=Path(calc_dir)

### Prepare the work directories for the calcualtions
angles=[]
for f in com_files:
   fname=f.stem
   angles.append(int(fname.rstrip(fname[-1])))

### Prepare the order of calculations
angles.sort()    
labels=[str(f)+"d" for f in angles]
multiplicities=[{"name": "Singlets", "M": 1, "NStates": 3, "Roots": 3}, 
                {"name": "Triplets", "M": 3, "NStates": 3, "Roots": 2}]

### Loop over labels and multipliticies and perform the spectrum workflow
pec={"Singlets": [], "Triplets": []}

for label in labels:
   for mult in multiplicities:
      fcom=Path(calc_dir/label)  
      energies=get_energies(fcom,mult["name"],mult["Roots"])
      #print(energies)
      if mult["name"] in pec.keys():
         pec[mult["name"]].append(energies)    

eref=pec["Singlets"][0][0]
au2eV=27.2114
#print(eref)

for mlt in pec.keys():
    fp=open(mlt+"_pec.dat","w")
    for (angle,enr) in zip(angles,pec[mlt]):
        enr=[(e-eref)*au2eV for e in enr]
        alst = [str(x) for x in [angle]+enr]
        fp.write("  ".join(alst)+"\n")
    fp.close()
