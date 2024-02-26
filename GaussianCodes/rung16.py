from pathlib import Path
import sys
import subprocess as subp
from g16io import read_com, write_com

com_dir=sys.argv[1]
calc_dir=sys.argv[2]
if len(sys.argv) > 3:
   nprocs=sys.argv[3]
else:
   nprocs=8


com_dir=Path(com_dir)
com_files = [f for f in com_dir.glob("*.com") if f.is_file()]

calc_dir=Path(calc_dir)

def get_geometry(label):

   infile=label+".com"
   infile=Path(com_dir/infile)
   inpdict=read_com(infile)
   
   return inpdict["geometry"]

def fill_inplines(label,mode,nprocs,nst,mult,geometry=[]):

   link_lines=["%mem=8GB"]
   link_lines.append("%NProcShared="+str(nprocs))
   if mode == 0:
     chkfname=label+".chk"
     link_lines.append("%chk="+chkfname)
     route_lines=["# B3LYP/6-31+G(d,p) SCRF=(Solvent=Acetonitrile,NonEquilibrium=Save) Guess=Mix"]
     title=["TC1 in ACN under TICT ground state LRPCM"]
   else:
     oldchkfname=label+".chk"
     chkfname=label+"_ss"+str(mode)+".chk"
     link_lines.append("%oldchk="+oldchkfname)
     link_lines.append("%chk="+chkfname)
     route_lines=[f"# B3LYP/6-31+G(d,p) TD(NStates={nst},Root={mode}) Geom=Check", "  SCRF=(Solvent=Acetonitrile,CorrectedLR,NonEquilibrium=Read)"]
     title=[f"TC1 in ACN under TICT state-specific cLRPCM for state {mode}"]
     
   if mult == 1:
     config="0 1"
   elif mult == 3:
     config="0 3"
   inpdict={
       "link_lines": link_lines,
       "route_lines": route_lines,
       "title": title,
       "config": config,
       "geometry": geometry
   }

   return inpdict

def launchgauss(workdir,inpfile,outfile,nprocs):

    cmdline=f"cd {workdir}; g16 < {inpfile} >& {outfile}"
    lp=subp.run(cmdline,shell=True,capture_output=True,text=True,check=True)
    if lp.returncode:
       raise Expection("Gaussian job not run!")
    else:
        return (lp.stdout, lp.stderr)

### Prepare the work directories for the calcualtions
angles=[]
for f in com_files:
   fname=f.stem
   angles.append(int(fname.rstrip(fname[-1])))
   #print(fname)
   pdir=Path(calc_dir/fname)
   pdir.mkdir(parents=True,exist_ok=True)
   singdir=Path(pdir/"Singlets") 
   singdir.mkdir(parents=True,exist_ok=True)
   tripdir=Path(pdir/"Triplets") 
   tripdir.mkdir(parents=True,exist_ok=True)

### Prepare the order of calculations
angles.sort()    
labels=[str(f)+"d" for f in angles]
multiplicities=[{"name": "Singlets", "M": 1, "NStates": 3, "Roots": 3}, 
                {"name": "Triplets", "M": 3, "NStates": 3, "Roots": 3}]

### Loop over labels and multipliticies and perform the spectrum workflow
for label in labels:
   geom=get_geometry(label)
   for mult in multiplicities:
      ### Create input file for lrpcm single-point calculation
      inpdata=fill_inplines(label,0,nprocs,mult["NStates"],mult["M"],geometry=geom)
      fcom=Path(calc_dir/label/mult["name"]/"clrpcm.com")  
      fout=Path(fcom.parent/"clrpcm.log")
      write_com(fcom,inpdata)
      wdir=fcom.parent
      (lstd,lserr)=launchgauss(wdir,fcom.name,fout.name,nprocs)
      ### Creat input files for state-specific CorrectedLR calculations
      for root in range(mult["Roots"]):
          inpdata=fill_inplines(label,root+1,nprocs,mult["NStates"],mult["M"])
          fcom=Path(calc_dir/label/mult["name"]/f"ss_pcm_root{root+1}.com")
          fout=Path(fcom.parent/f"ss_pcm_root{root+1}.log")
          wdir=fcom.parent
          write_com(fcom,inpdata)
          (lstd,lserr)=launchgauss(wdir,fcom.name,fout.name,nprocs)
