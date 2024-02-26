import sys
import pathlib as pl
import subprocess as subp
from g16io import *

fname=sys.argv[1]
fname=pl.Path(fname)
nprocs=10
d=read_com(fname)
geom=d["geometry"]
wdir=pl.Path("../test")
fcom=pl.Path(wdir/"trial.com")
write_com(fcom,d)
fout=pl.Path(wdir/"trial.log")

com_files = [f for f in fname.parent.glob("*.com") if f.is_file()]
for f in com_files:
    print(f.stem)
quit()

def launchgauss(workdir,inpfile,outfile,nprocs):
    1
    cmdline=f"cd {workdir}; g16 < {inpfile} >& {outfile}"
    lp=subp.run(cmdline,shell=True,capture_output=True,text=True,check=True)
    if lp.return_code:
       print(lp.stdout)
       print(lp.stderr)
       raise Expection("Gaussian job not run!")
    else:
       return (lp.stdout, lp.stderr)

(lstd,lserr)=launchgauss(wdir,fcom.name,fout.name,nprocs)
