def get_energies(fdir,mult,nroots)

   if mult == "Singlets":
      srchstr1="E(RB3LYP"
   else:
      srchstr1="E(UB3LYP"

   srchstr2="Total energy after correction"

   fgs=Path(fdir/mult/"clrpcm.log")
   #print(fgs)

   if fgs.exists():
      lines=fgs.read_text().splitlines()
      gslines=[line for line in lines if srchstr1 in line]
      gsener=float(gslines[-1].split()[4])
      #print(gsener)
   else:
      print("File not found!")

   exeners=[]
   for root in nroots:
       fex=Path(fdir/mult/f"sspcm_root{root+1}.log")
       if fex.exists():
          lines=[fex.read_text().splitlines()
          exlines=[line for line in lines if srchstr2 in line]
          ex=float(exlines[-1].split()[4])
          exeners.append(ex)

   return (gsener,exeners)
