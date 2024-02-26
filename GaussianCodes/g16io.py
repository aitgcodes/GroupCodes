def read_com(fin):
    
    #fpin=open(fin,"r")
    #lines=fpin.readlines()
    if fin.exists():
        lines=fin.read_text().splitlines()

    link_lines=[]
    route_lines=[]
    geometry=[]
    title=[]

    i = 0
    ### Read in link and route lines
    done = False
    while i <= len(lines) and not(done):
         line=lines[i].strip()
         if line.startswith("%"):
            link_lines.append(line.strip())
         elif line.startswith("#"):
            #print(line)
            while line:
                 print("yep",line)
                 route_lines.append(line)
                 i+=1 
                 line=lines[i].strip()
            done=True
         i+=1
    
    ### Read in title lines
    line=lines[i].strip()
    while line:
         title.append(line)
         i+=1
         line=lines[i].strip()
    i+=1
    ### Read in charge and multiplicity
    config=lines[i].strip()
    i+=1

    ### Read in geometry
    line=lines[i].strip()
    while line:
        geometry.append(line)
        i+=1
        line=lines[i].strip()

    inpdict={"link_lines": link_lines,
             "route_lines": route_lines,
             "title": title,
             "config": config,
             "geometry": geometry
            }

    return inpdict

def build_com(input_dict):
 
    lines=[]
    for line in input_dict["link_lines"]:
       lines.append(line)
    for line in input_dict["route_lines"]:
       lines.append(line)
    lines.append("")
    for line in input_dict["title"]:
       lines.append(line)
    lines.append("")
    lines.append(input_dict["config"])
    #lines.append("\n")
    for line in input_dict["geometry"]:
        lines.append(line)
    lines.append("")

    return(lines)

def write_com(fout,com_lines):

    lines=build_com(com_lines)
    fout.write_text("\n".join(lines)+"\n")
    #fp=open(fout,"w")
    #for line in com_lines:
    #    fp.write(line+"\n")
    #fp.close()
