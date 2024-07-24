import numpy as np
import sys
# from orderparams_latest import *

########################################################################

def choose_stack(arr,uax,j,k):
### uax,j,k expected in cyclic order
### indices of grid array assumed to be [z,x,y]
    stack=[]
    idxs=[]
    (nz,nx,ny) = arr[0].shape
    if uax == 0:
       idx=[(k,i,j) for i in range(nx)]
    elif uax == 1:
       idx=[(j,k,i) for i in range(ny)]
    else:
       idx=[(i,j,k) for i in range(nz)]

    for t in range(len(arr)):
        stack.append(np.array([arr[t][id] for id in idx]))
       #stack.append(arr[t][k,:,j])
    return stack

########################################################################################
def get_dihedral_from_home(arr,ref):     # [0,1,2,3]

    N=len(arr)
    narr=np.zeros(N)
    for i in range(1,N):
#        narr[i]=narr[i-1]+arr[(ref+i-1)%N]
        narr[i]=narr[i-1]+arr[(ref+i)%N]
    return narr
##############################################################################
#s1=choose_stack_a(B,0,0)
#print(s1[1])
#print(get_dihedral_from_home(s1[1],1))
######################################################################################

def get_domains_from_dihedral(arr_t, alpha, iref, N):
    import copy

    get_indices = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
    get_domain_sizes = lambda a, b: [j-i-1 for (i,j) in zip(a,b)]

    v=get_dihedral_from_home(arr_t,iref)
    v=np.cos(np.deg2rad(v))
    V=copy.deepcopy(v)
    V[1:] = V[1:]-np.cos(np.deg2rad(alpha))

    ## First get domain sizes of "up" dipoles
    V = np.sign(V)
    db = get_indices(-1,V)
    if len(db) == 0:
       dlens = [N]      
       dstarts = [0]
    else:
       db += [db[0]+N]   #### Periodic image of the first domain boundary
       dlens = get_domain_sizes(db[:-1],db[1:])
       dstarts = [(x+1)%N for x in db[:-1]] ### Starting site of domain

    ## Next get domain sizes of "down" dipoles
    db = get_indices(1,V)
    if len(db) == 0:
       dlens = [N]      
       dstarts = [0]
    else:
       db += [db[0]+N]   #### Periodic image of the first domain boundary
       dlens += get_domain_sizes(db[:-1],db[1:]) 
       dstarts += [(x+1)%N for x in db[:-1]] ### Starting site of domain

    domains = [(x,y) for (x,y) in zip(dstarts,dlens) if y > 0]

    return domains

def polar_op(arr,alpha):    # for a single chain/stack taking input [ [d1,d2,d3,d4] at t1, [d1,d2,d3,d4] at t2, ...]

    nt=len(arr) ## Total number of time steps
    N=len(arr[0]) ## Total number of sites along each direction

    # Array which keeps track of persistence times for each domain size
    # Updated when a particular domain changes (either start site or size)
    all_tp=[[] for i in range(N)]

    iref = 0
    ### domain size/length = no. of sites in a polar domain

    # Data corresponding to domains with the corresponding start site
    prev_length = np.zeros(N,dtype=int) ## domain sizes of previous frame
    tp = np.zeros(N,dtype=int) ## persistence time (number of steps domain persists)
  
    for arr_t in arr:

        domains = get_domains_from_dihedral(arr_t, alpha, iref, N)
        ### Only consider domains with non-zero length, i.e. at least one "up"
        ### dipole sandwiched between two "down" dipoles, or vice versa.

        sites=[i for i in range(N)]
        for (iref,dlen) in domains:
            if prev_length[iref] == dlen: ## domain with same starting site and size in last frame
               tp[iref] += 1
            elif tp[iref] > 0:       ## domain with same starting site but different size in last frame, resizing of domain with same leading site
               all_tp[prev_length[iref]-1].append(tp[iref])
               tp[iref] = 1
            else:                   ## domains that did not start from this site before
               tp[iref] = 1

            prev_length[iref] = dlen
            sites.remove(iref)

        ### Count sites that are no longer leading a domain they did before
        for iref in sites:
            if tp[iref] > 0:
               all_tp[prev_length[iref]-1].append(tp[iref])
               tp[iref] = 0
               prev_length[iref] = 0

    # for iref in range(N):
    #     all_tp[prev_length[iref]-1].append(tp[iref])
    # # print(tp)
    # print(tp, prev_length)

    return all_tp 

if __name__ == "__main__":

#######################################################################
### Read in the data and create a list of grids #######################
   A=np.loadtxt(sys.argv[1]) # Dihedral angle data at all frames with respect to nearest neighbor
   uax=int(sys.argv[2]) # unique axis
   data = [] # List of arrays (n x n x n) containing dihedral angle with respect to the nearest neighbor with each array denoting a frame
   (nt,nc)=A.shape
   n=round((nc-1)**(1./3)) #### Assuming the grid is cubic

   for t in range(len(A)):
       data.append(np.reshape(A[t,1:],(n,n,n)))

   OP=[[] for i in range(n)]

   for j in range(n):
       for k in range(n):
         s=choose_stack(data,uax,j,k) # List of lists, where each list contains dihedrals for a chain indexed by j,k for a specific frame
         OP_stack=polar_op(s,90.0)
         for i in range(n):
             OP[i] += OP_stack[i]

   # Write OP to file
   if len(sys.argv) > 3:
       import csv
       outfile = sys.argv[3]
       with open(outfile, "w") as f:
           wr = csv.writer(f, delimiter=" ")
           wr.writerows(OP)

   ### OP_stack will be a list of lists.
   ### The outer list is indexed by domain size (>0)
   ### and the inner list (for each size) has a collection of persistence times
   ### OP will be a similar list but over all stacks in a given plane.
   ### One can get, for instance, a distribution of persistence times for a given
   ### domain size. 
   for i in range(n):
       if len(OP[i]) != 0:
           print("Maximum and average persistence times for domain with {0:2d} sites are {1:4d} and {2:10.2f}, respectively".format(i+1,max(OP[i]),sum(OP[i])/float(len(OP[i]))))
       else:
           print("Domain with {0:2d} sites is empty".format(i+1))

   # from matplotlib import pyplot as plt

   # j = 3
   # fig, ax = plt.subplots(figsize =(10, 7))
   # lmin = 0.0
   # lmax = max(OP[j])
   # dl = (lmax-lmin)/float(10)
   # ax.hist(OP[j], bins = [lmin+i*dl for i in range(10)])

   # # Show plot
   # plt.show()
