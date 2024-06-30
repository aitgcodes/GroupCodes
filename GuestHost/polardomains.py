import numpy as np
import sys
from orderparams_latest import *

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

def polar_op(arr,alpha):    # for a single chain/stack taking input [ [d1,d2,d3,d4] at t1, [d1,d2,d3,d4] at t2, ...]

    import copy

    nt=len(arr) ## Total number of time steps
    N=len(arr[0]) ## Total number of sites along each direction
    all_tp=[[] for i in range(N)]

    iref = 0
    get_indices = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
    get_domain_sizes = lambda a, b: [j-i-1 for (i,j) in zip(a,b)]
    ### domain size/length = no. of sites in a polar domain

    prev_length = np.zeros(N,dtype=int)
    tp = np.zeros(N,dtype=int)
  
    for arr_t in arr:

        v=get_dihedral_from_home(arr_t,iref)
        v=np.cos(v)
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
        ### Only consider domains with non-zero length, i.e. at least one "up"
        ### dipole sandwiched between two "down" dipoles, or vice versa.

        sites=[i for i in range(N)]
        for (iref,dlen) in domains:
            if prev_length[iref] == dlen:
               tp[iref] += 1
            elif tp[iref] > 0:       ## resizing of domain with same leading site 
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

    return all_tp 

if __name__ == "__main__":

#######################################################################
### Read in the data and create a list of grids #######################
   A=np.loadtxt(sys.argv[1])
   uax=int(sys.argv[2])
   data = []
   (nt,nc)=A.shape
   n=round((nc-1)**(1./3)) #### Assuming the grid is cubic

   for t in range(len(A)):
       data.append(np.reshape(A[t,1:],(n,n,n)))

   OP=[[] for i in range(n)]

   for j in range(n):
       for k in range(n):
         s=choose_stack(data,uax,j,k)
         OP_stack=polar_op(s,90.0)
         for i in range(n):
             OP[i] += OP_stack[i]

   ### OP_stack will be a list of lists.
   ### The outer list is indexed by domain size (>0)
   ### and the inner list (for each size) has a collection of persistence times
   ### OP will be a similar list but over all stacks in a given plane.
   ### One can get, for instance, a distribution of persistence times for a given
   ### domain size. 
   for i in range(4):
       print("Maximum and average persistence times for domain with {0:2d} sites are {1:4d} and {2:10.2f}, respectively".format(i+1,max(OP[i]),sum(OP[i])/float(len(OP[i]))))

   from matplotlib import pyplot as plt

   j = 3
   fig, ax = plt.subplots(figsize =(10, 7))
   lmin = 0.0
   lmax = max(OP[j])
   dl = (lmax-lmin)/float(10)
   ax.hist(OP[j], bins = [lmin+i*dl for i in range(10)])

   # Show plot
   plt.show()
