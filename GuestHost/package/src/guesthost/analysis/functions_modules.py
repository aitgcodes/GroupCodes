import math
import numpy as np
from scipy import fftpack
from scipy import signal

def cord_extract_all_t(ind_t0,arr):
    cord=[]
    for i in range(len(arr)):
        cord.append(arr[i][ind_t0])
    cord=np.array(cord)
    return cord

def print_file(arr,dt,file_name):
    time=np.array([x*dt for x in range(len(arr[0]))])
    file=open(file_name,"w")
    for i in range(len(arr[0])):
        file.write("  {: .16E}".format(time[i]))
        for j in range(len(arr)):
            file.write("  {: .16E}".format(arr[j][i]))
        file.write("\n")
    file.close()

def put_in_one_array(arr,file_name):
    file=open(file_name,"w")
    arr=arr.reshape(1,arr.shape[0]*arr.shape[1])[0]
    for i in range(len(arr)):
            file.write("  {: .16E}".format(arr[i]))
            file.write(" \n")
    file.close()

def self_correlation_P1(B):
    N=len(B)
    C=[]
    for l in range(0,N,1):
        sumP1=0.0
        sumP1_D=0.0
        for j in range(0,N-l,1):
            sumP1=sumP1+np.dot(B[j],B[j+l])
            sumP1_D=sumP1_D+np.dot(B[j],B[j])
            #corrP1=sumP1/(N-l)
            corrP1=sumP1/sumP1_D
        C.append(corrP1)
    return np.array(C)

def L2_legendre(x):
    y=0.5*(3.0*(x**2.0)-1.0)
    return y

def self_correlation_P2(B):
    N=len(B)
    C=[]
    for l in range(0,N,1):
        sumP2=0.0
        sumP2_D=0.0
        for j in range(0,N-l,1):
            sumP2=sumP2+L2_legendre(np.dot(B[j],B[j+l]))
            sumP2_D=sumP2_D+L2_legendre(np.dot(B[j],B[j]))
#            corrP2=0.5*(3.0*(sumP2/(N-l))-1.0)
            corrP2=sumP2/sumP2_D
        C.append(corrP2)
    return np.array(C)

def autocorrelate(a,normalize=True):
    ### Expects 1-d array a of size N.
    ### If normalize set to True then divides by zero time component.
    ### Returns 1-d array only positive and zero lag components
    ### of size N (t=0 to t=N-1).
    nt = a.shape[0]
    lags = signal.correlation_lags(nt,nt)
    nlag = lags.size//2
    lags = lags[nlag:]
    lagden = np.array([(nt-lag) for lag in lags])
    ac = signal.correlate(a,a)
    ac = ac[nlag:]/lagden
    norm = max(ac[0],1e-8)
    if  normalize:
        ac = ac/norm                                                           
        #print(ac.size)
    return ac

def print_file_average(arr,dt,file_name):
    time=np.array([x*dt for x in range(len(arr[0]))])
    avg=np.mean(arr,axis=0)
    file=open(file_name,"w")
    for i in range(len(avg)):
        file.write("  {: .16E}".format(time[i]))
        file.write("  {: .16E}".format(avg[i]))
        file.write("\n")
    file.close()

def projection_of_point_on_plane(l1,l2,l3,P4):
    l1=np.array(l1)
    l2=np.array(l2)
    l3=np.array(l3)
    P4=np.array(P4)
    P12=l2-l1
    P13=l3-l1
    v=P4-l1
    w=np.cross(P12,P13)
    v_parallel=(np.dot(v,w)/(np.linalg.norm(w))**2.0)*w
    v_perpendicular=v-v_parallel
    p=l1+v_perpendicular
#    p=np.array(p)
    return p

def distance_of_point_on_plane(l1,l2,l3,P4):
    proj_P=projection_of_point_on_plane(l1,l2,l3,P4) 
    v=np.array(P4-proj_P)
    dis=np.linalg.norm(v)   
    s=P4-proj_P
    d=[dis,s]
    # return np.array(d)
    return d

def bandpass(fL,fH,b,s):
    N = int(np.ceil((4 / b)))
    if not N % 2: N += 1  # Make sure that N is odd.
    n = np.arange(N)
    # low-pass filter
    hlpf = np.sinc(2 * fH * (n - (N - 1) / 2.))
    hlpf *= np.blackman(N)
    hlpf = hlpf / np.sum(hlpf)
    # high-pass filter
    hhpf = np.sinc(2 * fL * (n - (N - 1) / 2.))
    hhpf *= np.blackman(N)
    hhpf = hhpf / np.sum(hhpf)
    hhpf = -hhpf
    hhpf[int((N - 1) / 2)] += 1
    # datageneration
    h = np.convolve(hlpf, hhpf)
    new_signal = np.array(np.convolve(s, h))
    return np.array(new_signal)

def FFT_transform(data,time_step_constant):
    data_fft=fftpack.fft(data)
    # And the power (sig_fft is of complex dtype)
    power = np.abs(data_fft)**2
    # The corresponding frequencies
    sample_freq = fftpack.fftfreq(data.size, time_step_constant)
    pos_mask = np.where(sample_freq > 0)
    freqs = sample_freq[pos_mask]
    peak_freq = freqs[power[pos_mask].argmax()]
    # rearranging 
    power1=power[:len(power)//2]
    sample_freq1=sample_freq[:len(sample_freq)//2]
    power2=power[len(power)//2:]
    sample_freq2=sample_freq[len(sample_freq)//2:]
    power=power2.tolist()+power1.tolist()
    sample_freq=sample_freq2.tolist()+sample_freq1.tolist()
    return [np.array([sample_freq,power]),peak_freq]

def MA_theta_phi(a,b,c,v):
    n=np.cross(a,b)
    F=np.array([[-n[1]/n[0],-n[2]/n[0]],[1.0,0.0],[0.0,1.0]])
    t=F @ np.linalg.inv(F.T@F) @ F.T
    P=(t@v.T).T
#    m1=a[0]*F[0][0]+a[1]*F[1][0]+a[2]*F[2][0]    #r=c1(basis1)+c2(basis2) >  r.a=0 >  c1m1+c2m2=0   >  c1/c2=(-m2/m1)=k
#    m2=a[0]*F[0][1]+a[1]*F[1][1]+a[2]*F[2][1]
    m1=np.dot(a,F.T[0])
    m2=np.dot(a,F.T[1])
    k=(-m2/m1)
#    R=[k*F[0][0]+F[0][1],k*F[1][0]+F[1][1],k*F[2][0]+F[2][1]]
    R=k*F.T[0]+F.T[1]
    r1=R/np.linalg.norm(R)
    v1=v/np.linalg.norm(v)
    c1=c/np.linalg.norm(c)
    a1=a/np.linalg.norm(a)
    p1=P/np.linalg.norm(P)
    theta=(180.0/np.pi)*np.arccos(np.dot(v1,c1))
#    phi=(180/pi)*atan2((dotpdt(p1,r1)),(dotpdt(p1,a1)))
    phi=(180/np.pi)*np.arctan2((np.dot(p1,r1)),(np.dot(p1,a1)))
    return [theta, phi]

def scatter_theta_phi(T,P,d):
    xaxis=np.array([float((d/2)+i) for i in range(-int(d),180+int(d),int(d))])
    yaxis=np.array([float((d/2)+i) for i in range(-180-int(d),180+int(d),int(d))])
    final=[]
    X=[]
    Y=[]
    for s in range(len(xaxis)):
        for t in range(len(yaxis)):
            count=0
            for i,j in zip(range(len(T)),range(len(P))):
                if (T[i]>(xaxis[s]-(d/2)) and T[i]<(xaxis[s+1]-(d/2))) and (P[j]>(yaxis[t]-(d/2)) and P[j]<(yaxis[t+1]-(d/2))):
                   count=count+1
            #final.append([xaxis[s],yaxis[t],count])
            X.append(xaxis[s])
            Y.append(yaxis[t])
            final.append(count)
    return np.array([X,Y,final])
            #print(xaxis[s], yaxis[t], count)

def dihedral2(p):
   b = p[:-1] - p[1:]
   b[0] *= -1
   v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
   # Normalize vectors
   v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
   b1 = b[1] / np.linalg.norm(b[1])
   x = np.dot(v[0], v[1])
   m = np.cross(v[0], b1)
   y = np.dot(m, v[1])
   return np.degrees(np.arctan2( y, x ))

def dihedral3(p):
   b = p[:-1] - p[1:]
   b[0] *= -1
   v = np.array( [np.cross(v,b[1]) for v in [b[0], b[2]] ] )
   # Normalize vectors
   v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
   return np.degrees(np.arccos( v[0].dot(v[1]) ))

def periodic_boundary_phi(P,d):
#    d=10.0
    yaxis=np.array([float((d/2)+i) for i in range(-180-int(d),180+int(d),int(d))])
    final=[]
    Y=[]
    for t in range(len(yaxis)):
        count=0
        for j in range(len(P)):
            if P[j]>(yaxis[t]-(d/2)) and P[j]<(yaxis[t+1]-(d/2)):
               count=count+1
        Y.append(yaxis[t])
        final.append(count)
    Y=Y[1:-1]
    final=final[1:-1]
    Pmin_ind=np.argmin(final)
    Pc=Y[Pmin_ind]+d/2
#    print(Pc)
#    NewP=[]
    for j in range(len(P)):
        if P[j]<Pc:
           P[j]=P[j]+360.0
    return (P,Pc)

def rpvecs(h):
    g = np.zeros((3,3),dtype=float)
    v = np.cross(h[1,:],h[2,:])
    vol = np.dot(h[0,:],v)
    #print(vol)
    lst = [[0,1,2],[1,2,0],[2,0,1]]
    for [i,j,k] in lst:
       # print(i,j,k)
        g[i,:] = np.cross(h[j,:],h[k,:])/vol
       # print(g[i,:])
    return g

def magnitude(v):
    vmag = np.sqrt(np.dot(v,v))
    return vmag

def proj2plane(iax,r,h):
    ### Calculate the reciprocal vector set to h
    ### Columns of h are the lattice vectors
    g = rpvecs(h)
    i = iax
    j = (iax+1)%3
    k = (iax-1)%3
    rad2deg = 180.0/np.pi
    rmag = magnitude(r)
    amag = magnitude(h[i,:])
    bmag = magnitude(h[j,:])
    cmag = magnitude(h[k,:]) 
    #### Calculate theta
    anga = math.acos(np.dot(r,h[i,:])/(rmag*amag))
    vprp = np.dot(r,g[i,:])*h[i,:]
    vpl = r - vprp
    vplmag = max(magnitude(vpl), 1e-8)
    phase = 1.0
    angb = np.dot(vpl,h[j,:])/(vplmag*bmag)
    angc = np.dot(vpl,g[k,:]) ## This is just the component of vector c in vector r
    if angc < 0: ## To ensure the phi measured counter-clockwise is kept positive
       phase = -1
    angb = phase*math.acos(angb)
    return (vpl,anga*rad2deg,angb*rad2deg)

