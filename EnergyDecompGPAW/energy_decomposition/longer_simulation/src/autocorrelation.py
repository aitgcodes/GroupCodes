from scipy import signal

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
