import guesthost.lattice as lat

# Create lattice                                                       
lattice = lat.HPLattice(motifs, cell, nx, ny, nz)                      

ind = (0, 0, 0)

# Compute local cell parameters                                        
params = lattice.ucell_localcellparameters(ind)               
                                                                       
# Compute MA orientation                                               
theta, phi = lattice.ucell_theta_phi(ind, dir=0)
                                                                       
# Compute order parameters                                             
eta_lat = lattice.ucell_η_lat(ind, dir_coup=0)       
eta_ma = lattice.ucell_η_MA(ind, dir=0)         
