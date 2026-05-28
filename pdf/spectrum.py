import numpy as np
from scipy.integrate import simpson, trapezoid
import pyccl as ccl
from tqdm import tqdm
import Functions as f
import sys

'''
This is the main code for calculating the spectrum DM(M,z,theta) (eq. 19) and vectorizing it,
as well as vectorizing the angular halo density and halo bias function.
'''

class data_generator():
    def __init__(self,cospar):
        self.cospar = cospar
        self.Om,Ob,s8,ns = cospar[:4]
        self.h0 = cospar[4]
        self.z_frb = cospar[5]
        Oc = self.Om-Ob
        hmd_200c = ccl.halos.MassDef200c
        self.cosmo = ccl.Cosmology(Omega_c=Oc, Omega_b=Ob, h=self.h0, sigma8=s8, n_s=ns)
        self.mf = ccl.halos.MassFuncTinker10(mass_def=hmd_200c)  #unit Mpc^-3
        self.bm = ccl.halos.HaloBiasTinker10(mass_def=hmd_200c)
        self.z_resolution = int(30)
        self.M_resolution = int(96)
        self.l_resolution = int(1000)
        self.int_linPow = 57.5 #placeholder for the second integral in (27)

    def get_rvir(self,M_halo,z):
        '''
        Virial radious
        Inputs:
            M_halo : array-like or float
                Halo mass  (Msun/h).
            z : array-like or float
                Redshift.
        Output:
            rvir : array-like or float
                Virial radius (Mpc/h).
                Comoving
        '''
        RHOC = 2.776e11
        DELTAVIR = 200.0
        Om = self.Om
        M_halo = np.asarray(M_halo)
        z = np.asarray(z)
        rhoc_of_z = RHOC*(Om*(1.0+z)**3.0 + (1.0-Om))/(1.0+z)**3.0
        rvir = (3.0*M_halo/(4.0*np.pi*DELTAVIR*rhoc_of_z))**(1.0/3.0)
        return rvir

    def halo_bias(self,M_halo,z):
        '''
        Halo bias function
        Inputs:
            halo mass Msun/h
            redshift
        '''
        return self.bm(self.cosmo,M_halo,1/(1+np.array(z)))

    def ang_halo_density(self,M_halo,z):
        '''
        angular halo density
        input:
        halo mass Msun/h
        redshift
        '''
        co_dist = ccl.comoving_angular_distance(self.cosmo, 1/(1+z))*self.h0 # Mpc/h
        H_z = 100*ccl.h_over_h0(self.cosmo,1/(1+z)) # h km/s/Mpc
        co_volume = co_dist**2/H_z *ccl.physical_constants.CLIGHT/1000 #Mpc^3/h^3
        halo_mass_function = self.mf(self.cosmo,M_halo/self.h0,1/(1+np.array(z))) /self.h0**3 # h^3/Mpc^3
        return co_volume*halo_mass_function 

    def vec_spec(self,par,start,end):
        '''
        Y(M,z,l) vectorization, parallel run for different mass segments (start,end)
        '''
        spectrum = np.zeros((self.z_resolution,(end-start),1000))
        redshift = np.linspace(0.0001,self.z_frb,self.z_resolution)
        mass = np.logspace(8,16,self.M_resolution)
        for i_z,z in enumerate(tqdm(redshift)):
            j = 0
            for i_M in range(start,end):
                sp = f.modelG(z,mass[i_M],self.cospar,par)
                rmax =  5*self.get_rvir(mass[i_M],z)
                dom = np.linspace(0,rmax,1000)
                y_max = np.zeros(1000)
                y_max[:998] = np.sqrt(rmax**2 - dom[:998]**2)
                y_max[999] = 0
                x3_val = np.full(1000,0)
                spectrum[i_z,j,:] = sp.frb_host.dispersion_measure_in_halo_3d(dom,-y_max,x3_val)*(1+z)
                j+=1
        np.save(f"10/spectrum_{start}.npy",spectrum)

    def Yvec(self):
        '''
        angular halo density and halo bias vectorization 
        '''
        M_min,M_max = 8,16
        z_max = self.z_frb
        M_values= np.logspace(M_min,M_max,self.M_resolution) 
        z_values= np.linspace(0.0001,z_max,self.z_resolution)
        ahd = np.zeros((self.M_resolution,self.z_resolution))
        hb = np.zeros((self.M_resolution,self.z_resolution))
        
        #vectorizing the angular halo density
        for i,z in enumerate(z_values):
            ahd[:,i] = self.ang_halo_density(M_values,z)
        np.save("ahd.npy", ahd)
        #vectorizing the halo bias
        for i,z in enumerate(z_values):
            hb[:,i] = self.halo_bias(M_values,z)
        np.save("hb.npy", hb)


#parallel run for Y


if __name__ == "__main__":
    # Get start and end indices from command-line arguments
    cosmoparams = sys.argv[1:7]
    cospar = [float(i) for i in cosmoparams]
    params = sys.argv[7:16]
    par = [float(i) for i in params]
    
    start_idx_s = int(sys.argv[16])
    end_idx_s = int(sys.argv[17])

   # Perform computation for the segment
    data_generator(cospar).vec_spec(par,start_idx_s,end_idx_s)
