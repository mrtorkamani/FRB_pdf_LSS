import numpy as np
from scipy.integrate import simpson, trapezoid
import pyccl as ccl
from tqdm import tqdm
from scipy.interpolate import interp1d
import Functions as f
import sys

'''
This is the main code for calculating the spectrum Y(M,z,l) (eq. 29) and vectorizing it,
'''

class data_generator():
    def __init__(self,cospar):
        self.cospar = cospar
        self.Om,Ob,s8,ns = cospar[:4]
        self.h0 = cospar[4]
        self.z_frb = cospar[5]
        l_max = interp1d([.1,0.7,1.5,3,4,5],[1,1,.1,.1,.1,.1],kind='linear')
        self.l_max = l_max(self.z_frb)
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
        #return rvir
        
        hmd_200c = ccl.halos.MassDef200c
        return hmd_200c.get_radius(self.cosmo,M_halo/self.h0,1.0/(z + 1)) *self.h0*(1+z)
    
    def Y(self,start,end):
        '''
        Y(M,z,l) vectorization, parallel run for different mass segments (start,end)'''
        data = np.load("10/spectrum.npy")
        Y = np.zeros(((end-start),self.M_resolution,self.l_resolution),dtype=complex)
        redshift = np.linspace(0.0001,self.z_frb,self.z_resolution)
        mass = np.logspace(8,16,self.M_resolution)
        l_val= np.linspace(0,1,self.l_resolution)
        l_trans = l_val**5
        l_values = l_trans*self.l_max
        l_values = l_values[None, None, :]  # Shape: (1, 1, L)
        i =0
        for i_z in tqdm(range(start,end)):
            pred = data[i_z,:,:,None]
            rmax = 5*self.get_rvir(mass,redshift[i_z])
            dom = np.linspace(0,rmax,1000)
            d_A = ccl.comoving_angular_distance(self.cosmo,1/(1+np.array(redshift[i_z])))*self.h0
            theta = np.arctan(dom.T/d_A)
            theta = theta[..., None]
            integrand = 2*np.pi* theta *(np.exp(1j*l_values* pred)-1)
            integral = simpson(integrand,x=theta, axis=1)
            Y[i,:,:] = integral
            i+=1
        np.save(f"10/Y_{start}.npy",Y)


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
    data_generator(cospar).Y(start_idx_s,end_idx_s)
