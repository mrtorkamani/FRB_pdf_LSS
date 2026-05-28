import numpy as np
from scipy.integrate import simpson, trapezoid
import pyccl as ccl
from tqdm import tqdm
from astropy import units as u
from astropy import constants
from scipy.interpolate import interp1d


path = "10/"

class One_Point_PDF():
    def __init__(self,cospar):
        self.cospar = cospar
        Om,Ob,s8,ns = cospar[:4]
        self.h0 = cospar[4]
        self.z_frb = cospar[5]
        Oc = Om-Ob
        hmd_200c = ccl.halos.MassDef200c
        self.cosmo = ccl.Cosmology(Omega_c=Oc, Omega_b=Ob, h=self.h0, sigma8=s8, n_s=ns)
        self.mf = ccl.halos.MassFuncTinker10(mass_def=hmd_200c)  #unit Mpc^-3
        self.bm = ccl.halos.HaloBiasTinker10(mass_def=hmd_200c)
        self.z_resolution = int(30)
        self.M_resolution = int(96)
        self.l_resolution = int(1000)
        self.int_linPow = 57.5 #placeholder for the second integral in (27)

    def p(self):
        '''
        Sum over all halo masses and redshifts 
        (eq 31)
        '''
        M_min,M_max = 8,16
        z_max = self.z_frb
        M_values= np.logspace(M_min,M_max,self.M_resolution) 
        z_values= np.linspace(0.0001,z_max,self.z_resolution)
        
        Y = np.load(path+"Y.npy")
        ahd =np.load(path+"ahd.npy")
        integrand =  ahd[:,:,None]*Y
        integral_M = simpson(integrand, x=np.log10(M_values), axis=1)
        
        integral = simpson(integral_M, x=z_values, axis=0)
        
        return np.exp(integral)
    def alpha(self):
        '''
        alpha function
        (eq 35)
        Output:
            alpha : array-like
                alpha function (unitless)
        '''
        M_min,M_max = 8,16
        M_values= np.logspace(M_min,M_max,self.M_resolution) 
        
        Y = np.load(path+"Y.npy")
        ahd =np.load(path+"ahd.npy")
        hb =np.load(path+"hb.npy")
        #integral = np.zeros((self.z_resolution,self.l_resolution),dtype=complex)
        
        integrand =  hb[:,:,None]*ahd[:,:,None]*Y
        integral_M = simpson(integrand, x=np.log10(M_values), axis=1)
        return integral_M
    def growth_factor_integral(self):
        '''
        First integral (eq 47)
        '''
        z_max = self.z_frb
        z_values = np.linspace(0.0001, z_max, self.z_resolution)
        growth_factors = ccl.growth_factor(self.cosmo, 1 / (1 + z_values))
        h_values = 100 * ccl.h_over_h0(self.cosmo, 1 / (1 + z_values))
        integrand = (h_values[:, None] * self.alpha()**2 * growth_factors[:, None]**2)
        integral = simpson(integrand, x=z_values, axis=0)
        return integral
        
    def p_cl(self):
        '''
        clustered pdf from unclustered pdf
        (eq47)
        '''
        return self.p() *np.exp(1/2 * self.growth_factor_integral()*self.int_linPow/(ccl.physical_constants.CLIGHT/1000))

    ####################MEAN##########
    '''
    This part is for calculating the mean DM, which is used for the x-axis of the final PDF.
    '''

    def weight(self,z):
        Ob = self.cospar[1]
        H0 = self.h0 * 100 * u.km / u.s / u.Mpc
        m2 = 1/u.m/u.m
        prefac = 3.0*constants.c.value / \
        (constants.m_p.value*constants.G.value*8.0*np.pi) * \
        (H0.to(1/u.s)).value * \
        m2.to(u.parsec/u.cm**3).value*Ob
        f_IGM = 0.9  # keep constant for now at redshifts < 1. Can be calculated from https://github.com/FRBs/FRB/blob/main/frb/dm/igm.py as f_diffuse
        f_He = 0.24  # Helium fraction
        f_H = 1. - f_He
        f_e = (f_H + 1/2. * f_He)  # electron fraction
        result = f_e * f_IGM * prefac*(1+z)/ccl.h_over_h0(self.cosmo,1/(1+z))
        return result
    def mean(self,z):
        zval = np.linspace(0,z,1000)
        integrand = self.weight(zval)
        integral = simpson(integrand,x=zval)
        return integral
    #####################MEAN############

def final_result(cospar):
    '''
    This is the final function for calculating the PDF of DM. inverse fourier transform of p_cl (eq. 47) to get the PDF in real space.
    Output:
    ----------
    pdf: 2D array-like
        The final PDF as a function of DM.
    pcl: 2D array-like
        P_cl.
    '''
    data = One_Point_PDF(cospar).p_cl()
    d = np.linspace(0,1,1000)
    y = d**5
    l_max = interp1d([.1,.7,1.5,3,4,5],[1,1,.1,.1,.1,.1],kind='linear')
    redshift = cospar[5]
    d = y*l_max(redshift)
    firstp = np.zeros((2,1000),dtype=complex)
    firstp[0] =d    
    firstp[1] =data

    spline = interp1d(d,data,kind='linear')
    R =  np.linspace(0,l_max(redshift),100000)
    p1 = spline(R)

    plambda = np.zeros((2,100000),dtype=complex)
    plambda[0]=R
    plambda[1]=p1

    final = []
    mean = One_Point_PDF(cospar).mean(redshift)
    dom = np.linspace(0,3*mean,1000)

    for i in dom:
        integrand = p1*np.exp(-1j*R*i)
        integral = simpson(integrand,x=R)
        final.append(integral.real)
    norm_final = np.trapz(final,dom)
    final /= norm_final
    final = np.array(final)
    pdf = np.zeros((2,1000))
    pdf[0] = dom
    pdf[1] = final
    return pdf, firstp
