import sys, os
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import astropy.io.fits as fits
import camb


def generate_power_spectra_CAMB(lmax, r=10**(-2), H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06, ns=0.965, lens_potential_accuracy=1, nt=0, ntrun=0, type_power='total', typeless_bool=False):
    """
    Return [Cl^TT, Cl^BB, Cl^EE, Cl^TE], in the format [ell_max, number_of_correlation]
    """
    pars = camb.CAMBparams(max_l_tensor=lmax, parameterization='tensor_param_indeptilt')
    pars.WantTensors = True

    pars.Accuracy.AccurateBB = True
    pars.Accuracy.AccuratePolarization = True
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, tau=tau)
    pars.InitPower.set_params(As=2e-9, ns=ns, r=r, parameterization='tensor_param_indeptilt', nt=nt, ntrun=ntrun) # May have modifications to add depending on the user's version of CAMB
    pars.max_eta_k_tensor = lmax + 100  # 15000  # 100

    pars.set_for_lmax(lmax, lens_potential_accuracy=lens_potential_accuracy)

    print("Calculating spectra from CAMB !")
    results = camb.get_results(pars)

    powers = results.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True, lmax=lmax)    
    if typeless_bool: # Return all spectra computed
        return powers
    return powers[type_power] # Return spectra corresponding to type_power



def save_fits_file_from_array(array_to_save, name_column_array, directory_path='', outname=''):
    """ Save arrays with corresponding name_column_array in directory_path+outname
        Note that the first index of the array_to_save will correspond to the number of columns recorded in the fits file

        To save c_ells, the array_to_save should be a vector in 1 dimension
        name_column_array must be a list of the name of the array, even if there is only 1 array
    """
    
    print('Shape arrau', array_to_save.shape)
    print('Recording 3-maps in', directory_path, outname)
    if directory_path=='':
        raise Exception('NO DIRECTORY PATH GIVEN TO SAVE FITS FILE !')
    
    dimension_array = array_to_save.shape[0]

    list_hdu = [fits.PrimaryHDU(array_to_save)]
    if len(array_to_save.shape) == 1:
        list_hdu.append(fits.BinTableHDU.from_columns([fits.Column(name=name_column_array[0], array=array_to_save, format='1D')]))
    else :
        list_hdu.append(fits.BinTableHDU.from_columns([fits.Column(name=name_column_array[i], array=array_to_save[i,:], format='1D') for i in range(dimension_array)]))

    hdu_list = fits.HDUList(list_hdu)
    # Header = fits.Header()


    if directory_path[-1] != '/':
        directory_path += '/'
    endfile = ''
    if outname[-5:] != '.fits':
        endfile = '.fits'
    hdu_list.writeto(directory_path+outname+endfile,overwrite=True)

# directory_path='/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions'
# outname = 'c_ell_file_lmax_4'

# save_fits_file_from_array(c_ell_final, ['cell'], directory_path=directory_path, outname=outname)
# c_ell_final = np.array([1.02198999e+03, 4.83707131e+02, 2.74501326e+02, 1.75106785e+02, # TT
#        3.35966904e-02, 2.31558395e-02, 1.30274054e-02, 6.30205108e-03, # EE
#        1.62109069e-04, 8.43358385e-05, 4.52237384e-05, 2.45132511e-05, # BB
#        2.73161593e+00, 1.57054793e+00, 9.05205777e-01, 5.26001677e-01]) # TE
