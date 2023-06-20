import numpy as np
import healpy as hp
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import pymaster as nmt

def save_fits_file_from_array(array_to_save, name_column_array, directory_path='', outname=''):
    """ Save arrays with corresponding name_column_array in directory_path+outname
        Note that the first index of the array_to_save will correspond to the number of columns recorded in the fits file

        To save c_ells, the array_to_save should be a vector in 1 dimension
        name_column_array must be a list of the name of the array, even if there is only 1 array
    """
    
    print('Shape array', array_to_save.shape)
    print('Recording file in', directory_path, outname)
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

def save_full_covariance_matrix(directory_path, outname, covariance_dict, list_correl=['T','E','B']):

    first_letter = list_correl[0]
    first_correl = first_letter+first_letter+'_'+first_letter+first_letter
    lmax = covariance_dict[first_correl].shape[0]

    full_covariance = np.zeros((len(list_correl)*lmax, len(list_correl)*lmax))

    for i in range(len(list_correl)):
        for j in range(len(list_correl)):
            correl = list_correl[i] + list_correl[i] + '_' + list_correl[j] + list_correl[j]
            full_covariance[i*lmax:(i+1)*lmax, j*lmax:(j+1)*lmax] = covariance_dict[correl]
    
    save_fits_file_from_array(full_covariance.ravel(order='C'), ['covariance_matrix'], directory_path=directory_path, outname=outname)


