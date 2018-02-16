#!/usr/bin/python# My imports
from __future__ import division
import numpy as np
import argparse
import urllib
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.units import Quantity
from astroquery.vizier import Vizier


def output(header=False, parameters=None):
    """Create the output file 'synthresults.dat'

    Input
    -----
    overwrite - Overwrite the file
    header    - Only use True if this is for the file to be created
    """

    if header:
        hdr = ['star_name, star_teff, star_erteff, star_logg, star_erlogg, star_metal, dMcal']
        with open('star_gaia_hip_parallax.dat', 'w') as output:
            output.write('\t'.join(hdr)+'\n')
    else:
        with open('star_gaia_hip_parallax.dat', 'a') as output:
            output.write('\t'.join(map(str, parameters))+'\n')

def get_mass_radius(star, vmag, er_vmag, parallax, er_parallax, temp, er_temp, metal, er_metal):
    """Returns mass, radius, and age from padova interface.
       Enter star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
       Parallax is in mas.
       Be careful with the constrains in age!
    """
    url = 'http://stev.oapd.inaf.it/cgi-bin/param_1.3'
    #These are the parameters in the webpage to tune
    form_data = {'param_version': '1.3', 'star_name': star, 'star_teff': int(temp), 'star_sigteff': int(er_temp), 'star_feh': round(metal,2), 'star_sigfeh': round(er_metal,2),  'star_vmag': round(vmag,3), 'star_sigvmag': round(er_vmag,3), 'star_parallax': round(parallax,3), 'star_sigparallax': round(er_parallax,3),  'isoc_kind': 'parsec_CAF09_v1.1', 'kind_interp': '1', 'kind_tpagb': '0', 'kind_pulsecycle': '0', 'kind_postagb': '-1', 'imf_file': 'tab_imf/imf_chabrier_lognormal.dat', 'sfr_file': 'tab_sfr/sfr_const_z008.dat', 'sfr_minage': '0.1e9', 'sfr_maxage': '12.0e9', 'flag_sismo': '0', 'photsys_file': 'tab_mag_odfnew/tab_mag_ubvrijhk.dat', 'submit_form': 'Submit' }
    print('Parameters for Padova interface.')
    print(form_data)
    urllib.urlretrieve(url, 'parameters.html', lambda x,y,z:0, urllib.urlencode(form_data))

    #write results
    with open('parameters.html') as f:
        line = f.readlines()[15]

    line = line.replace(' .<p>Results for ', '')
    line = line.replace('Mass=', '')
    line = line.replace('&plusmn;', ' , ')
    line = line.replace(' ','')
    line = line.replace('<i>R</i>=','')
    line = line.replace(':Age=',',')
    line = line.replace('<i>M</i>&#9737','')
    line = line.replace('<i>R</i>&#9737','')
    line = line.replace('log<i>g</i>=','')
    line = line.replace('(cgs)','')
    line = line.replace('Gyr','')
    line = line.split(',')

    name =     str(line[0])
    age =      float(line[1])
    erage =    float(line[2])
    mass =     float(line[3])
    ermass =   float(line[4])
    logg_p =   float(line[5])
    erlogg_p = float(line[6])
    radius =   float(line[7])
    erradius = float(line[8])
    print name, ' done.'
    return mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p

def gaia_hip_paralax_query(star):
    '''Query objects by name and get the gaia parallax'''

    #result_table = Simbad.query_object(object)
    customSimbad = Simbad()
    customSimbad.add_votable_fields('flux(V)')
    customSimbad.get_votable_fields()
    result_table = customSimbad.query_object(star)
    ra   = result_table[0][1]
    dec  = result_table[0][2]
    vmag = result_table[0][-1]
    vmag = round(vmag, 3)
    c = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.degree))
    radius = Quantity(0.001388888888888889, u.deg)
    j = Gaia.cone_search_async(c, radius)
    r = j.get_results()
    p = r['parallax']
    p_er = r['parallax_error']
    try:
        plx = round(p[0], 3)
        plx_er = round(p_er[0], 3)
    except IndexError as err:
        plx, plx_er = ('nan', 'nan')
    v = Vizier(catalog='I/311/hip2')
    result = v.query_region(c, radius=radius)
    try:
        hip    = result['I/311/hip2'][0][6]
        hip_er = result['I/311/hip2'][0][7]
    except TypeError as err:
        hip, hip_er = ('nan', 'nan')
    return vmag, plx, plx_er, hip, hip_er

def mass_radius_age(star, teff, erteff, logg, erlogg, feh, erfeh):
    '''Query objects by name and get the gaia parallax'''

    vmag, p, p_er, hip, hip_er = gaia_hip_paralax_query(x)
    if p != 'nan':
        parallax = p
        er_parallax = p_er
    if p == 'nan' and hip != 'nan':
        parallax = hip
        er_parallax = hip_er
    if p == 'nan' and hip == 'nan':
        parallax = 'nan'
        er_parallax = 'nan'
    mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p = get_mass_radius(star, vmag, er_vmag, parallax, er_parallax, teff, er_teff, feh, er_feh)
    return mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p

def bcflow(teff):
    """Table 1 from Torres 2010: Bolometric Corrections by Flower (1996) as a Function of Temperature:
    BC_V = a + b(log T_eff) + c(log T_eff)^2 + sdot sdot sdot
    Coefficient	log T_eff < 3.70	3.70 < log T_eff < 3.90	   log T_eff>3.90
    a	       -0.190537291496456E+05	-0.370510203809015E+05	-0.118115450538963E+06
    b	        0.155144866764412E+05	 0.385672629965804E+05	 0.137145973583929E+06
    c	       -0.421278819301717E+04	-0.150651486316025E+05	-0.636233812100225E+05
    d	        0.381476328422343E+03	 0.261724637119416E+04	 0.147412923562646E+05
    e	                    ... 	    -0.170623810323864E+03	-0.170587278406872E+04
    f	                    ... 	            ... 	         0.788731721804990E+02
    """

    a = [-0.190537291496456e+05, -0.370510203809015e+05, -0.118115450538963e+06]
    b = [0.155144866764412e+05, 0.385672629965804e+05, 0.137145973583929e+06]
    c = [-0.421278819301717e+04, -0.150651486316025e+05, -0.636233812100225e+05]
    d = [0.381476328422343e+03, 0.261724637119416e+04, 0.147412923562646e+05]
    e = [-0.170623810323864e+03, -0.170587278406872e+04]
    f = [0.788731721804990e+02]

    lteff= np.log10(teff)
    if lteff < 3.7:
        bc = a[0] + (b[0]*lteff) + (c[0]*(lteff**2)) + (d[0]*(lteff**3))
    elif (lteff >= 3.7) and (lteff < 3.9):
        bc = a[1] + (b[1]*lteff) + (c[1]*(lteff**2)) + (d[1]*(lteff**3)) + (e[0]*(lteff**4))
    elif lteff >= 3.9:
        bc = a[2] + (b[2]*lteff) + (c[2]*(lteff**2)) + (d[2]*(lteff**3)) + (e[1]*(lteff**4)) + (f[0]*(lteff)**5)
    return bc

def logg_trigomonetric(teff, mass, v, bc, par, dpar, dteff, dmass):
    """Calculate the trigonometric logg and error"""
    #np.geterr()
    if mass == 'nan':
        logg, dlogg = 'nan', 'nan'

    else:
        e = 2.718281828
        logg  = 4.44 + np.log10(mass) + (4.0*np.log10(teff/5777.)) + (0.4*(v + bc)) + (2.0*np.log10(par/1000.0)) + 0.108
        logg  = np.round(logg,2)
        dlogg = np.sqrt(((dmass*np.log10(e))/mass)**2 + ((4.*dteff*np.log10(e))/teff)**2 + ((2.*0.05*np.log10(e))/par)**2)
        dlogg = np.round(dlogg,2)
    return logg, dlogg


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='set flag if there is information on the parallax and V mag.')
    parser.add_argument('star_file', type=str, help ='Parallaxes and Vmags are included by default')
    parser.add_argument('-o', '--out', help ='Parallaxes and Vmags are included by default', choices=['plx', 'grav'])
    args = parser.parse_args()

    star = np.genfromtxt(args.star_file, dtype=None, delimiter='\t', comments='#', usecols=range(7), names=['star', 'teff', 'erteff', 'logg', 'erlogg', 'feh', 'erfeh'])
    star_name    = star['star']
    star_teff    = star['teff']
    star_erteff  = star['erteff']
    star_logg    = star['logg']
    star_erlogg  = star['erlogg']
    star_metal   = star['feh']
    star_ermetal = star['erfeh']

for i, x in enumerate(star_name):
    if args.out == 'plx':
        vmag, p, p_er, hip, hip_er = gaia_hip_paralax_query(x)
        print(x, vmag, p, p_er, hip, hip_er)
        parameters = [x, vmag, p, p_er, hip, hip_er]
        output(parameters=parameters)
    if args.out == 'grav':
        vmag, p, p_er, hip, hip_er = gaia_hip_paralax_query(x)
        er_vmag = 0.03
        mass, ermass, radius, erradius, age, erage, logg_p, erlogg_p = get_mass_radius(star, vmag, er_vmag, parallax, er_parallax, teff, er_teff, feh, er_feh)
        logg, dlogg = logg_trigomonetric(teff, mass, v, bc, par, dpar, dteff, dmass)
        parameters = [x, vmag, p, p_er, hip, hip_er]
        output(parameters=parameters)
