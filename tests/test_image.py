import numpy as np
import exoscene.image as image
import matplotlib.pyplot as plt
import pkg_resources
import os
import astropy.io.fits as fits
import astropy.units as u

## Load CGI PSF model
hlc_psf_path =  pkg_resources.resource_filename('exoscene', 'data/cgi_hlc_psf')
psf_cube_fname = os.path.join(hlc_psf_path, 'hlc_os11_psfs_oversampled.fits')
psf_r_fname = os.path.join(hlc_psf_path, 'hlc_os11_psfs_radial_offsets.fits')
psf_angle_fname = os.path.join(hlc_psf_path, 'hlc_os11_psfs_azimuth_offsets.fits')

psf_cube = fits.getdata(psf_cube_fname)
psf_hdr = fits.getheader(psf_cube_fname)
hires_pixscale_as = psf_hdr['PIXAS'] * u.arcsec
hires_pixscale_LoD = psf_hdr['PIXLAMD']
        
r_offsets_LoD = fits.getdata(psf_r_fname)
r_offsets_as = r_offsets_LoD * hires_pixscale_as / hires_pixscale_LoD
angles = fits.getdata(psf_angle_fname)

# Save radial and angle offsets of PSF cube
r_offsets_LoD = fits.getdata(psf_r_fname)
r_offsets_as = r_offsets_LoD * hires_pixscale_as / hires_pixscale_LoD
angles = fits.getdata(psf_angle_fname)

# Save PSF cube and array center
offset_psfs = psf_cube
Np = offset_psfs.shape[-1]
cx = Np // 2 # Array center in zero-based indices

# print(f"angles[0]: {angles[0]}")
# print(f"r_offsets_as[10]: {r_offsets_as[10]}")
# plt.imshow(psf_cube[0,10],origin='lower')
# plt.colorbar()
# plt.show()
# plt.close()

deltax_as, deltay_as = 0.250, 0.250

planet_psf = image.get_hires_psf_at_xy_os11(
    offset_psfs, r_offsets_as.value, angles,
    hires_pixscale_as.value, deltax_as, deltay_as)

# plt.imshow(planet_psf,origin='lower')
# plt.colorbar()
# plt.show()
# plt.close()

def test_get_hires_psf_at_xy_os11():

    n_offsets = 4
    n_angles = 4
    pixscale_as = 0.01 
    imsize = [101,101]
    cen_px = (imsize[0] - 1) / 2

    offsets_as = np.linspace(0.1,0.4,n_offsets)
    angles = np.linspace(0,270,n_angles)

    # Normalized fake "PSFs"
    up = np.array([[0,0,0],
                   [1,1,1],
                   [0,1,0]],dtype=np.float64) / 4
    right = np.array([[0,1,0],
                      [0,1,1],
                      [0,1,0]],dtype=np.float64) / 4
    down = np.array([[0,1,0],
                     [1,1,1],
                     [0,0,0]],dtype=np.float64) / 4
    left = np.array([[0,1,0],
                     [1,1,0],
                     [0,1,0]],dtype=np.float64) / 4
            
    offax_psf_cube = np.zeros((n_angles,n_offsets,*imsize))

    for i_angle, angle in enumerate(angles):

        angle = angle % 360.

        if angle < 45:
            psf = right
        elif angle < 135:
            psf = up
        elif angle < 225:
            psf = left
        elif angle < 315:
            psf = down
        else:
            psf = right

        for i_offset, offset_as in enumerate(offsets_as):

            x_as = offset_as * np.cos(angle * 2 * np.pi / 360.)
            y_as = offset_as * np.sin(angle * 2 * np.pi / 360.)
            
            x_pix = x_as / pixscale_as
            y_pix = y_as / pixscale_as

            x_start= int(x_pix - 1 + cen_px)
            x_end = int(x_pix + 2 + cen_px)
            y_start= int(y_pix - 1 + cen_px)
            y_end = int(y_pix + 2 + cen_px)
            
            offax_psf_cube[i_angle,i_offset,y_start:y_end,x_start:x_end] = psf

        
    print(f"angles[0]: {angles[0]}")
    print(f"offset_as[2]: {offsets_as[2]}")
    plt.imshow(offax_psf_cube[0,2],origin='lower')
    plt.colorbar()
    plt.show()
    plt.close()
    print(f"angles[1]: {angles[1]}")
    print(f"offset_as[2]: {offsets_as[2]}")
    plt.imshow(offax_psf_cube[1,2],origin='lower')
    plt.colorbar()
    plt.show()
    plt.close()
    print(f"angles[2]: {angles[2]}")
    print(f"offset_as[2]: {offsets_as[2]}")
    plt.imshow(offax_psf_cube[2,2],origin='lower')
    plt.colorbar()
    plt.show()
    plt.close()
    print(f"angles[3]: {angles[3]}")
    print(f"offset_as[2]: {offsets_as[2]}")
    plt.imshow(offax_psf_cube[3,2],origin='lower')
    plt.colorbar()
    plt.show()
    plt.close()

    # Show that generation and shifts occur correctly
    


test_get_hires_psf_at_xy_os11()