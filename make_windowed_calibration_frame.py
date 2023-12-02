"""
Take a full frame calibration frame and a set of subframe
coordinates in FITS standard and convert produce a windowed
calibration frame for process
"""
import argparse as ap
import jastro.housekeeping as jhk

# pylint: disable=invalid-name

def arg_parse():
    """
    Parse command line arguments
    """
    p = ap.ArgumentParser()
    p.add_argument('calib_frame',
                   help='calibration filename',
                   type=str)
    p.add_argument('frame_coords',
                   help='frame coords in FITS format, e.g. [x1:x2,y1:y2]',
                   type=str)
    return p.parse_args()

if __name__ == "__main__":
    args = arg_parse()
    # get frame coords in numpy format
    sf = args.frame_coords.split(',')
    sfx = sf[0].split(":")
    sfy = sf[1].split(':')
    x1 = int(sfx[0][1:]) - 1
    x2 = int(sfx[1])
    y1 = int(sfy[0]) - 1
    y2 = int(sfy[1][:-1])

    # load the calibration frame
    data, hdr = jhk.load_fits_image(args.calib_frame)

    # apply the window
    data_sf = data[y1:y2, x1:x2]

    # remake the fits header
    new_hdr = {}
    new_hdr['CAM-WIND'] = args.frame_coords
    new_hdr['EXPTIME'] = hdr['EXPTIME']

    # write out the windowed calibration frame
    new_filename = f"{args.calib_frame.split('.f')[0]}-X{x1+1}_{x2}-Y{y1+1}_{y2}.fits"
    jhk.write_fits_image(new_filename, data_sf, new_hdr, clobber=True)
