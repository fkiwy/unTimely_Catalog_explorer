import os
from os.path import exists
import math
import sys
import tempfile
import traceback
import subprocess
import numpy as np
from PIL import Image, ImageOps, ImageDraw
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy.time import Time
from astropy.utils.data import download_file
# from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp


ALTERNATE_SCAN = 0
SEPARATE_SCAN = 1
MERGE_SCAN = 2


def search_by_coordinates(target_ra, target_dec, box_size=100, finder_charts=False, overlays=False, overlay_color='green', overlay_labels=False,
                          overlay_label_color='red', neowise_contrast=3, show_result_table_in_browser=True, save_result_table=False, result_table_format='ascii',
                          result_table_extension='dat', open_finder_charts=False, finder_charts_format='pdf', animated_gif=False, scan_dir_mode=ALTERNATE_SCAN,
                          directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300):

    class ImageBucket:
        def __init__(self, data, x, y, band, year_obs, wcs, overlay_label, overlay_ra=None, overlay_dec=None):
            self.data = data
            self.x = x
            self.y = y
            self.band = band
            self.year_obs = year_obs
            self.wcs = wcs
            self.overlay_label = overlay_label
            self.overlay_ra = overlay_ra
            self.overlay_dec = overlay_dec

    def process_image_data(hdu):
        try:
            data = hdu.data
            wcs, shape = find_optimal_celestial_wcs([hdu], frame='icrs')
            data, _ = reproject_interp(hdu, wcs, shape_out=shape)
            position = SkyCoord(ra*u.deg, dec*u.deg)
            cutout = Cutout2D(data, position, box_size*u.arcsec, wcs=wcs, mode='partial')
            data = cutout.data
            wcs = cutout.wcs
            x, y = wcs.world_to_pixel(position)
            return data, x, y, wcs
        except Exception:
            print('A problem occurred while creating an image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
            print(traceback.format_exc())
            return None, 0, 0, None

    def create_rgb_image(r, g, b):
        try:
            xmax, ymax = g.shape
            vmin, vmax = get_min_max(g, lo=neowise_contrast, hi=100-neowise_contrast)
            image = Image.fromarray(make_lupton_rgb(r, g, b, minimum=vmin, stretch=vmax-vmin, Q=0)).resize((zoom*xmax, zoom*ymax), Image.NONE)
            image = image.transpose(Image.FLIP_TOP_BOTTOM)
            image = ImageOps.invert(image)
            return image
        except Exception:
            print('A problem occurred while creating a color image for object ra={ra}, dec={dec}'.format(ra=ra, dec=dec))
            print(traceback.format_exc())
            return None

    def plot_image(image_bucket, img_idx):
        try:
            data = image_bucket.data
            x = image_bucket.x
            y = image_bucket.y
            band = image_bucket.band
            year_obs = image_bucket.year_obs
            wcs = image_bucket.wcs
            overlay_label = image_bucket.overlay_label
            overlay_ra = image_bucket.overlay_ra
            overlay_dec = image_bucket.overlay_dec

            ax = fig.add_subplot(rows, cols, img_idx, projection=wcs)
            ax.plot(x, y, 'ro', fillstyle='none', markersize=7, markeredgewidth=0.2)
            ax.plot(x, y, 'ro', fillstyle='none', markersize=0.2, markeredgewidth=0.2)
            ax.text(0.04, 0.90, band, color='black', fontsize=1.8, transform=ax.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
            ax.text(0.04, 0.06, year_obs, color='black', fontsize=1.8, transform=ax.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
            ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))

            if overlays:
                ax.scatter(overlay_ra, overlay_dec, transform=ax.get_transform('icrs'), s=1.0,
                           edgecolor=overlay_color, facecolor='none', linewidths=0.2)

            if overlay_labels:
                for i in range(len(overlay_label)):
                    ax.text(overlay_ra[i], overlay_dec[i], overlay_label[i], transform=ax.get_transform('icrs'),
                            color=overlay_label_color, size=1)  # ha='center', va='center',

            vmin, vmax = get_min_max(data)
            ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
            ax.axis('off')
        except Exception:
            print('A problem occurred while plotting an image for object ra={ra}, dec={dec}, band={band}'.format(ra=ra, dec=dec, band=band))
            print(traceback.format_exc())

    def create_lupton_rgb(data):
        vmin, vmax = get_min_max(data)
        stretch = 1 if vmax-vmin == 0 else vmax-vmin
        return make_lupton_rgb(data, data, data, minimum=vmin, stretch=stretch, Q=0)

    def get_min_max(data, lo=5, hi=95):
        med = np.nanmedian(data)
        mad = np.nanmedian(abs(data - med))
        dev = np.nanpercentile(data, hi) - np.nanpercentile(data, lo)
        vmin = med - 2.0 * mad
        vmax = med + 2.0 * dev
        return vmin, vmax

    def get_neowise_image(ra, dec, epoch, band, size):
        download_url = 'http://byw.tools/cutout?ra={ra}&dec={dec}&size={size}&band={band}&epoch={epoch}'
        download_url = download_url.format(ra=ra, dec=dec, size=size, band=band, epoch=epoch)
        try:
            return fits.open(download_file(download_url, cache=cache, show_progress=show_progress, timeout=timeout))
        except Exception:
            return None

    def get_year(mjd):
        time = Time(mjd, scale='utc', format='mjd')
        return time.ymdhms['year']

    def get_epoch(mjd):
        time = Time(mjd, scale='utc', format='mjd')
        year = time.ymdhms['year']
        month = time.ymdhms['month']
        return str(year) + '/' + str(month)

    def create_obj_name(ra, dec, precision=6):
        ra = round(ra, precision)
        dec = round(dec, precision)
        ra_str = str(ra)
        dec_str = str(dec) if dec < 0 else '+' + str(dec)
        return ra_str + dec_str

    def start_file(filename):
        if sys.platform == 'win32':
            os.startfile(filename)
        else:
            opener = 'open' if sys.platform == 'darwin' else 'evince'
            subprocess.call([opener, filename])

    def calculate_magnitude(flux):
        if flux <= 0:
            return np.nan
        return 22.5 - 2.5 * math.log10(flux)

    def box_contains_target(box_center_ra, box_center_dec, target_ra, target_dec, box_size):
        d = 5
        if (target_dec > d and box_center_dec < d) or (target_dec < d and box_center_dec > d):
            return False, 0, 0

        ra = math.radians(target_ra)
        dec = math.radians(target_dec)
        ra0 = math.radians(box_center_ra)
        dec0 = math.radians(box_center_dec)

        delta_ra = math.degrees(abs(ra0 - ra) * math.cos(dec0))
        if d < delta_ra < 360 - d:
            return False, 0, 0

        # World to pixel
        box_center = box_size/2 + 0.5
        cosc = math.sin(dec0) * math.sin(dec) + math.cos(dec0) * math.cos(dec) * math.cos(ra - ra0)
        x = (math.cos(dec) * math.sin(ra - ra0)) / cosc
        y = (math.cos(dec0) * math.sin(dec) - math.sin(dec0) * math.cos(dec) * math.cos(ra - ra0)) / cosc
        scale = 3600 / pixel_scale
        x = math.degrees(x) * -scale
        y = math.degrees(y) * scale
        x += box_center
        y += box_center
        y = box_size - y
        x -= 0.5
        y += 0.5

        """ This is much too slow!
        w = WCS(naxis=2)
        w.wcs.crpix = [box_center, box_center]
        w.wcs.crval = [box_center_ra, box_center_dec]
        w.wcs.cunit = ["deg", "deg"]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cdelt = [-0.000763888888889, 0.000763888888889]
        w.array_shape = [box_size, box_size]
        x, y = w.world_to_pixel(SkyCoord(target_ra*u.deg, target_dec*u.deg))
        """

        # Distance to closest edge
        if x > box_center:
            x = box_size - x

        if y > box_center:
            y = box_size - y

        # Check if box contains target
        match = True
        if np.isnan(x) or np.isnan(y) or x < 0 or y < 0:
            match = False

        return match, x, y

    def find_catalog_entries(file_path, file_number):
        hdul = fits.open(base_url + file_path.replace('./', ''), cache=cache, show_progress=show_progress, timeout=timeout)

        data = hdul[1].data
        hdul.close()

        tot_cat_entries = len(data)

        coords_w1 = []
        coords_w2 = []

        object_number = 0

        print('Number of catalog entries scanned for file ' + file_path + ':')

        for i in range(len(data)):
            row = data[i]

            catalog_ra = row['ra']
            catalog_dec = row['dec']

            match, _, _ = box_contains_target(target_ra, target_dec, catalog_ra, catalog_dec, img_size)

            if match:
                object_number += 1
                object_label = str(file_number) + '.' + str(object_number)

                target_coords = SkyCoord(target_ra*u.deg, target_dec*u.deg)
                catalog_coords = SkyCoord(catalog_ra*u.deg, catalog_dec*u.deg)
                target_dist = target_coords.separation(catalog_coords).arcsec

                band = row['band']

                """
                All fluxes are in Vega nanomaggies (nMgy; Finkbeiner+ 2004AJ....128.2577F), so that the Vega magnitude of a source is given by 22.5-2.5log(flux).
                The absolute calibration is ultimately inherited from AllWISE through the calibration of Meisner+ (2017AJ....154..161M).
                This inheritance depends on details of the PSF normalization at large radii, which is uncertain.
                To improve the agreement between unWISE and AllWISE fluxes, we recommend subtracting 4mmag from unWISE W1 and 32mmag from unWISE W2 fluxes.
                """
                # Calculate Vega magnitude from flux
                mag_corr = 4 if band == 1 else 32
                mag = calculate_magnitude(row['flux'] - mag_corr)
                if np.isnan(mag):
                    dmag = np.nan
                else:
                    dmag = calculate_magnitude(row['dflux'] - mag_corr)
                    dmag = dmag/1000

                result_table.add_row((
                    object_label,
                    target_dist,
                    row['x'],
                    row['y'],
                    row['flux'],
                    row['dx'],
                    row['dy'],
                    row['dflux'],
                    row['qf'],
                    row['rchi2'],
                    row['fracflux'],
                    row['fluxlbs'],
                    row['dfluxlbs'],
                    row['fwhm'],
                    row['spread_model'],
                    row['dspread_model'],
                    row['fluxiso'],
                    row['xiso'],
                    row['yiso'],
                    row['sky'],
                    row['ra'],
                    row['dec'],
                    row['coadd_id'],
                    row['band'],
                    row['unwise_detid'],
                    row['nm'],
                    row['primary'],
                    row['flags_unwise'],
                    row['flags_info'],
                    row['EPOCH'],
                    row['FORWARD'],
                    row['MJDMIN'],
                    row['MJDMAX'],
                    row['MJDMEAN'],
                    mag,
                    dmag
                ))

                if band == 1:
                    coords_w1.append((object_label, row['ra'], row['dec']))
                else:
                    coords_w2.append((object_label, row['ra'], row['dec']))

            if i > 0 and i % 10000 == 0:
                print(str(i) + '/' + str(tot_cat_entries))

        print(str(tot_cat_entries) + '/' + str(tot_cat_entries))

        return coords_w1, coords_w2

    # ---------------------------------------
    # Code for search_by_coordinates function
    # ---------------------------------------
    # base_url = 'https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/'
    base_url = 'http://unwise.me/data/neo7/untimely-catalog/'  # faster, no timeouts!

    os.chdir(directory)

    pixel_scale = 2.75

    img_size = math.ceil(box_size / pixel_scale)

    index_file = 'untimely_index-neo7.fits'
    if exists(index_file):
        hdul = fits.open(index_file)
    else:
        hdul = fits.open(base_url + index_file + '.gz', cache=cache, show_progress=show_progress, timeout=timeout)
        hdul.writeto(index_file)

    data = hdul[1].data
    hdul.close()

    tot_entries = len(data)

    file_series = []
    tile_catalog_files = None

    prev_coadd_id = ''

    print('Number of index entries scanned:')

    for i in range(len(data)):
        row = data[i]

        # band = row['BAND']
        # epoch = row['EPOCH']
        coadd_id = row['COADD_ID']
        catalog_filename = row['CATALOG_FILENAME']
        tile_center_ra = row['RA']
        tile_center_dec = row['DEC']

        # print(band, coadd_id, epoch, catalog_filename, tile_center_ra, tile_center_dec)

        match, x, y = box_contains_target(tile_center_ra, tile_center_dec, target_ra, target_dec, 2048)

        if match:
            if coadd_id != prev_coadd_id:
                if tile_catalog_files:
                    file_series.append(tile_catalog_files)
                tile_catalog_files = []
                xy = x * y
                tile_catalog_files.append(xy)

            tile_catalog_files.append(catalog_filename)

        prev_coadd_id = coadd_id

        if i > 0 and i % 10000 == 0:
            print(str(i) + '/' + str(tot_entries))

    print(str(tot_entries) + '/' + str(tot_entries))

    file_series.append(tile_catalog_files)

    file_series.sort(key=lambda x: x[0], reverse=True)

    print(file_series)

    w1_overlays = []
    w2_overlays = []

    if len(file_series) > 0:
        catalog_files = file_series[0]

        result_table = Table(names=(
            'object_label',
            'target_dist',
            'x',
            'y',
            'flux',
            'dx',
            'dy',
            'dflux',
            'qf',
            'rchi2',
            'fracflux',
            'fluxlbs',
            'dfluxlbs',
            'fwhm',
            'spread_model',
            'dspread_model',
            'fluxiso',
            'xiso',
            'yiso',
            'sky',
            'ra',
            'dec',
            'coadd_id',
            'band',
            'unwise_detid',
            'nm',
            'primary',
            'flags_unwise',
            'flags_info',
            'EPOCH',
            'FORWARD',
            'MJDMIN',
            'MJDMAX',
            'MJDMEAN',
            'mag',
            'dmag'
        ), dtype=('S', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
                  'f', 'f', 'f', 'f', 'S', 'i', 'S', 'i', 'i', 'i', 'i', 'i', 'i', 'f', 'f', 'f', 'f', 'f')
        )

        for i in range(len(catalog_files)):
            if i == 0:
                continue
            coords_w1, coords_w2 = find_catalog_entries(catalog_files[i], i)
            if len(coords_w1) > 0:
                w1_overlays.append(coords_w1)
            if len(coords_w2) > 0:
                w2_overlays.append(coords_w2)

        # result_table.sort('target_dist')
        # result_table.pprint_all()

        if save_result_table:
            result_file_name = 'unTimely_Catalog_search results_' + create_obj_name(target_ra, target_dec) + '.' + result_table_extension
            result_table.write(result_file_name, format=result_table_format, overwrite=True)

        if show_result_table_in_browser:
            result_table.show_in_browser(jsviewer=True)

    if finder_charts:
        # Prepare for plotting the unWISE images
        fig = plt.figure()
        fig.set_figheight(5)
        fig.set_figwidth(5)
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.5)

        cols = 6
        rows = 12

        ra = target_ra
        dec = target_dec

        # --------------------------------------------
        # Create W1 images with corresponding overlays
        # --------------------------------------------

        # Collect W1 images
        images = []
        scans = []
        imageW1 = get_neowise_image(ra, dec, epoch=0, band=1, size=img_size)
        if imageW1:
            hdu = imageW1[0]
            header = hdu.header
            meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
            time = Time(meanmjd, scale='utc', format='mjd')
            prev_year = time.ymdhms['year']
            prev_month = time.ymdhms['month']
            scan_change = False
            print('Downloading W1 images ...')
            print('Ascending scans for ' + str(prev_year) + ':')

            for i in range(0, 100, 1):
                try:
                    imageW1 = get_neowise_image(ra, dec, epoch=i, band=1, size=img_size)
                    if not imageW1:
                        break

                    hdu = imageW1[0]
                    header = hdu.header
                    meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
                    time = Time(meanmjd, scale='utc', format='mjd')
                    year = time.ymdhms['year']
                    month = time.ymdhms['month']
                    if year in (2011, 2013):
                        continue

                    if year == prev_year:
                        if month - prev_month > 4 and not scan_change:
                            scan_change = True
                            images.append(scans)
                            scans = []
                            print('Descending scans for ' + str(year) + ':')
                    else:
                        scan_change = False
                        images.append(scans)
                        scans = []
                        print('Ascending scans for ' + str(year) + ':')

                    scans.append((hdu, meanmjd))
                    print('  ' + str(year) + '/' + str(month))

                    prev_year = year
                    prev_month = month
                except Exception:
                    print('A problem occurred while creating WISE time series for object ra={ra}, dec={dec}, epoch={epoch}, band=1'
                          .format(ra=ra, dec=dec, epoch=i))
                    print(traceback.format_exc())

            images.append(scans)

        # Create W1 epoch images
        epochs = []
        for scans in images:
            size = len(scans) + 1
            hdu = scans[0][0]
            header = hdu.header
            data = hdu.data
            mjd = scans[0][1]
            for scan in scans:
                data += scan[0].data
                mjd += scan[1]
            hdu = fits.PrimaryHDU(data=data/size, header=header)
            epochs.append((hdu, get_epoch(mjd/size)))
        images = epochs

        # Process W1 images
        w1_images = []
        for i in range(min(len(images), len(w1_overlays))):
            image = images[i]
            w1, x, y, wcs = process_image_data(image[0])
            overlay_label = [coords[0] for coords in w1_overlays[i]]
            overlay_ra = [coords[1] for coords in w1_overlays[i]]
            overlay_dec = [coords[2] for coords in w1_overlays[i]]
            w1_images.append(ImageBucket(w1, x, y, 'W1', image[1], wcs, overlay_label, overlay_ra, overlay_dec))

        # Plot W1 images
        img_idx = 0
        for image_bucket in w1_images:
            img_idx += 1
            plot_image(image_bucket, img_idx)

        r = img_idx % cols
        if r > 0:
            img_idx += cols - r

        # --------------------------------------------
        # Create W2 images with corresponding overlays
        # --------------------------------------------

        # Collect W2 images
        images = []
        scans = []
        imageW2 = get_neowise_image(ra, dec, epoch=0, band=2, size=img_size)
        if imageW2:
            hdu = imageW2[0]
            header = hdu.header
            meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
            time = Time(meanmjd, scale='utc', format='mjd')
            prev_year = time.ymdhms['year']
            prev_month = time.ymdhms['month']
            scan_change = False
            print('Downloading W2 images ...')
            print('Ascending scans for ' + str(prev_year) + ':')

            for i in range(0, 100, 1):
                try:
                    imageW2 = get_neowise_image(ra, dec, epoch=i, band=2, size=img_size)
                    if not imageW2:
                        break

                    hdu = imageW2[0]
                    header = hdu.header
                    meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
                    time = Time(meanmjd, scale='utc', format='mjd')
                    year = time.ymdhms['year']
                    month = time.ymdhms['month']
                    if year in (2011, 2013):
                        continue

                    if year == prev_year:
                        if month - prev_month > 4 and not scan_change:
                            scan_change = True
                            images.append(scans)
                            scans = []
                            print('Descending scans for ' + str(year) + ':')
                    else:
                        scan_change = False
                        images.append(scans)
                        scans = []
                        print('Ascending scans for ' + str(year) + ':')

                    scans.append((hdu, meanmjd))
                    print('  ' + str(year) + '/' + str(month))

                    prev_year = year
                    prev_month = month
                except Exception:
                    print('A problem occurred while creating WISE time series for object ra={ra}, dec={dec}, epoch={epoch}, band=2'
                          .format(ra=ra, dec=dec, epoch=i))
                    print(traceback.format_exc())

            images.append(scans)

        # Create W2 epoch images
        epochs = []
        for scans in images:
            size = len(scans) + 1
            hdu = scans[0][0]
            header = hdu.header
            data = hdu.data
            mjd = scans[0][1]
            for scan in scans:
                data += scan[0].data
                mjd += scan[1]
            hdu = fits.PrimaryHDU(data=data/size, header=header)
            epochs.append((hdu, get_epoch(mjd/size)))
        images = epochs

        # Process W2 images
        w2_images = []
        for i in range(min(len(images), len(w2_overlays))):
            image = images[i]
            w2, x, y, wcs = process_image_data(image[0])
            overlay_label = [coords[0] for coords in w2_overlays[i]]
            overlay_ra = [coords[1] for coords in w2_overlays[i]]
            overlay_dec = [coords[2] for coords in w2_overlays[i]]
            w2_images.append(ImageBucket(w2, x, y, 'W2', image[1], wcs, overlay_label, overlay_ra, overlay_dec))

        # Plot W2 images
        for image_bucket in w2_images:
            img_idx += 1
            plot_image(image_bucket, img_idx)

        filename = 'unTimely_Catalog_finder_charts_' + create_obj_name(ra, dec) + '.' + finder_charts_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=finder_charts_format)
        plt.close()

        if open_finder_charts:
            start_file(filename)

        if animated_gif:
            w1_images_plus_overlays = []
            for i in range(min(len(w1_images), len(w1_overlays))):
                w1_images_plus_overlays.append((w1_images[i], w1_overlays[i]))

            w2_images_plus_overlays = []
            for i in range(min(len(w2_images), len(w2_overlays))):
                w2_images_plus_overlays.append((w2_images[i], w2_overlays[i]))

            w1_images = w1_images_plus_overlays
            w2_images = w2_images_plus_overlays

            w1_reordred = []
            w2_reordred = []

            # Separate scan directions
            if scan_dir_mode == SEPARATE_SCAN:
                for i in range(0, len(w1_images), 2):
                    w1_reordred.append(w1_images[i])
                for i in range(1, len(w1_images), 2):
                    w1_reordred.append(w1_images[i])

                for i in range(0, len(w2_images), 2):
                    w2_reordred.append(w2_images[i])
                for i in range(1, len(w2_images), 2):
                    w2_reordred.append(w2_images[i])

                w1_images = w1_reordred
                w2_images = w2_reordred

            # Merge scan directions
            if scan_dir_mode == MERGE_SCAN:
                for i in range(0, len(w1_images), 2):
                    w1_asc = w1_images[i][0].data
                    w1_asc_y = w1_images[i][0].year_obs
                    w1_asc_o = w1_images[i][1]
                    try:
                        w1_des = w1_images[i+1][0].data
                        w1_des_o = w1_images[i+1][1]
                    except IndexError:
                        w1_des = w1_asc
                        w1_des_o = w1_asc_o
                    w1_images[i][0].data = (w1_asc+w1_des)/2
                    w1_images[i][0].year_obs = w1_asc_y[0:4]
                    w1_asc_o.extend(w1_des_o)
                    w1_reordred.append(w1_images[i])

                for i in range(0, len(w2_images), 2):
                    w2_asc = w2_images[i][0].data
                    w2_asc_y = w2_images[i][0].year_obs
                    w2_asc_o = w2_images[i][1]
                    try:
                        w2_des = w2_images[i+1][0].data
                        w2_des_o = w2_images[i+1][1]
                    except IndexError:
                        w2_des = w2_asc
                        w2_des_o = w2_asc_o
                    w2_images[i][0].data = (w2_asc+w2_des)/2
                    w2_images[i][0].year_obs = w2_asc_y[0:4]
                    w2_asc_o.extend(w2_des_o)
                    w2_reordred.append(w2_images[i])

                w1_images = w1_reordred
                w2_images = w2_reordred

            # Draw settings
            zoom = 10
            duration = 200
            stroke_width = 3
            circle_radius = 50
            point_radius = 2
            overlay_radius = 5
            red = (255, 0, 0)
            green = (0, 128, 0)

            # Create animated GIF - W1 with overlays
            images = []

            for i in range(len(w1_images)):
                w1_bucket = w1_images[i][0]
                w1 = w1_bucket.data
                year_obs = w1_bucket.year_obs
                wcs = w1_bucket.wcs

                # Create RGB image
                rgb_image = create_rgb_image(w1, w1, w1)
                rgb_image.info['duration'] = duration

                # Draw a crosshair
                w, h = rgb_image.size
                cx = w/2
                cy = h/2
                draw = ImageDraw.Draw(rgb_image)
                draw.arc((cx-circle_radius, cy-circle_radius, cx+circle_radius, cy+circle_radius),
                         start=0, end=360, fill=red, width=stroke_width)
                draw.arc((cx-point_radius, cy-point_radius, cx+point_radius, cy+point_radius),
                         start=0, end=360, fill=red, width=stroke_width)

                # Draw catalog overlays
                w1_overlays = w1_images[i][1]
                overlay_ra = [coords[1] for coords in w1_overlays]
                overlay_dec = [coords[2] for coords in w1_overlays]
                for i in range(len(overlay_ra)):
                    world = SkyCoord(overlay_ra[i]*u.deg, overlay_dec[i]*u.deg)
                    x, y = wcs.world_to_pixel(world)
                    x += 0.5
                    y += 0.5
                    x *= zoom
                    y *= zoom
                    y = h - y
                    draw.arc((x-overlay_radius, y-overlay_radius, x+overlay_radius, y+overlay_radius),
                             start=0, end=360, fill=green, width=stroke_width)

                # Draw epoch text
                draw.text((10, 10), 'W1 ' + year_obs, red)

                images.append(rgb_image)

            filename = 'Animated_time_series_w1_' + create_obj_name(ra, dec) + '.gif'
            images[0].save(filename, save_all=True, append_images=images[1:], loop=0)

            if open_finder_charts:
                start_file(filename)

            # Create animated GIF - W2 with overlays
            images = []

            for i in range(len(w2_images)):
                w2_bucket = w2_images[i][0]
                w2 = w2_bucket.data
                year_obs = w2_bucket.year_obs
                wcs = w2_bucket.wcs

                # Create RGB image
                rgb_image = create_rgb_image(w2, w2, w2)
                rgb_image.info['duration'] = duration

                # Draw a crosshair
                w, h = rgb_image.size
                cx = w/2
                cy = h/2
                draw = ImageDraw.Draw(rgb_image)
                draw.arc((cx-circle_radius, cy-circle_radius, cx+circle_radius, cy+circle_radius),
                         start=0, end=360, fill=red, width=stroke_width)
                draw.arc((cx-point_radius, cy-point_radius, cx+point_radius, cy+point_radius),
                         start=0, end=360, fill=red, width=stroke_width)

                # Draw catalog overlays
                w2_overlays = w2_images[i][1]
                overlay_ra = [coords[1] for coords in w2_overlays]
                overlay_dec = [coords[2] for coords in w2_overlays]
                for i in range(len(overlay_ra)):
                    world = SkyCoord(overlay_ra[i]*u.deg, overlay_dec[i]*u.deg)
                    x, y = wcs.world_to_pixel(world)
                    x += 0.5
                    y += 0.5
                    x *= zoom
                    y *= zoom
                    y = h - y
                    draw.arc((x-overlay_radius, y-overlay_radius, x+overlay_radius, y+overlay_radius),
                             start=0, end=360, fill=green, width=stroke_width)

                # Draw epoch text
                draw.text((10, 10), 'W2 ' + year_obs, red)

                images.append(rgb_image)

            filename = 'Animated_time_series_w2_' + create_obj_name(ra, dec) + '.gif'
            images[0].save(filename, save_all=True, append_images=images[1:], loop=0)

            if open_finder_charts:
                start_file(filename)

            # Create animated GIF - W1+W2 without overlays
            images = []
            for i in range(len(w1_images)):
                w1_bucket = w1_images[i][0]
                w2_bucket = w2_images[i][0]
                w1 = w1_bucket.data
                w2 = w2_bucket.data
                year_obs = w1_bucket.year_obs

                # Create RGB image
                rgb_image = create_rgb_image(w1, (w1+w2)/2, w2)
                rgb_image.info['duration'] = duration

                # Draw a crosshair
                w, h = rgb_image.size
                cx = w/2
                cy = h/2
                draw = ImageDraw.Draw(rgb_image)
                draw.arc((cx-circle_radius, cy-circle_radius, cx+circle_radius, cy+circle_radius),
                         start=0, end=360, fill=red, width=stroke_width)
                draw.arc((cx-point_radius, cy-point_radius, cx+point_radius, cy+point_radius),
                         start=0, end=360, fill=red, width=stroke_width)

                # Draw epoch text
                draw.text((10, 10), 'W1+W2 ' + year_obs, red)

                images.append(rgb_image)

            filename = 'Animated_time_series_' + create_obj_name(ra, dec) + '.gif'
            images[0].save(filename, save_all=True, append_images=images[1:], loop=0)

            if open_finder_charts:
                start_file(filename)
