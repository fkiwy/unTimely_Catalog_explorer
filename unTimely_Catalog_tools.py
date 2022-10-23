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


class unTimelyCatalogExplorer:

    # Scan direction modes for image blinks
    ALTERNATE_SCAN = 0
    SEPARATE_SCAN = 1
    MERGE_SCAN = 2

    def __init__(self, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
                 catalog_base_url='https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/',
                 catalog_index_file='untimely_index-neo7.fits'):
        """
        Creates an unTimelyCatalogExplorer instance with the given parameters

        Parameters
        ----------
        directory : str, optional
            Directory where the finder charts should be saved. The default is tempfile.gettempdir().
        cache : bool, optional
            Whether to cache the downloaded files. The default is True.
        show_progress : bool, optional
            Whether to show the file download progress. The default is True.
        timeout : int, optional
            Timeout for remote requests in seconds. The default is 300.
        catalog_base_url : str, optional
            Base URL to access the unTimely Catalog. The default is 'https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/'.
        catalog_index_file : str, optional
            Catalog index file name. The default is 'untimely_index-neo7.fits'.

        Returns
        -------
        An unTimelyCatalogExplorer instance.

        """
        self.catalog_base_url = catalog_base_url
        self.catalog_index_file = catalog_index_file
        self.cache = cache
        self.show_progress = show_progress
        self.timeout = timeout
        self.result_table = None
        self.w1_images = None
        self.w2_images = None
        self.w1_overlays = None
        self.w2_overlays = None
        self.pixel_scale = 2.75
        os.chdir(directory)

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

    def process_image_data(self, hdu, ra, dec, box_size):
        data = hdu.data
        wcs, shape = find_optimal_celestial_wcs([hdu], frame='icrs')
        data, _ = reproject_interp(hdu, wcs, shape_out=shape)
        position = SkyCoord(ra*u.deg, dec*u.deg)
        cutout = Cutout2D(data, position, box_size*u.arcsec, wcs=wcs, mode='partial')
        data = cutout.data
        wcs = cutout.wcs
        x, y = wcs.world_to_pixel(position)
        return data, x, y, wcs

    def create_rgb_image(self, r, g, b, image_contrast, image_zoom):
        xmax, ymax = g.shape
        vmin, vmax = self.get_min_max(g, image_contrast)
        image = Image.fromarray(make_lupton_rgb(r, g, b, minimum=vmin, stretch=vmax-vmin, Q=0)).resize((image_zoom*xmax, image_zoom*ymax), Image.NONE)
        image = image.transpose(Image.FLIP_TOP_BOTTOM)
        image = ImageOps.invert(image)
        return image

    def plot_image(self, image_bucket, fig, rows, cols, img_idx, overlays, overlay_color, overlay_labels, overlay_label_color, image_contrast):
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

        if overlays or overlay_labels:
            ax.scatter(overlay_ra, overlay_dec, transform=ax.get_transform('icrs'), s=1.0,
                       edgecolor=overlay_color, facecolor='none', linewidths=0.2)

        if overlay_labels:
            for i in range(len(overlay_label)):
                ax.text(overlay_ra[i], overlay_dec[i], overlay_label[i], transform=ax.get_transform('icrs'),
                        color=overlay_label_color, size=1)

        vmin, vmax = self.get_min_max(data, image_contrast)
        ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
        ax.axis('off')

    def get_min_max(self, data, image_contrast):
        lo = image_contrast
        hi = 100-image_contrast
        med = np.nanmedian(data)
        mad = np.nanmedian(abs(data - med))
        dev = np.nanpercentile(data, hi) - np.nanpercentile(data, lo)
        vmin = med - 2.0 * mad
        vmax = med + 2.0 * dev
        return vmin, vmax

    def get_neowise_image(self, ra, dec, epoch, band, size):
        download_url = 'http://byw.tools/cutout?ra={ra}&dec={dec}&size={size}&band={band}&epoch={epoch}'
        download_url = download_url.format(ra=ra, dec=dec, size=size, band=band, epoch=epoch)
        try:
            return fits.open(download_file(download_url, cache=self.cache, show_progress=self.show_progress, timeout=self.timeout))
        except Exception:
            return None

    def get_epoch(self, mjd):
        time = Time(mjd, scale='utc', format='mjd')
        year = time.ymdhms['year']
        month = time.ymdhms['month']
        return str(year) + '/' + str(month)

    def create_obj_name(self, ra, dec, precision=6):
        ra = round(ra, precision)
        dec = round(dec, precision)
        ra_str = str(ra)
        dec_str = str(dec) if dec < 0 else '+' + str(dec)
        return ra_str + dec_str

    def create_j_designation(self, ra, dec):
        return 'J' + SkyCoord(ra*u.degree, dec*u.degree).to_string('hmsdms', sep='', precision=2).replace(' ', '')

    def start_file(self, filename):
        if sys.platform == 'win32':
            os.startfile(filename)
        else:
            opener = 'open' if sys.platform == 'darwin' else 'evince'
            subprocess.call([opener, filename])

    def calculate_magnitude(self, flux):
        if flux <= 0:
            return np.nan
        return 22.5 - 2.5 * math.log10(flux)

    def box_contains_target(self, box_center_ra, box_center_dec, target_ra, target_dec, box_size):
        # Pre-filtering on ra and dec to avoid cases not well handled by the world to pixel solution
        d = 1  # Tile size in degrees: 4048 * 2.75 / 3600 = 1.564 deg (1.564 / 2 = 0.782 ~ 1 deg)
        if abs(box_center_dec - target_dec) > d:
            return False, 0, 0
        if -d < target_dec < d and d < abs(box_center_ra - target_ra) < 360 - d:
            return False, 0, 0

        # World to pixel
        ra = math.radians(target_ra)
        dec = math.radians(target_dec)
        ra0 = math.radians(box_center_ra)
        dec0 = math.radians(box_center_dec)
        cosc = math.sin(dec0) * math.sin(dec) + math.cos(dec0) * math.cos(dec) * math.cos(ra - ra0)
        x = (math.cos(dec) * math.sin(ra - ra0)) / cosc
        y = (math.cos(dec0) * math.sin(dec) - math.sin(dec0) * math.cos(dec) * math.cos(ra - ra0)) / cosc
        scale = 3600 / self.pixel_scale
        x = math.degrees(x) * -scale
        y = math.degrees(y) * scale
        x += box_size/2
        y += box_size/2
        y = box_size - y

        """ Too slow!
        w = WCS(naxis=2)
        w.wcs.crpix = [box_size/2, box_size/2]
        w.wcs.crval = [box_center_ra, box_center_dec]
        w.wcs.cunit = ['deg', 'deg']
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        w.wcs.cdelt = [-0.000763888888889, 0.000763888888889]
        w.array_shape = [box_size, box_size]
        x, y = w.world_to_pixel(SkyCoord(target_ra*u.deg, target_dec*u.deg))
        """

        # Distance to closest edge
        if x > box_size/2:
            x = box_size - x

        if y > box_size/2:
            y = box_size - y

        # Check if box contains target
        match = True
        if np.isnan(x) or np.isnan(y) or x < 0 or y < 0:
            match = False

        return match, x, y

    def find_catalog_entries(self, file_path, file_number, target_ra, target_dec, img_size, result_table):
        hdul = fits.open(self.catalog_base_url + file_path.replace('./', ''), cache=self.cache, show_progress=self.show_progress, timeout=self.timeout)

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

            match, _, _ = self.box_contains_target(target_ra, target_dec, catalog_ra, catalog_dec, img_size)

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
                mag = self.calculate_magnitude(row['flux'] - mag_corr)
                if np.isnan(mag):
                    dmag = np.nan
                else:
                    dmag = self.calculate_magnitude(row['dflux'] - mag_corr)
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

    def search_by_coordinates(self, target_ra, target_dec, box_size=100, show_result_table_in_browser=False, save_result_table=True,
                              result_table_format='ascii', result_table_extension='dat'):
        """
        Search the catalog by coordinates (box search).

        Parameters
        ----------
        ra : float
            Right ascension in decimal degrees.
        dec : float
            Declination in decimal degrees.
        box_size : int, optional
            Image size in arcseconds. The default is 100.
        show_result_table_in_browser : bool, optional
            Whether to show the result table in your browser (columns can be sorted). The default is False.
        save_result_table : bool, optional
            Whether to save the result table to the directory specified in the constructor ``unTimelyCatalogExplorer(directory=)``. The default is True.
        result_table_format : str, optional
            Result table output format. The default is 'ascii'.
        result_table_extension : str, optional
            Result table file extension. The default is 'dat'.

        Returns
        -------
        result_table : astropy.table.table.Table
            Result table containing the catalog entries located within a field of view of the specified size at the given coordinates.

        """
        self.target_ra = target_ra
        self.target_dec = target_dec
        self.box_size = box_size
        self.img_size = int(round(box_size / self.pixel_scale, 0))

        if exists(self.catalog_index_file):
            hdul = fits.open(self.catalog_index_file)
        else:
            hdul = fits.open(self.catalog_base_url + self.catalog_index_file + '.gz', cache=self.cache, show_progress=self.show_progress, timeout=self.timeout)
            hdul.writeto(self.catalog_index_file)

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

            match, x, y = self.box_contains_target(tile_center_ra, tile_center_dec, target_ra, target_dec, 2048)

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

        self.w1_overlays = []
        self.w2_overlays = []

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
                'epoch',
                'forward',
                'mjdmin',
                'mjdmax',
                'mjdmean',
                'mag',
                'dmag'
            ), dtype=('S', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
                      'f', 'f', 'f', 'f', 'S', 'i', 'S', 'i', 'i', 'i', 'i', 'i', 'i', 'f', 'f', 'f', 'f', 'f'),
                units=('', 'arcsec', 'pix', 'pix', 'nMgy', 'pix', 'pix', 'nMgy', '', '', '', 'nMgy', 'nMgy', 'pix',
                       '', '', '', '', '', 'nMgy', 'deg', 'deg', '', '', '', '', '', '', '', '', '', 'd', 'd', 'd',
                       'mag', 'mag'),
                descriptions=('Unique object label for a specific result set to retrieve the corresponding source on the finder charts',
                              'Angular distance to the target coordinates',
                              'x coordinate',
                              'y coordinate',
                              'Vega flux',
                              'x uncertainty',
                              'y uncertainty',
                              'formal flux uncertainty',
                              'PSF-weighted fraction of good pixels',
                              'PSF-weighted average chi2',
                              'PSF-weighted fraction of flux from this source',
                              'FWHM of PSF at source location',
                              'local-background-subtracted flux',
                              'formal fluxlbs uncertainty',
                              'SExtractor-like source size parameter',
                              'uncertainty in spread_model',
                              'flux derived from linear least squares fit to neighbor-subtracted image; significant difference from ordinary flux indicates a convergence issue',
                              'x coordinate derived from linear least squares fit to neighbor-subtracted image; significant difference from ordinary x indicates a convergence issue',
                              'y coordinate derived from linear least squares fit to neighbor-subtracted image; significant difference from ordinary y indicates a convergence issue',
                              'residual sky at source location',
                              'R.A.',
                              'decl.',
                              'unWISE/AllWISE coadd_id of source',
                              '1 for W1, 2 for W2',
                              'detection ID, unique in catalog',
                              'number of images in coadd at source',
                              'source located in primary region of coadd',
                              'unWISE flags at source location',
                              'additional flags at source location',
                              'unWISE epoch number',
                              "boolean, were input frames acquired pointing forward (1) or backward (0) along Earth's orbit",
                              'MJD value of earliest contributing exposure',
                              'MJD value of latest contributing exposure',
                              'mean of MJDMIN and MJDMAX',
                              'Vega magnitude given by 22.5-2.5log(flux)',
                              'magnitude uncertainty'
                              )
            )

            for i in range(len(catalog_files)):
                if i == 0:
                    continue
                coords_w1, coords_w2 = self.find_catalog_entries(catalog_files[i], i, target_ra, target_dec, self.img_size, result_table)
                if len(coords_w1) > 0:
                    self.w1_overlays.append(coords_w1)
                if len(coords_w2) > 0:
                    self.w2_overlays.append(coords_w2)

            # result_table.sort('target_dist')
            # result_table.pprint_all()

            self.result_table = result_table

            if save_result_table:
                result_file_name = 'unTimely_Catalog_search results_' + self.create_obj_name(target_ra, target_dec) + '.' + result_table_extension
                result_table.write(result_file_name, format=result_table_format, overwrite=True)

            if show_result_table_in_browser:
                result_table.show_in_browser(jsviewer=True)

            return result_table

    def create_finder_charts(self, overlays=True, overlay_color='green', overlay_labels=False, overlay_label_color='red',
                             image_contrast=3, open_file=False, file_format='pdf'):
        """
        Create finder charts for W1 and W2 at each epoch with overplotted catalog positions (overlays)

        Parameters
        ----------
        overlays : bool, optional
            Whether to plot W1 and W2 catalog positions on the finder charts (overlays). The default is True.
        overlay_color : str, optional
            Overlay color. The default is 'green'.
        overlay_labels : bool, optional
            Whether to plot catalog entry labels on the finder charts (tied to overlay circles). The default is False.
        overlay_label_color : str, optional
            Label color. The default is 'red'.
        image_contrast : int, optional
            Contrast of W1 and W2 images. The default is 3.
        open_file : bool, optional
            Whether to open the saved finder charts automatically. The default is False.
        file_format : str, optional
            Output file format: pdf, png, eps, etc.. The default is 'pdf'.

        Raises
        ------
        Exception
            If method ``search_by_coordinates`` has not been called first.

        Returns
        -------
        None.

        """
        if self.w1_overlays is None and self.w2_overlays is None:
            raise Exception('Method ``search_by_coordinates`` must be called first.')

        # Figure settings
        fig = plt.figure()
        fig.set_figheight(5)
        fig.set_figwidth(5)
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.5)

        cols = 6
        rows = 12

        ra = self.target_ra
        dec = self.target_dec
        img_size = self.img_size
        w1_overlays = self.w1_overlays
        w2_overlays = self.w2_overlays

        self.image_contrast = image_contrast
        self.open_file = open_file
        self.file_format = file_format

        # Collect W1 images
        images = []
        scans = []
        imageW1 = self.get_neowise_image(ra, dec, epoch=0, band=1, size=self.img_size)
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
                    imageW1 = self.get_neowise_image(ra, dec, epoch=i, band=1, size=self.img_size)
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
            epochs.append((hdu, self.get_epoch(mjd/size)))
        images = epochs

        # Process W1 images
        self.w1_images = []
        for i in range(min(len(images), len(w1_overlays))):
            image = images[i]
            w1, x, y, wcs = self.process_image_data(image[0], ra, dec, self.box_size)
            overlay_label = [coords[0] for coords in w1_overlays[i]]
            overlay_ra = [coords[1] for coords in w1_overlays[i]]
            overlay_dec = [coords[2] for coords in w1_overlays[i]]
            self.w1_images.append(self.ImageBucket(w1, x, y, 'W1', image[1], wcs, overlay_label, overlay_ra, overlay_dec))

        # Plot W1 images
        img_idx = 0
        for image_bucket in self.w1_images:
            img_idx += 1
            self.plot_image(image_bucket, fig, rows, cols, img_idx, overlays, overlay_color, overlay_labels, overlay_label_color, image_contrast)

        r = img_idx % cols
        if r > 0:
            img_idx += cols - r

        # Collect W2 images
        images = []
        scans = []
        imageW2 = self.get_neowise_image(ra, dec, epoch=0, band=2, size=img_size)
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
                    imageW2 = self.get_neowise_image(ra, dec, epoch=i, band=2, size=img_size)
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
            epochs.append((hdu, self.get_epoch(mjd/size)))
        images = epochs

        # Process W2 images
        self.w2_images = []
        for i in range(min(len(images), len(w2_overlays))):
            image = images[i]
            w2, x, y, wcs = self.process_image_data(image[0], ra, dec, self.box_size)
            overlay_label = [coords[0] for coords in w2_overlays[i]]
            overlay_ra = [coords[1] for coords in w2_overlays[i]]
            overlay_dec = [coords[2] for coords in w2_overlays[i]]
            self.w2_images.append(self.ImageBucket(w2, x, y, 'W2', image[1], wcs, overlay_label, overlay_ra, overlay_dec))

        # Plot W2 images
        for image_bucket in self.w2_images:
            img_idx += 1
            self.plot_image(image_bucket, fig, rows, cols, img_idx, overlays, overlay_color, overlay_labels, overlay_label_color, image_contrast)

        # Info text
        coords = SkyCoord(ra*u.deg, dec*u.deg)
        info_idx = math.ceil(img_idx / cols) * cols + 1
        fontsize = 2.2
        ax = fig.add_subplot(rows, cols, info_idx)
        ax.text(0.05, 0.70, r'$\alpha$ = ' + str(round(coords.ra.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.text(0.05, 0.55, r'$\delta$ = ' + str(round(coords.dec.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.text(0.05, 0.40, '$l$ = ' + str(round(coords.galactic.l.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.text(0.05, 0.25, '$b$ = ' + str(round(coords.galactic.b.value, 6)), fontsize=fontsize, transform=ax.transAxes)
        ax.axis('off')

        # Info text cont'd
        hmsdms = coords.to_string('hmsdms', sep=':', precision=2)
        hms = hmsdms[0:11]
        dms = hmsdms[12:24] if dec < 0 else hmsdms[13:24]
        ax = fig.add_subplot(rows, cols, info_idx + 1)
        ax.text(0, 0.72, '(' + hms + ')', fontsize=fontsize, transform=ax.transAxes)
        ax.text(0, 0.57, '(' + dms + ')', fontsize=fontsize, transform=ax.transAxes)
        ax.text(0, 0.42, 'Size = ' + str(int(self.box_size)) + ' arcsec', fontsize=fontsize, transform=ax.transAxes)
        ax.text(0, 0.27, 'North up, East left', fontsize=fontsize, transform=ax.transAxes)
        ax.axis('off')

        filename = 'unTimely_Catalog_finder_charts_' + self.create_obj_name(ra, dec) + '.' + file_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        if open_file:
            self.start_file(filename)

    def create_ligh_curves(self, photometry_radius=5, yticks=None, open_file=None, file_format=None):
        """
        Create light curves using W1 and W2 photometry of all available epochs.

        Parameters
        ----------
        photometry_radius : float, optional
            Radius to search for the photometry used to create the light curves. The default is 5.
        yticks : tuple, optional
            Tuple containing y axis tick values. The default is pyplot's automatic tick allocation.
        open_file : bool, optional
            Whether to open the saved light curves automatically. The default is None (value given by method ``create_finder_charts`` will be used).
        file_format : bool, optional
            Output file format: pdf, png, eps, etc.. The default is None (value given by method ``create_finder_charts`` will be used).

        Raises
        ------
        Exception
            If method ``search_by_coordinates`` has not been called first.

        Returns
        -------
        None.

        """
        if self.result_table is None:
            raise Exception('Method ``search_by_coordinates`` must be called first.')

        ra = self.target_ra
        dec = self.target_dec
        result_table = self.result_table

        if open_file is None:
            open_file = self.open_file
        if file_format is None:
            file_format = self.file_format

        # Filter by target distance
        mask = result_table['target_dist'] <= photometry_radius
        phot_table = result_table[mask]

        if (len(phot_table) == 0):
            print('No photometry found in specified radius (default is 5 arcsec) to create any light curves.')
            return

        # Get W1 photometry
        mask = phot_table['band'] == 1
        phot_table_w1 = phot_table[mask]

        # Get W2 photometry
        mask = phot_table['band'] == 2
        phot_table_w2 = phot_table[mask]

        # Plot light curves
        x1 = Time(phot_table_w1['mjdmean'], format='mjd').jyear
        y1 = phot_table_w1['mag']
        x2 = Time(phot_table_w2['mjdmean'], format='mjd').jyear
        y2 = phot_table_w2['mag']
        plt.figure(figsize=(8, 4))
        plt.title(self.create_j_designation(ra, dec))
        plt.plot(x1, y1, lw=1, linestyle='--', markersize=3, marker='o', label='W1')
        plt.plot(x2, y2, lw=1, linestyle='--', markersize=3, marker='o', label='W2')
        plt.xticks(range(2010, 2021, 1))
        if yticks:
            plt.yticks(yticks)
        plt.xlabel('Year')
        plt.ylabel('Magnitude (mag)')
        plt.legend(loc='best')
        plt.gca().invert_yaxis()
        plt.grid(color='grey', alpha=0.5, linestyle='-.', linewidth=0.2, axis='both')

        # Save light curves plot
        filename = 'unTimely_Catalog_light_curves_' + self.create_obj_name(ra, dec) + '.' + file_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        # Open saved file
        if open_file:
            self.start_file(filename)

    def create_image_blinks(self, blink_duration=300, image_zoom=10, image_contrast=None, scan_dir_mode=ALTERNATE_SCAN, display_blinks=False):
        """
        Create W1 and W2 image blinks with overplotted catalog positions in GIF format.

        Parameters
        ----------
        blink_duration : int, optional
            Duration each image is shown in milliseconds. The default is 200.
        image_zoom : int, optional
            Scaling factor to be applied on W1 and W2 images. The default is 10.
        image_contrast : int, optional
            Contrast of W1 and W2 images. The default is None (value given by method ``create_finder_charts`` will be used).
        scan_dir_mode : int, optional
            Order in which the image epochs are displayed. The default is ALTERNATE_SCAN.
            - ALTERNATE_SCAN : epoch0asc, epoch0desc, epoch1asc, epoch1desc, ...
            - SEPARATE_SCAN : epoch0asc, epoch1asc, ... epoch0desc, epoch1desc, ...
            - MERGE_SCAN : epoch0asc+epoch0desc, epoch1asc+epoch1desc, ...
        display_blinks : bool, optional
            Whether to display the image blinks in your system's media player. The default is False.

        Raises
        ------
        Exception
            If method ``create_finder_charts`` has not been called first.

        Returns
        -------
        None.

        """
        if self.w1_images is None and self.w2_images is None:
            raise Exception('Method ``create_finder_charts`` must be called first.')

        ra = self.target_ra
        dec = self.target_dec
        w1_images = self.w1_images
        w2_images = self.w2_images
        w1_overlays = self.w1_overlays
        w2_overlays = self.w2_overlays

        if image_contrast is None:
            image_contrast = self.image_contrast

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
        if scan_dir_mode == self.SEPARATE_SCAN:
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
        if scan_dir_mode == self.MERGE_SCAN:
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

            if np.all(w1 == 0):
                continue

            year_obs = w1_bucket.year_obs
            wcs = w1_bucket.wcs

            # Create RGB image
            rgb_image = self.create_rgb_image(w1, w1, w1, image_contrast, image_zoom)

            rgb_image.info['duration'] = blink_duration

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
                x *= image_zoom
                y *= image_zoom
                y = h - y
                draw.arc((x-overlay_radius, y-overlay_radius, x+overlay_radius, y+overlay_radius),
                         start=0, end=360, fill=green, width=stroke_width)

            # Draw epoch text
            draw.text((10, 10), 'W1 ' + year_obs, red)

            images.append(rgb_image)

        filename = 'Animated_time_series_w1_' + self.create_obj_name(ra, dec) + '.gif'
        images[0].save(filename, save_all=True, append_images=images[1:], loop=0)

        if display_blinks:
            self.start_file(filename)

        # Create animated GIF - W2 with overlays
        images = []

        for i in range(len(w2_images)):
            w2_bucket = w2_images[i][0]
            w2 = w2_bucket.data

            if np.all(w2 == 0):
                continue

            year_obs = w2_bucket.year_obs
            wcs = w2_bucket.wcs

            # Create RGB image
            rgb_image = self.create_rgb_image(w2, w2, w2, image_contrast, image_zoom)

            rgb_image.info['duration'] = blink_duration

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
                x *= image_zoom
                y *= image_zoom
                y = h - y
                draw.arc((x-overlay_radius, y-overlay_radius, x+overlay_radius, y+overlay_radius),
                         start=0, end=360, fill=green, width=stroke_width)

            # Draw epoch text
            draw.text((10, 10), 'W2 ' + year_obs, red)

            images.append(rgb_image)

        filename = 'Animated_time_series_w2_' + self.create_obj_name(ra, dec) + '.gif'
        images[0].save(filename, save_all=True, append_images=images[1:], loop=0)

        if display_blinks:
            self.start_file(filename)

        # Create animated GIF - W1+W2 without overlays
        images = []
        for i in range(len(w1_images)):
            w1_bucket = w1_images[i][0]
            w2_bucket = w2_images[i][0]
            w1 = w1_bucket.data
            w2 = w2_bucket.data

            if np.all(w1 == 0) or np.all(w2 == 0):
                continue

            year_obs = w1_bucket.year_obs

            # Create RGB image
            rgb_image = self.create_rgb_image(w1, (w1+w2)/2, w2, image_contrast, image_zoom)

            rgb_image.info['duration'] = blink_duration

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

        filename = 'Animated_time_series_' + self.create_obj_name(ra, dec) + '.gif'
        images[0].save(filename, save_all=True, append_images=images[1:], loop=0)

        if display_blinks:
            self.start_file(filename)
