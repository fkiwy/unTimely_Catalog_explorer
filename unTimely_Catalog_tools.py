import os
from os.path import exists
import sys
import math
import time
import certifi
import warnings
import tempfile
import traceback
import subprocess
import multiprocessing
import numpy as np
import pandas as pd
from urllib import parse
from PIL import Image, ImageOps, ImageDraw, ImageFont
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.stats import sigma_clip
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

    def __init__(self, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, allow_insecure=False, suppress_console_output=False, ignore_warnings=True,
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
        allow_insecure : bool, optional
            Whether to allow downloading files over a TLS/SSL connection even when the server
            certificate verification failed. The default is False.
        suppress_console_output : bool, optional
            Whether to suppress all console output except error messages. The default is False.
        catalog_base_url : str, optional
            Base URL to access the unTimely Catalog. The default is 'https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/'.
        catalog_index_file : str, optional
            Catalog index file name. The default is 'untimely_index-neo7.fits'.

        Returns
        -------
        An unTimelyCatalogExplorer instance.

        """
        if ignore_warnings:
            warnings.simplefilter('ignore', category=Warning)
        self.catalog_base_url = catalog_base_url
        self.catalog_index_file = catalog_index_file
        self.cache = cache
        self.show_progress = False if suppress_console_output else show_progress
        self.timeout = timeout
        self.allow_insecure = allow_insecure
        self.suppress_console_output = suppress_console_output
        self.open_file = False
        self.file_format = 'pdf'
        self.result_table = None
        self.w1_images = None
        self.w2_images = None
        self.w1_overlays = None
        self.w2_overlays = None
        self.pixel_scale = 2.75
        self.unwise_flags = {
            0: 'bright star core and wings',
            1: 'PSF-based diffraction spike',
            2: 'optical ghost',
            3: 'first latent',
            4: 'second latent',
            5: 'circular halo',
            6: 'bright star saturation',
            7: 'geometric diffraction spike'
        }
        self.unwise_info_flags = {
            0: 'bright source off coadd edge',
            1: 'in large galaxy in HyperLeda',
            2: 'in M31 or Magellanic Cloud',
            3: 'may contain bright star center',
            4: 'may be affected by saturation',
            5: 'nebulosity may be present',
            6: 'deblending discouraged here',
            7: 'only "sharp" sources here'
        }
        plt.rcParams.update({'font.size': 8, 'font.family': 'Arial'})
        os.chdir(directory)
        certifi.where()

    class ImageBucket:
        def __init__(self, data, x, y, band, year_obs, wcs, overlay_label, overlay_ra=None, overlay_dec=None, forward=None):
            self.data = data
            self.x = x
            self.y = y
            self.band = band
            self.year_obs = year_obs
            self.wcs = wcs
            self.overlay_label = overlay_label
            self.overlay_ra = overlay_ra
            self.overlay_dec = overlay_dec
            self.forward = forward

    def decompose_flags(self, flag, flag_dict):
        flag = int(flag)
        bits = []
        descr = []

        for i in range(0, 64):
            x = 1
            if (flag & (x << i)) > 0:
                bits.append(i)
                descr.append(flag_dict[i])

        bits = ','.join(list(map(str, bits)))
        descr = ','.join(list(map(str, descr)))
        return bits, descr

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
                ax.annotate(overlay_label[i], (overlay_ra[i], overlay_dec[i]), xycoords=ax.get_transform('icrs'),
                            annotation_clip=True, color=overlay_label_color, size=1)

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

    def get_l1b_photometry(self, ra, dec, radius):
        query_url = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'

        payload = {
            'catalog': 'allwise_p3as_mep',
            'spatial': 'cone',
            'objstr': ' '.join([str(ra), str(dec)]),
            'radius': str(radius),
            'radunits': 'arcsec',
            'outfmt': '1',
            'selcols': 'w1mpro_ep,w1sigmpro_ep,w2mpro_ep,w2sigmpro_ep,mjd,qi_fact,saa_sep,moon_masked'
        }
        r = download_file(query_url + parse.urlencode(payload), cache=self.cache, show_progress=self.show_progress, timeout=self.timeout, allow_insecure=self.allow_insecure)
        allwise = Table.read(r, format='ascii')

        payload = {
            'catalog': 'neowiser_p1bs_psd',
            'spatial': 'cone',
            'objstr': ' '.join([str(ra), str(dec)]),
            'radius': str(radius),
            'radunits': 'arcsec',
            'outfmt': '1',
            'selcols': 'w1mpro,w1sigmpro,w2mpro,w2sigmpro,mjd,qi_fact,saa_sep,moon_masked,qual_frame'
        }
        r = download_file(query_url + parse.urlencode(payload), cache=self.cache, show_progress=self.show_progress, timeout=self.timeout, allow_insecure=self.allow_insecure)
        neowise = Table.read(r, format='ascii')

        # Apply quality constraints
        """
        allwise = allwise[
            (allwise['qi_fact'] > 0.9) *
            (allwise['saa_sep'] > 0) *
            (allwise['moon_masked'] == '0000')
        ]

        neowise = neowise[
            (neowise['qi_fact'] > 0.9) *
            (neowise['saa_sep'] > 0) *
            (neowise['moon_masked'] == '00') *
            (neowise['qual_frame'] > 0)
        ]
        """

        if len(allwise) > 0:
            allwise.rename_column('w1mpro_ep', 'w1mpro')
            allwise.rename_column('w1sigmpro_ep', 'w1sigmpro')
            allwise.rename_column('w2mpro_ep', 'w2mpro')
            allwise.rename_column('w2sigmpro_ep', 'w2sigmpro')
            allwise['qual_frame'] = np.nan
            neowise = vstack([allwise, neowise])

        if len(neowise) > 0:
            neowise['mjd'].unit = 'd'
            neowise['year'] = np.round(Time(neowise['mjd'], format='mjd').jyear, 7)
            neowise.sort('mjd')

        return neowise

    def std_error(self, data):
        return np.ma.std(data) / np.ma.sqrt(len(data))

    def get_neowise_image(self, ra, dec, epoch, band, size):
        download_url = 'http://byw.tools/cutout?ra={ra}&dec={dec}&size={size}&band={band}&epoch={epoch}'
        download_url = download_url.format(ra=ra, dec=dec, size=size, band=band, epoch=epoch)
        try:
            return fits.open(download_file(download_url, cache=self.cache, show_progress=self.show_progress, timeout=self.timeout, allow_insecure=self.allow_insecure))
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

    def disable_print(self):
        self.stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def enable_print(self):
        sys.stdout = self.stdout

    def printout(self, *args):
        if not self.suppress_console_output:
            print(' '.join([str(arg) for arg in args]))

    def print_result_table_info(self):
        info_table = Table(names=['Name', 'Type', 'Unit', 'Description'], dtype=['S', 'S', 'S', 'S'])
        info_table['Name'].format = '<'
        info_table['Type'].format = '<'
        info_table['Unit'].format = '<'
        info_table['Description'].format = '<'
        for colname in self.result_table.colnames:
            col = self.result_table[colname]
            info_table.add_row([col.name, str(col.dtype), str(col.unit), col.description])
        info_table.pprint_all()

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

    def find_catalog_entries(self, file_path, file_number, target_ra, target_dec, box_size, cone_radius, nearest_neighbor, epoch, forward):
        if not self.show_progress:
            self.disable_print()
        hdul = fits.open(download_file(self.catalog_base_url + file_path.replace('./', ''), cache=self.cache,
                         show_progress=self.show_progress, timeout=self.timeout, allow_insecure=self.allow_insecure))
        if not self.show_progress:
            self.enable_print()

        data = hdul[1].data
        hdul.close()

        table = Table(data)
        target_coords = SkyCoord([target_ra*u.deg], [target_dec*u.deg])
        catalog_coords = SkyCoord(table['ra'], table['dec'], unit='deg')
        target_dist = target_coords.separation(catalog_coords).arcsec
        table.add_column(target_dist, name='target_dist')
        table.sort('target_dist')

        if cone_radius:
            box_size = 2 * cone_radius / self.pixel_scale
        else:
            box_size = box_size / self.pixel_scale

        coords_w1 = []
        coords_w2 = []

        object_number = 0

        result_table = Table(names=(
            'source_label',
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
            'dmag',
            'flags_unwise_bits',
            'flags_unwise_descr',
            'flags_info_bits',
            'flags_info_descr'
        ), dtype=('S', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
                  'f', 'f', 'f', 'f', 'S', 'i', 'S', 'i', 'i', 'i', 'i', 'i', 'i', 'f', 'f', 'f', 'f', 'f',
                  'S', 'S', 'S', 'S'),
            units=('', 'arcsec', 'pix', 'pix', 'nMgy', 'pix', 'pix', 'nMgy', '', '', '', 'nMgy', 'nMgy', 'pix',
                   '', '', '', '', '', 'nMgy', 'deg', 'deg', '', '', '', '', '', '', '', '', '', 'd', 'd', 'd',
                   'mag', 'mag', '', '', '', ''),
            descriptions=('Unique source label within a specific result set that can be used to retrieve the corresponding source on the finder charts',
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
                          'Vega magnitude given by 22.5-2.5log10(flux)',
                          'magnitude uncertainty',
                          'unWISE flags bits',
                          'unWISE flags description',
                          'info flags bits',
                          'info flags description')
        )

        for row in table:
            catalog_ra = row['ra']
            catalog_dec = row['dec']
            target_dist = row['target_dist']

            if cone_radius and target_dist > cone_radius:
                continue

            match, _, _ = self.box_contains_target(target_ra, target_dec, catalog_ra, catalog_dec, box_size + 2)

            if match:
                object_number += 1
                source_label = str(file_number) + '.' + str(object_number)
                band = row['band']
                """
                From Schlafly et al. 2019:
                    Fluxes and corresponding uncertainties are given in linear flux units, specifically, in Vega nanomaggies (nMgy; Finkbeiner et al. 2004).
                    The corresponding Vega magnitudes are given by 22.5-2.5log10(flux).
                    The agreement between unWISE and AllWISE magnitudes can be improved by subtracting 4 mmag and 32 mmag from W1 and W2.
                """
                # Calculate Vega magnitude from flux
                # flux_corr = 4 if band == 1 else 32
                flux = row['flux']  # - flux_corr
                mag = self.calculate_magnitude(flux)
                if np.isnan(mag):
                    dmag = np.nan
                else:
                    dflux = row['dflux']
                    mag_upper = self.calculate_magnitude(flux - dflux)
                    mag_lower = self.calculate_magnitude(flux + dflux)
                    dmag = (mag_upper - mag_lower) / 2

                flags_unwise_bits, flags_unwise_descr = self.decompose_flags(row['flags_unwise'], self.unwise_flags)
                flags_info_bits, flags_info_descr = self.decompose_flags(row['flags_info'], self.unwise_info_flags)

                result_table.add_row((
                    source_label,
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
                    dmag,
                    flags_unwise_bits,
                    flags_unwise_descr,
                    flags_info_bits,
                    flags_info_descr
                ))

                if band == 1:
                    coords_w1.append((source_label, row['ra'], row['dec']))
                else:
                    coords_w2.append((source_label, row['ra'], row['dec']))

                if nearest_neighbor:
                    break

        return coords_w1, coords_w2, epoch, forward, result_table

    def search_by_coordinates(self, target_ra, target_dec, box_size=100, cone_radius=None, nearest_neighbor=False, show_result_table_in_browser=False,
                              save_result_table=True, result_table_format='ascii', result_table_extension='dat', multi_processing=False):
        """
        Search the catalog by coordinates (box search).

        Parameters
        ----------
        ra : float
            Right ascension in decimal degrees.
        dec : float
            Declination in decimal degrees.
        box_size : int, optional
            Box search size and/or image size in arcseconds. The default is 100.
        cone_radius : int, optional
            Cone search radius in arcseconds. If specified, a cone search will be performed (instead of a box search) around the given coordinates within the given radius.
            However, the value of the ``box_size`` parameter still defines the image size of the finder charts and image blinks.
        nearest_neighbor : bool, optional
            Whether to include only the nearest detection to the target coordinates in the result table. The default is False.
        show_result_table_in_browser : bool, optional
            Whether to show the result table in your browser (columns can be sorted). The default is False.
        save_result_table : bool, optional
            Whether to save the result table to the directory specified in the constructor ``unTimelyCatalogExplorer(directory=)``. The default is True.
        result_table_format : str, optional
            Result table output format. The default is 'ascii'.
        result_table_extension : str, optional
            Result table file extension. The default is 'dat'.
        multi_processing : bool, optional
            Whether to allow multi-processing when downloading and scanning unTimely catalog files (faster but higher CPU usage). The default is False.
            ``unTimelyCatalogExplorer`` methods must be called within following ``if`` statement: ``if __name__ == '__main__':``

        Returns
        -------
        result_table : astropy.table.table.Table
            Result table containing the catalog entries located within a field of view of the specified size at the given coordinates.

        """
        start_time = time.time()

        self.target_ra = target_ra
        self.target_dec = target_dec
        self.box_size = box_size
        self.img_size = round(box_size / self.pixel_scale)
        self.result_table_format = result_table_format
        self.result_table_extension = result_table_extension

        if exists(self.catalog_index_file):
            hdul = fits.open(self.catalog_index_file)
        else:
            hdul = fits.open(download_file(self.catalog_base_url + self.catalog_index_file + '.gz', cache=self.cache,
                             show_progress=self.show_progress, timeout=self.timeout, allow_insecure=self.allow_insecure))
            hdul.writeto(self.catalog_index_file)

        data = hdul[1].data
        hdul.close()

        table = Table(data)

        file_series = []
        tile_catalog_files = None

        prev_coadd_id = table[0]['COADD_ID']

        self.printout('Scanning catalog index file ...')

        for row in table:
            epoch = row['EPOCH']
            forward = row['FORWARD']
            coadd_id = row['COADD_ID']
            catalog_filename = row['CATALOG_FILENAME']
            tile_center_ra = row['RA']
            tile_center_dec = row['DEC']

            match, x, y = self.box_contains_target(tile_center_ra, tile_center_dec, target_ra, target_dec, 2048)

            if match:
                if coadd_id != prev_coadd_id:
                    if tile_catalog_files:
                        file_series.append(tile_catalog_files)
                    tile_catalog_files = []
                    xy = x * y
                    tile_catalog_files.append(xy)

                tile_catalog_files.append((catalog_filename, epoch, forward))

            prev_coadd_id = coadd_id

        file_series.append(tile_catalog_files)

        file_series.sort(key=lambda x: x[0], reverse=True)

        self.w1_overlays = []
        self.w2_overlays = []

        if len(file_series) > 0:
            catalog_files = file_series[0]

            self.printout('Scanning individual catalog files ...')

            if multi_processing:
                tasks = []
                for i in range(1, len(catalog_files)):
                    catalog_filename = catalog_files[i][0]
                    epoch = catalog_files[i][1]
                    forward = catalog_files[i][2]
                    tasks.append((catalog_filename, i, target_ra, target_dec, box_size, cone_radius, nearest_neighbor, epoch, forward))

                # Multi-processing
                # print('Number of CPUs:', multiprocessing.cpu_count())
                pool = multiprocessing.Pool()
                results = pool.starmap(self.find_catalog_entries, tasks)
                pool.close()
                pool.join()

                coords_w1 = [row[0] for row in results]
                coords_w2 = [row[1] for row in results]
                epoch = [row[2] for row in results]
                forward = [row[3] for row in results]
                tables = [row[4] for row in results]

                for i in range(len(coords_w1)):
                    if len(coords_w1[i]) > 0:
                        self.w1_overlays.append((coords_w1[i], epoch[i], forward[i]))
                for i in range(len(coords_w2)):
                    if len(coords_w2[i]) > 0:
                        self.w2_overlays.append((coords_w2[i], epoch[i], forward[i]))

            else:
                tables = []
                for i in range(1, len(catalog_files)):
                    catalog_filename = catalog_files[i][0]
                    epoch = catalog_files[i][1]
                    forward = catalog_files[i][2]
                    self.printout(catalog_filename)
                    coords_w1, coords_w2, epoch, forward, result_table = self.find_catalog_entries(catalog_filename, i, target_ra, target_dec, box_size,
                                                                                                   cone_radius, nearest_neighbor, epoch, forward)
                    if len(coords_w1) > 0:
                        self.w1_overlays.append((coords_w1, epoch, forward))
                    if len(coords_w2) > 0:
                        self.w2_overlays.append((coords_w2, epoch, forward))
                    tables.append(result_table)

            self.result_table = vstack(tables)

            if save_result_table:
                result_file_name = 'unTimely_Catalog_search results_' + self.create_obj_name(target_ra, target_dec) + '.' + result_table_extension
                self.result_table.write(result_file_name, format=result_table_format, overwrite=True)

            if show_result_table_in_browser:
                self.result_table.show_in_browser(jsviewer=True)

            end_time = time.time()
            elapsed_time = end_time - start_time
            self.printout('Execution time: {:.2f} seconds'.format(elapsed_time))

            return self.result_table

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

        self.printout('Creating finder charts ...')

        # Figure settings
        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)
        plt.subplots_adjust(wspace=0, hspace=0.05, right=0.275)

        cols = 6
        rows = 30

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
        for w1_overlay in w1_overlays:
            try:
                coords_w1 = w1_overlay[0]
                epoch = w1_overlay[1]
                forward = w1_overlay[2]
                imageW1 = self.get_neowise_image(ra, dec, epoch, 1, img_size)
                if not imageW1:
                    break
                hdu = imageW1[0]
                header = hdu.header
                meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
                images.append((hdu, self.get_epoch(meanmjd), coords_w1, forward))
            except Exception:
                print('A problem occurred while creating WISE time series for object ra={ra}, dec={dec}, epoch={epoch}, band=1'
                      .format(ra=ra, dec=dec, epoch=epoch))
                print(traceback.format_exc())

        # Process W1 images
        self.w1_images = []
        for image in images:
            coords_w1 = image[2]
            forward = image[3]
            w1, x, y, wcs = self.process_image_data(image[0], ra, dec, self.box_size)
            overlay_label = [coords[0] for coords in coords_w1]
            overlay_ra = [coords[1] for coords in coords_w1]
            overlay_dec = [coords[2] for coords in coords_w1]
            self.w1_images.append(self.ImageBucket(w1, x, y, 'W1', image[1], wcs, overlay_label, overlay_ra, overlay_dec, forward))

        # Plot W1 images
        img_idx = 0
        for image_bucket in self.w1_images:
            if np.all(image_bucket.data == 0):
                continue
            img_idx += 1
            self.plot_image(image_bucket, fig, rows, cols, img_idx, overlays, overlay_color, overlay_labels, overlay_label_color, image_contrast)

        r = img_idx % cols
        if r > 0:
            img_idx += cols - r

        # Collect W2 images
        images = []
        for w2_overlay in w2_overlays:
            try:
                coords_w2 = w2_overlay[0]
                epoch = w2_overlay[1]
                forward = w2_overlay[2]
                imageW2 = self.get_neowise_image(ra, dec, epoch, 2, img_size)
                if not imageW2:
                    break
                hdu = imageW2[0]
                header = hdu.header
                meanmjd = (header['MJDMIN']+header['MJDMAX'])/2
                images.append((hdu, self.get_epoch(meanmjd), coords_w2, forward))
            except Exception:
                print('A problem occurred while creating WISE time series for object ra={ra}, dec={dec}, epoch={epoch}, band=2'
                      .format(ra=ra, dec=dec, epoch=epoch))
                print(traceback.format_exc())

        # Process W2 images
        self.w2_images = []
        for image in images:
            coords_w2 = image[2]
            forward = image[3]
            w2, x, y, wcs = self.process_image_data(image[0], ra, dec, self.box_size)
            overlay_label = [coords[0] for coords in coords_w2]
            overlay_ra = [coords[1] for coords in coords_w2]
            overlay_dec = [coords[2] for coords in coords_w2]
            self.w2_images.append(self.ImageBucket(w2, x, y, 'W2', image[1], wcs, overlay_label, overlay_ra, overlay_dec, forward))

        # Plot W2 images
        for image_bucket in self.w2_images:
            if np.all(image_bucket.data == 0):
                continue
            img_idx += 1
            self.plot_image(image_bucket, fig, rows, cols, img_idx, overlays, overlay_color, overlay_labels, overlay_label_color, image_contrast)

        # Info text
        coords = SkyCoord(ra*u.deg, dec*u.deg)
        info_idx = math.ceil(img_idx / cols) * cols + 1
        fontsize = 2.6
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

    def create_light_curves(self, photometry_radius=5, yticks=None, open_file=None, file_format=None, overplot_l1b_phot=False, bin_l1b_phot=False,
                            save_l1b_phot=False, variability_threshold=0.1, legend_location='best', plot_statistics=False):
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
        overplot_l1b_phot : bool, optional
            Whether to overplot L1b photometry. The default is False.
        bin_l1b_phot : bool, optional
            Whether to bin L1b photometry by sky pass and plot the median magnitude. The default is False.
        save_l1b_phot : bool, optional
            Whether to save the L1b photometry table to the directory specified in the constructor ``unTimelyCatalogExplorer(directory=)``. The default is False.
        variability_threshold: float, optional
            The source is considered as variable if max_magnitude - mean_magnitude >= variability_threshold. The default is 0.1.
        legend_location : str, optional
            Matplotlib legend location string ('upper left', 'upper right', 'lower left', 'lower right', etc.). The default is 'best'.
        plot_statistics : bool, optional
            Whether to plot magnitude statistics below the light curves figure. The default is False.

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

        self.printout('Creating light curves ...')

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
            print(f'No unTimely photometry found in specified radius ({photometry_radius} arcsec) '
                  f'at given coordinates ({self.target_ra} {self.target_dec}).')
            return

        # Get W1 photometry
        mask = phot_table['band'] == 1
        phot_table_w1 = self.filter_table_by_min_distance_per_epoch(phot_table[mask])

        # Get W2 photometry
        mask = phot_table['band'] == 2
        phot_table_w2 = self.filter_table_by_min_distance_per_epoch(phot_table[mask])

        # Create plot
        plt.figure(figsize=(8, 4))
        plt.title(self.create_j_designation(ra, dec))

        yr1 = Time(phot_table_w1['mjdmean'], format='mjd').jyear
        yr2 = Time(phot_table_w2['mjdmean'], format='mjd').jyear
        w1 = phot_table_w1['mag'].value
        w2 = phot_table_w2['mag'].value
        e_w1 = phot_table_w1['dmag'].value
        e_w2 = phot_table_w2['dmag'].value

        plt.errorbar(yr1, w1, e_w1, lw=0.5, ms=2, marker='o', capsize=1.5, capthick=0.3, elinewidth=0.3, label='unTimely W1', zorder=2, c='blue')
        plt.errorbar(yr2, w2, e_w2, lw=0.5, ms=2, marker='o', capsize=1.5, capthick=0.3, elinewidth=0.3, label='unTimely W2', zorder=3, c='red')

        if plot_statistics:
            stats, peaks = self.create_light_curve_stats(yr1, w1, e_w1, 'unTimely W1', 'blue', variability_threshold)
            text = '\n unTimely W1 - Statistics: ' + stats + ' Notable peaks: ' + peaks
            plt.text(0, -0.1, text, ha='left', va='top', fontsize=5, transform=plt.gca().transAxes)

            stats, peaks = self.create_light_curve_stats(yr2, w2, e_w2, 'unTimely W2', 'red', variability_threshold)
            text = '\n unTimely W2 - Statistics: ' + stats + ' Notable peaks: ' + peaks
            plt.text(0, -0.15, text, ha='left', va='top', fontsize=5, transform=plt.gca().transAxes)

        if overplot_l1b_phot:
            phot_table = self.get_l1b_photometry(ra, dec, photometry_radius)

            if (len(phot_table) == 0):
                print(f'No L1b photometry found in specified radius ({photometry_radius} arcsec) '
                      f'at given coordinates ({self.target_ra} {self.target_dec}).')
            else:
                if save_l1b_phot:
                    result_file_name = 'l1b_phot_table_' + self.create_obj_name(ra, dec) + '.' + self.result_table_extension
                    phot_table.write(result_file_name, format=self.result_table_format, overwrite=True)

                sigma, maxiters = 3, None
                year = phot_table['year']

                if bin_l1b_phot:
                    yr = []
                    w1 = []
                    w2 = []
                    e_w1 = []
                    e_w2 = []

                    if len(phot_table) > 0:
                        year_bin = np.trunc(year / 0.5)
                        grouped = phot_table.group_by(year_bin)

                        for group in grouped.groups:
                            w1_clipped = sigma_clip(group['w1mpro'], sigma=sigma, maxiters=maxiters)
                            w2_clipped = sigma_clip(group['w2mpro'], sigma=sigma, maxiters=maxiters)
                            yr.append(np.ma.median(group['year']))
                            w1.append(np.ma.median(w1_clipped))
                            w2.append(np.ma.median(w2_clipped))
                            e_w1.append(self.std_error(w1_clipped))
                            e_w2.append(self.std_error(w2_clipped))

                    plt.errorbar(yr, w1, e_w1, lw=0.5, ms=2, marker='o', capsize=1.5, capthick=0.3, elinewidth=0.3, label='L1b median W1', zorder=0, c='lightskyblue')
                    plt.errorbar(yr, w2, e_w2, lw=0.5, ms=2, marker='o', capsize=1.5, capthick=0.3, elinewidth=0.3, label='L1b median W2', zorder=1, c='pink')

                    if plot_statistics:
                        stats, peaks = self.create_light_curve_stats(yr, w1, e_w1, 'L1b median W1', 'lightskyblue', variability_threshold)
                        text = '\n L1b median W1 - Statistics: ' + stats + ' Notable peaks: ' + peaks
                        plt.text(0, -0.2, text, ha='left', va='top', fontsize=5, transform=plt.gca().transAxes)

                        stats, peaks = self.create_light_curve_stats(yr, w2, e_w2, 'L1b median W2', 'pink', variability_threshold)
                        text = '\n L1b median W2 - Statistics: ' + stats + ' Notable peaks: ' + peaks
                        plt.text(0, -0.25, text, ha='left', va='top', fontsize=5, transform=plt.gca().transAxes)
                else:
                    w1_clipped = sigma_clip(phot_table['w1mpro'], sigma=sigma, maxiters=maxiters)
                    w2_clipped = sigma_clip(phot_table['w2mpro'], sigma=sigma, maxiters=maxiters)
                    plt.scatter(year, w1_clipped, s=3, label='L1b W1', zorder=0, c='lightskyblue')
                    plt.scatter(year, w2_clipped, s=3, label='L1b W2', zorder=1, c='pink')

        if yticks:
            plt.yticks(yticks)
        plt.xlabel('Year')
        plt.ylabel('Magnitude (mag)')
        plt.gca().invert_yaxis()
        plt.grid(color='grey', alpha=0.3, linestyle='-.', linewidth=0.2, axis='both')
        plt.tick_params(axis='both', length=0)
        lgnd = plt.legend(loc=legend_location)
        lgnd.get_frame().set_boxstyle('Square')

        # Save light curves plot
        filename = 'unTimely_Catalog_light_curves_' + self.create_obj_name(ra, dec) + '.' + file_format
        plt.savefig(filename, dpi=600, bbox_inches='tight', format=file_format)
        plt.close()

        # Open saved file
        if open_file:
            self.start_file(filename)

    def filter_table_by_min_distance_per_epoch(self, data_table):
        # Initialize an empty table to store the filtered results
        filtered_table = Table(names=data_table.colnames, dtype=data_table.dtype)

        # Iterate over unique epochs
        for epoch in np.unique(data_table['epoch']):
            # Select rows corresponding to the current epoch
            epoch_rows = data_table[data_table['epoch'] == epoch]

            # Find the row with the smallest distance
            min_distance_row = epoch_rows[np.argmin(epoch_rows['target_dist'])]

            # Append the row with the smallest distance to the filtered table
            filtered_table.add_row(min_distance_row)

        return filtered_table

    def create_light_curve_stats(self, time, magnitude, error, band, color, variability_threshold):
        """
        self.printout(band + ':')

        # Mean quartiles and differences
        q1, q2, q3 = self.calculate_mean_quartiles(magnitude)
        q3_q2 = '{:.3f}'.format(q3 - q2)
        q2_q1 = '{:.3f}'.format(q2 - q1)
        q3_q1 = '{:.3f}'.format(q3 - q1)
        q1 = '{:.3f}'.format(q1)
        q2 = '{:.3f}'.format(q2)
        q3 = '{:.3f}'.format(q3)
        self.printout(f'mean(Q1)={q1}, Q2={q2}, mean(Q3)={q3}, mean(Q3)-q2={q3_q2}, Q2-mean(Q1))={q2_q1}, mean(Q3)-mean(Q1)={q3_q1}')

        # Skewness and Kurtosis
        from scipy.stats import skew, kurtosis
        skewness = '{:.3f}'.format(skew(magnitude))
        kurtosis = '{:.3f}'.format(kurtosis(magnitude))
        self.printout("Skewness:", skewness)
        self.printout("Kurtosis:", kurtosis)

        self.printout()
        """

        # Calculate statistics
        min_magnitude = np.min(magnitude)
        max_magnitude = np.max(magnitude)
        mean_magnitude = np.mean(magnitude)
        # median_magnitude = np.median(magnitude)
        # std = np.std(magnitude)
        weighted_std = self.weighted_std(magnitude, error)
        var_index = self.variability_index(magnitude, error)
        # neumann_ratio = self.von_neumann_ratio(magnitude)
        # iqr = self.calculate_iqr(magnitude)

        stats = {
            'min': '{:.3f}'.format(min_magnitude),
            'max': '{:.3f}'.format(max_magnitude),
            'max-min': '{:.3f}'.format(max_magnitude - min_magnitude),
            'mean': '{:.3f}'.format(mean_magnitude),
            # 'median': '{:.3f}'.format(median_magnitude),
            # 'std': '{:.3f}'.format(std),
            'weighted std': '{:.3f}'.format(weighted_std),
            'variability index': '{:.3f}'.format(var_index),
            # 'von Neumann ratio': '{:.3f}'.format(neumann_ratio),
            # 'interquartile range': '{:.3f}'.format(iqr)
        }

        stats = ', '.join(f'{key}={value}' for key, value in stats.items())

        if max_magnitude - mean_magnitude < variability_threshold:
            self.printout(f'{band} magnitude fluctuations are below the specified/default threshold (variability_threshold={variability_threshold}).')
            return stats, 'none'

        # Find peaks in the magnitude data
        from scipy.signal import find_peaks
        data = pd.DataFrame({'time': time, 'magnitude': magnitude})
        peaks, _ = find_peaks(data['magnitude']*-1, height=-mean_magnitude+variability_threshold)

        # Get time and magnitude values at peak positions
        peak_times = data['time'].iloc[peaks]
        peak_magnitudes = data['magnitude'].iloc[peaks]

        plt.scatter(peak_times, peak_magnitudes, s=20, facecolors='none', edgecolors=color, linewidth=0.5, label=f'{band} Peaks')

        peak_dates = [Time(time, format='jyear').datetime.date().strftime('%Y-%m-%d') for time in peak_times]
        peak_magnitudes = peak_magnitudes.tolist()

        peaks = {}
        for date, magnitude in zip(peak_dates, peak_magnitudes):
            peaks[date] = '{:.3f}'.format(magnitude)

        peaks = ', '.join(f'{key}: {value}' for key, value in peaks.items())

        return stats, peaks

    def variability_index(self, magnitudes, errors):
        magnitudes = np.array(magnitudes)
        errors = np.array(errors)

        # Calculate the weighted mean magnitude
        weighted_mean = np.sum(magnitudes / errors ** 2) / np.sum(1.0 / errors ** 2)

        # Calculate the variability index
        var_index = np.sqrt(np.sum(((magnitudes - weighted_mean) / errors) ** 2) / (len(magnitudes) - 1))

        return var_index

    def weighted_std(self, values, weights):
        # Calculate the weighted mean
        weighted_mean = np.average(values, weights=weights)

        # Calculate the weighted sum of squares of differences
        weighted_sum_squares_diff = np.sum(weights * (values - weighted_mean) ** 2)

        # Calculate the sum of weights
        sum_weights = np.sum(weights)

        # Calculate the weighted standard deviation
        weighted_std = np.sqrt(weighted_sum_squares_diff / sum_weights)

        return weighted_std

    def von_neumann_ratio(self, values):
        # Calculate the differences between successive values
        differences = np.diff(values)

        # Calculate the mean square successive difference
        mean_square_successive_diff = np.mean(differences**2)

        # Calculate the distribution variance
        distribution_variance = np.var(values)

        # Calculate the von Neumann ratio
        von_neumann_ratio = mean_square_successive_diff / distribution_variance

        return von_neumann_ratio

    def calculate_iqr(self, values):
        # Calculate the first quartile (Q1)
        q1 = np.percentile(values, 25)

        # Calculate the third quartile (Q3)
        q3 = np.percentile(values, 75)

        # Calculate the Interquartile Range (IQR)
        iqr = q3 - q1

        return iqr

    def calculate_mean_quartiles(self, values):
        q1 = np.percentile(values, 25)
        q2 = np.percentile(values, 50)
        q3 = np.percentile(values, 75)
        mean_q1 = np.mean([x for x in values if x <= q1])
        mean_q3 = np.mean([x for x in values if x >= q3])
        return mean_q1, q2, mean_q3

    def create_image_blinks(self, overlays=False, blink_duration=300, image_zoom=10, image_contrast=None, separate_scan_dir=False, display_blinks=False):
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
        separate_scan_dir : bool, optional
            Whether to separate sky scans into forward and backward directions. The default is False.
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

        self.printout('Creating image blinks ...')

        ra = self.target_ra
        dec = self.target_dec
        w1_images = self.w1_images
        w2_images = self.w2_images

        # Separate scan directions
        if separate_scan_dir:
            # Separate W1 scans
            w1_forward = []
            w1_backward = []
            for w1_bucket in w1_images:
                if (w1_bucket.forward == 1):
                    w1_forward.append(w1_bucket)
                else:
                    w1_backward.append(w1_bucket)
            w1_images = w1_backward + w1_forward

            # Separate W2 scans
            w2_forward = []
            w2_backward = []
            for w2_bucket in w2_images:
                if (w2_bucket.forward == 1):
                    w2_forward.append(w2_bucket)
                else:
                    w2_backward.append(w2_bucket)
            w2_images = w2_backward + w2_forward

        # Draw settings
        stroke_width = 3
        circle_radius = 50
        point_radius = 2
        overlay_radius = 5
        red = (255, 0, 0)
        green = (0, 128, 0)

        try:
            font = ImageFont.truetype('arial.ttf', 12)
        except OSError:
            font = ImageFont.load_default()

        # Create animated GIF - W1 with overlays
        images = []

        for w1_bucket in w1_images:
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
            if overlays:
                overlay_ra = w1_bucket.overlay_ra
                overlay_dec = w1_bucket.overlay_dec
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
            draw.text((10, 10), 'W1 ' + year_obs, red, font=font)

            images.append(rgb_image)

        if images:
            filename = 'Animated_time_series_w1_' + self.create_obj_name(ra, dec) + '.gif'
            images[0].save(filename, save_all=True, append_images=images[1:], loop=0)
            if display_blinks:
                self.start_file(filename)

        # Create animated GIF - W2 with overlays
        images = []

        for w2_bucket in w2_images:
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
            if overlays:
                overlay_ra = w2_bucket.overlay_ra
                overlay_dec = w2_bucket.overlay_dec
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
            draw.text((10, 10), 'W2 ' + year_obs, red, font=font)

            images.append(rgb_image)

        if images:
            filename = 'Animated_time_series_w2_' + self.create_obj_name(ra, dec) + '.gif'
            images[0].save(filename, save_all=True, append_images=images[1:], loop=0)
            if display_blinks:
                self.start_file(filename)

        # Create animated GIF - W1+W2 without overlays
        images = []
        for i in range(min(len(w1_images), len(w2_images))):
            w1_bucket = w1_images[i]
            w2_bucket = w2_images[i]
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
            draw.text((10, 10), 'W1+W2 ' + year_obs, red, font=font)

            images.append(rgb_image)

        if images:
            filename = 'Animated_time_series_' + self.create_obj_name(ra, dec) + '.gif'
            images[0].save(filename, save_all=True, append_images=images[1:], loop=0)
            if display_blinks:
                self.start_file(filename)
