import os
import sys
import math
import certifi
import warnings
import tempfile
import traceback
import subprocess
import numpy as np
import pandas as pd
import hpgeom as hp
import pyarrow.fs
import pyarrow.compute
import pyarrow.dataset
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
import astropy.units as u
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp


class unTimelyCatalogExplorer:

    def __init__(self, directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, allow_insecure=False, suppress_console_output=False, ignore_warnings=True):
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

        Returns
        -------
        An unTimelyCatalogExplorer instance.

        """
        if ignore_warnings:
            warnings.simplefilter('ignore', category=Warning)
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

        return allwise, neowise

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

    def init_result_table(self):
        return Table(names=(
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

    def create_result_table(self, df, target_ra, target_dec, nearest_neighbor):
        table = Table.from_pandas(df)

        grouped = table.group_by(['band', 'EPOCH'])

        target_coords = SkyCoord([target_ra*u.deg], [target_dec*u.deg])

        self.w1_overlays = []
        self.w2_overlays = []

        tables = []

        group_number = 0

        for group in grouped.groups:
            result_table = self.init_result_table()

            catalog_coords = SkyCoord(group['ra'], group['dec'], unit='deg')
            target_dist = target_coords.separation(catalog_coords).arcsec
            group.add_column(target_dist, name='target_dist')
            group.sort('target_dist')

            coords_w1 = []
            coords_w2 = []

            group_number += 1
            object_number = 0

            for row in group:
                object_number += 1
                source_label = str(group_number) + '.' + str(object_number)

                """
                From Schlafly et al. 2019:
                    Fluxes and corresponding uncertainties are given in linear flux units, specifically, in Vega nanomaggies (nMgy; Finkbeiner et al. 2004).
                    The corresponding Vega magnitudes are given by 22.5-2.5log10(flux).
                    The agreement between unWISE and AllWISE magnitudes can be improved by subtracting 4 mmag and 32 mmag from W1 and W2.
                """
                # Calculate Vega magnitude from flux
                # flux_corr = 4 if row['band'] == 1 else 32
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
                    row['target_dist'],
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

                if row['band'] == 1:
                    coords_w1.append((source_label, row['ra'], row['dec']))
                else:
                    coords_w2.append((source_label, row['ra'], row['dec']))

                if nearest_neighbor:
                    break

            tables.append(result_table)

            if group[0]['band'] == 1:
                self.w1_overlays.append((coords_w1, group[0]['EPOCH'],  group[0]['FORWARD']))
            else:
                self.w2_overlays.append((coords_w2, group[0]['EPOCH'],  group[0]['FORWARD']))

        return vstack(tables)

    def search_by_coordinates(self, target_ra, target_dec, box_size=100, cone_radius=None, nearest_neighbor=False, show_result_table_in_browser=False,
                              save_result_table=True, result_table_format='ascii', result_table_extension='dat'):
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

        Returns
        -------
        result_table : astropy.table.table.Table
            Result table containing the catalog entries located within a field of view of the specified size at the given coordinates.

        """
        self.printout('Gathering data (may take several minutes) ...')

        self.target_ra = target_ra
        self.target_dec = target_dec
        self.box_size = box_size
        self.img_size = round(box_size / self.pixel_scale)

        size = box_size * u.arcsec.to(u.deg)

        dec_offset = size / 2
        dec_min, dec_max = target_dec - dec_offset, target_dec + dec_offset  # deg
        dec_mean = (dec_min + dec_max) / 2
        ra_offset = math.degrees(math.radians(size) / math.cos(math.radians(dec_mean))) / 2
        ra_min, ra_max = target_ra - ra_offset, target_ra + ra_offset  # deg
        polygon_corners = [(ra_min, dec_min), (ra_min, dec_max), (ra_max, dec_max), (ra_max, dec_min)]

        K = 5  # HEALPix order at which the dataset is partitioned
        nside = hp.order_to_nside(K)
        polygon_pixels = hp.query_polygon(
            nside=nside,
            a=[corner[0] for corner in polygon_corners],
            b=[corner[1] for corner in polygon_corners],
            nest=True,
            inclusive=True
        )

        fs = pyarrow.fs.S3FileSystem(region='us-west-2')
        bucket = 'nasa-irsa-wise'
        catalog_root = f'{bucket}/unwise/neo7/catalogs/time_domain/healpix_k{K}/unwise-neo7-time_domain-healpix_k{K}.parquet'
        parquet_ds = pyarrow.dataset.parquet_dataset(f'{catalog_root}/_metadata', filesystem=fs, partitioning='hive')

        region_tbl = parquet_ds.to_table(
            filter=(pyarrow.compute.field(f'healpix_k{K}').isin(polygon_pixels)
                    & (pyarrow.compute.field('ra') > ra_min)
                    & (pyarrow.compute.field('ra') < ra_max)
                    & (pyarrow.compute.field('dec') > dec_min)
                    & (pyarrow.compute.field('dec') < dec_max))
        )

        if cone_radius:
            radius = cone_radius * u.arcsec.to(u.deg)
            target_coord = SkyCoord(ra=target_ra * u.degree, dec=target_dec * u.degree)
            region_coords = SkyCoord(ra=region_tbl['ra'] * u.degree, dec=region_tbl['dec'] * u.degree)
            region_tbl = region_tbl.filter(target_coord.separation(region_coords).degree < radius)

        if len(region_tbl) == 0:
            self.printout(f'No sources found for given coordinates={target_ra} {target_dec}, box size={box_size}, cone radius={cone_radius}')
            self.result_table = self.init_result_table()
        else:
            df = region_tbl.to_pandas()
            self.result_table = self.create_result_table(df, target_ra, target_dec, nearest_neighbor)

        if save_result_table:
            result_file_name = 'unTimely_Catalog_search results_' + self.create_obj_name(target_ra, target_dec) + '.' + result_table_extension
            self.result_table.write(result_file_name, format=result_table_format, overwrite=True)

        if show_result_table_in_browser:
            self.result_table.show_in_browser(jsviewer=True)

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
        if self.result_table is None:
            raise Exception('Method ``search_by_coordinates`` must be called first.')

        if len(self.result_table) == 0:
            return

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
                            variability_threshold=0.1, legend_location='best', plot_statistics=False):
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

        if len(self.result_table) == 0:
            return

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
            print(f'No photometry found in specified radius ({photometry_radius} arcsec) '
                  f'at given coordinates ({self.target_ra} {self.target_dec}) to create light curves.')
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
            allwise, neowise = self.get_l1b_photometry(ra, dec, photometry_radius)
            allwise['mjd'].unit = 'd'
            neowise['mjd'].unit = 'd'

            allwise_year = Time(allwise['mjd'], format='mjd').jyear
            neowise_year = Time(neowise['mjd'], format='mjd').jyear

            sigma = 3
            maxiters = None

            if bin_l1b_phot:
                yr = []
                w1 = []
                w2 = []
                e_w1 = []
                e_w2 = []

                if len(allwise) > 0:
                    allwise.add_column(allwise_year, name='year')
                    year_bin = np.trunc(allwise_year / 0.5)
                    grouped = allwise.group_by(year_bin)

                    for group in grouped.groups:
                        w1_clipped = sigma_clip(group['w1mpro_ep'], sigma=sigma, maxiters=maxiters)
                        w2_clipped = sigma_clip(group['w2mpro_ep'], sigma=sigma, maxiters=maxiters)
                        yr.append(np.ma.median(group['year']))
                        w1.append(np.ma.median(w1_clipped))
                        w2.append(np.ma.median(w2_clipped))
                        e_w1.append(self.std_error(w1_clipped))
                        e_w2.append(self.std_error(w2_clipped))

                if len(neowise) > 0:
                    neowise.add_column(neowise_year, name='year')
                    year_bin = np.trunc(neowise_year / 0.5)
                    grouped = neowise.group_by(year_bin)

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
                w1_clipped = sigma_clip(allwise['w1mpro_ep'], sigma=sigma, maxiters=maxiters)
                w2_clipped = sigma_clip(allwise['w2mpro_ep'], sigma=sigma, maxiters=maxiters)
                plt.scatter(allwise_year, w1_clipped, s=3, zorder=0, c='lightskyblue')
                plt.scatter(allwise_year, w2_clipped, s=3, zorder=1, c='pink')

                w1_clipped = sigma_clip(neowise['w1mpro'], sigma=sigma, maxiters=maxiters)
                w2_clipped = sigma_clip(neowise['w2mpro'], sigma=sigma, maxiters=maxiters)
                plt.scatter(neowise_year, w1_clipped, s=3, label='L1b W1', zorder=0, c='lightskyblue')
                plt.scatter(neowise_year, w2_clipped, s=3, label='L1b W2', zorder=1, c='pink')

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
        if self.result_table is None:
            raise Exception('Method ``create_finder_charts`` must be called first.')

        if len(self.result_table) == 0:
            return

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
