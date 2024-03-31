# unTimely_Catalog_explorer

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

<a href="https://ascl.net/2211.005"><img src="https://img.shields.io/badge/ascl-2211.005-blue.svg?colorB=262255" alt="ascl:2211.005" /></a>

A search and visualization tool for the unTimely Catalog[[1]](#1), a full-sky, time-domain unWISE catalog.

The tool allows to:
- search the catalog by coordinates (box or cone search) ```search_by_coordinates(ra, dec)```,
- create finder charts for W1 and W2 at each epoch with overplotted catalog positions (overlays) ```create_finder_charts()```,
- create light curves using W1 and W2 photometry of all available epochs ```create_light_curves()```,
- create W1 and W2 image blinks with overplotted catalog positions in GIF format ```create_image_blinks()```.

## Module dependencies

The Python Standard Library, NumPy, Matplotlib, Pillow (PIL Fork), Certifi, Scipy, Astropy and Reproject (which is an Astropy affiliated package)

## Installation

The code can be installed as follows:
```
git clone https://github.com/fkiwy/unTimely_Catalog_explorer.git
cd unTimely_Catalog_explorer
python setup.py install
```

## Example usage

Most of the parameters can be omitted as they have default values (see [API documentation](#apidoc) for more details).

An instance of the unTimelyCatalogExplorer class has to be created first ```ucx = unTimelyCatalogExplorer()```.

The ```search_by_coordinates``` method must always be called before all others.

The ```create_image_blinks``` method depends on the results of the ```create_finder_charts``` method.

```python
from unTimely_Catalog_tools import unTimelyCatalogExplorer
import tempfile

ucx = unTimelyCatalogExplorer(directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300, suppress_console_output=False,
                              catalog_base_url='https://unwise.me/data/neo7/untimely-catalog/',
                              catalog_index_file='untimely_index-neo7.fits')

result_table = ucx.search_by_coordinates(26.9783833, 23.6616914, box_size=100, cone_radius=None, show_result_table_in_browser=False,
                                         save_result_table=True, result_table_format='ascii.ipac', result_table_extension='dat')

# Do whatever you want with the result table data here
print(result_table.info)

ucx.create_finder_charts(overlays=True, overlay_color='green', overlay_labels=False, overlay_label_color='red',
                         image_contrast=5, open_file=False, file_format='pdf')

ucx.create_light_curves(photometry_radius=2, yticks=None, open_file=False, file_format='png', overplot_l1b_phot=True, bin_l1b_phot=True)

ucx.create_image_blinks(blink_duration=300, image_zoom=10, image_contrast=5, separate_scan_dir=False, display_blinks=False)
```

Note that the default output directory for the finder charts, image blinks and light curve plots is your system temp folder (given by ``tempfile.gettempdir()``).

## Example output

### Result table:
```
|source_label|target_dist|        x|        y|      flux|          dx|          dy|     dflux|        qf|     rchi2|  fracflux|   fluxlbs|  dfluxlbs|     fwhm|  spread_model|dspread_model|   fluxiso|          xiso|          yiso|         sky|       ra|      dec|coadd_id|band|          unwise_detid|  nm|primary|flags_unwise|flags_info|epoch|forward|   mjdmin|   mjdmax|  mjdmean|       mag|        dmag|flags_unwise_bits|flags_unwise_descr|flags_info_bits|flags_info_descr|
|        char|      float|    float|    float|     float|       float|       float|     float|     float|     float|     float|     float|     float|    float|         float|        float|     float|         float|         float|       float|    float|    float|    char|long|                  char|long|   long|        long|      long| long|   long|    float|    float|    float|     float|       float|             char|              char|           char|            char|
|            |     arcsec|      pix|      pix|      nMgy|         pix|         pix|      nMgy|          |          |          |      nMgy|      nMgy|      pix|              |             |          |              |              |        nMgy|      deg|      deg|        |    |                      |    |       |            |          |     |       |        d|        d|        d|       mag|         mag|                 |                  |               |                |
|        null|       null|     null|     null|      null|        null|        null|      null|      null|      null|      null|      null|      null|     null|          null|         null|      null|          null|          null|        null|     null|     null|    null|null|                  null|null|   null|        null|      null| null|   null|     null|     null|     null|      null|        null|             null|              null|           null|            null|
          1.1  0.87831515  279.6703 356.38605  2937.9165  0.014060003 0.0139493905  24.824362        1.0  3.6525035 0.97762567  2932.1765   25.30369 2.5033457   0.0012370348 0.00031610663  2938.0347    -0.00686241    0.017446637   -0.6429966  26.97835 23.661934 0264p242    1 0264p242w1o0003427e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  13.829902  0.009174345                                                                       
          1.2    6.800647  278.3556 358.60797  295.64267    0.1356889   0.13655952   23.57525        1.0  1.1573884 0.36883256   290.9348  24.039568 2.5038018    0.004512131  0.0028988284  295.75702     0.03663178    -0.11140295  -0.40257493 26.976494 23.660936 0264p242    1 0264p242w1o0030758e000   10       1            0          0     0       0 55213.777 55216.355 55215.066  16.323082   0.08676343                                                                       
          1.3   10.734045 276.22272 354.01282  139.53488    0.2820289    0.2787143  22.819721        1.0  0.3663603  0.5460553  136.47597  23.282099 2.4959884  -0.0023428798   0.006013545  139.55301    0.093110986   -0.051865812     -0.50946 26.980318 23.659294 0264p242    1 0264p242w1o0030751e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  17.138292   0.17917164                                                                       
          1.4   12.830857   282.567 359.72888  250.25577   0.16098073   0.15887521  23.201101        1.0  1.8302455  0.7132105   241.6268  23.660118  2.502495    0.015723884  0.0036360174  250.27533   -0.039657693    -0.03137485   -0.5397417 26.975573 23.664156 0264p242    1 0264p242w1o0026389e000   11       1            0          0     0       0 55213.777 55216.355 55215.066   16.50404   0.10094801                                                                       
          1.5   14.474302  277.3332 361.20578   642.6594  0.061838116  0.061267607  23.219679 0.99999994  0.5048074  0.8716953   635.2518  23.687193 2.5000417   0.0020625591  0.0013212693   642.6394    0.013772808    0.041024134  -0.27924097 26.974323 23.660162 0264p242    1 0264p242w1o0003390e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  15.480048  0.039245375                                                                       
          1.6   16.850346 273.78732  358.9083   94.68198   0.41980743   0.41060808  22.843729        1.0 0.46448976  0.5553647  86.752884   23.30474 2.4963253   0.0048452616   0.008926477   94.59819    -0.16292378    -0.05268831  -0.24089655 26.976227 23.657448 0264p242    1 0264p242w1o0030746e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  17.559332    0.2672214                                                                       
          1.7   24.311607 270.75854 354.26965  186.80829   0.21058543   0.20848638  22.858791  0.9999999 0.40986335 0.85064834  182.78297  23.324059   2.49984    0.002097249   0.004466529  186.80855   0.0067806835    0.018961437  -0.20955364 26.980085  23.65512 0264p242    1 0264p242w1o0003307e000   11       1            0          0     0       0 55213.777 55216.355 55215.066   16.82151   0.13352522                                                                       
          1.8   25.183413 270.22202 357.05176 117.079636    0.3422278   0.33152387   23.25296  0.9999999  0.5889581 0.74251676 112.012596  23.722586 2.4961278   -0.009067297  0.0075610518  117.14969    -0.08486951    -0.15145485  -0.10735986 26.977764 23.654718 0264p242    1 0264p242w1o0026349e000   10       1            0          0     0       0 55213.777 55216.355 55215.066  17.328796   0.21854028                                                                       
          1.9   26.611568 272.67844 349.33646   350.6362  0.112922154  0.111713514  23.004868        1.0  0.7470472  0.9601253  341.51666   23.47193  2.502432  0.00046521425  0.0023680562  350.63937   -0.007589876   0.0136228325  -0.45855537 26.984205 23.656572 0264p242    1 0264p242w1o0003336e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  16.137857  0.071336426                                                                       
         1.10   26.739912 279.12183  346.6229  195.84848   0.20177622   0.19973092  23.160677 0.99999994 0.37018228  0.9257498  199.33455   23.62958 2.5000474   0.0062518716   0.004382402  195.82774   0.0010374679   -0.012010553   -1.0438267  26.98649 23.661486 0264p242    1 0264p242w1o0003413e000   10       1            0          0     0       0 55213.777 55216.355 55215.066  16.770199   0.12900075                                                                       
```
![Full result table](Example%20output/unTimely_Catalog_search%20results_26.978383%2B23.661691.dat)

**Column description:**
```
       Name          Type   Unit                                                                       Description                                                                     
------------------ ------- ------ -----------------------------------------------------------------------------------------------------------------------------------------------------
source_label       S32     None   Unique source label within a specific result set that can be used to retrieve the corresponding source on the finder charts                          
target_dist        float32 arcsec Angular distance to the target coordinates                                                                                                           
x                  float32 pix    x coordinate                                                                                                                                         
y                  float32 pix    y coordinate                                                                                                                                         
flux               float32 nMgy   Vega flux                                                                                                                                            
dx                 float32 pix    x uncertainty                                                                                                                                        
dy                 float32 pix    y uncertainty                                                                                                                                        
dflux              float32 nMgy   formal flux uncertainty                                                                                                                              
qf                 float32 None   PSF-weighted fraction of good pixels                                                                                                                 
rchi2              float32 None   PSF-weighted average chi2                                                                                                                            
fracflux           float32 None   PSF-weighted fraction of flux from this source                                                                                                       
fluxlbs            float32 nMgy   FWHM of PSF at source location                                                                                                                       
dfluxlbs           float32 nMgy   local-background-subtracted flux                                                                                                                     
fwhm               float32 pix    formal fluxlbs uncertainty                                                                                                                           
spread_model       float32 None   SExtractor-like source size parameter                                                                                                                
dspread_model      float32 None   uncertainty in spread_model                                                                                                                          
fluxiso            float32 None   flux derived from linear least squares fit to neighbor-subtracted image; significant difference from ordinary flux indicates a convergence issue     
xiso               float32 None   x coordinate derived from linear least squares fit to neighbor-subtracted image; significant difference from ordinary x indicates a convergence issue
yiso               float32 None   y coordinate derived from linear least squares fit to neighbor-subtracted image; significant difference from ordinary y indicates a convergence issue
sky                float32 nMgy   residual sky at source location                                                                                                                      
ra                 float32 deg    R.A.                                                                                                                                                 
dec                float32 deg    decl.                                                                                                                                                
coadd_id           S32     None   unWISE/AllWISE coadd_id of source                                                                                                                    
band               int32   None   1 for W1, 2 for W2                                                                                                                                   
unwise_detid       S32     None   detection ID, unique in catalog                                                                                                                      
nm                 int32   None   number of images in coadd at source                                                                                                                  
primary            int32   None   source located in primary region of coadd                                                                                                            
flags_unwise       int32   None   unWISE flags at source location                                                                                                                      
flags_info         int32   None   additional flags at source location                                                                                                                  
epoch              int32   None   unWISE epoch number                                                                                                                                  
forward            int32   None   boolean, were input frames acquired pointing forward (1) or backward (0) along Earth's orbit                                                         
mjdmin             float32 d      MJD value of earliest contributing exposure                                                                                                          
mjdmax             float32 d      MJD value of latest contributing exposure                                                                                                            
mjdmean            float32 d      mean of MJDMIN and MJDMAX                                                                                                                            
mag                float32 mag    Vega magnitude given by 22.5-2.5log10(flux)                                                                                                            
dmag               float32 mag    magnitude uncertainty                                                                                                                                
flags_unwise_bits  S32     None   unWISE flags bits                                                                                                                                    
flags_unwise_descr S32     None   unWISE flags description                                                                                                                             
flags_info_bits    S32     None   info flags bits                                                                                                                                      
flags_info_descr   S32     None   info flags description                                                                                                                               
```

### Finder charts:
![Finder charts](Example%20output/unTimely_Catalog_finder_charts_26.978383%2B23.661691.png)

### Light curves:
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691.png)
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691_L1b.png)
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691_median_L1b.png)
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691_stats.png)

### Image blinks:
![Image blinks - variable](Example%20output/Animated_time_series_w1_26.978383%2B23.661691.gif) | ![Image blinks - color](Example%20output/Animated_time_series_26.978383%2B23.661691.gif)
![Image blinks - high PM](Example%20output/Animated_time_series_w2_133.79476-7.245146.gif)

## <a id="apidoc">API documentation</a>

## <kbd>class</kbd> `unTimelyCatalogExplorer`

### <kbd>constructor</kbd> `__init__`
```python
__init__(directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
         catalog_base_url='https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/',
         catalog_index_file='untimely_index-neo7.fits'):
```
Creates an unTimelyCatalogExplorer instance with the given parameters

#### <ins>Parameters</ins>
- directory : str, optional  
    Directory where the finder charts, image blinks and light curve plots should be saved. The default is tempfile.gettempdir().
- cache : bool, optional  
    Whether to cache the downloaded files. The default is True.
- show_progress : bool, optional  
    Whether to show the file download progress. The default is True.
- timeout : int, optional  
    Timeout for remote requests in seconds. The default is 300.
- allow_insecure : bool, optional  
    Whether to allow downloading files over a TLS/SSL connection even when the server certificate verification failed. The default is False.
- suppress_console_output : bool, optional  
    Whether to suppress all console output except error messages. The default is False.
- catalog_base_url : str, optional  
    Base URL to access the unTimely Catalog. The default is 'https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/'.
- catalog_index_file : str, optional  
    Catalog index file name. The default is 'untimely_index-neo7.fits'.

#### <ins>Returns</ins>
An unTimelyCatalogExplorer instance.

---

### <kbd>method</kbd> `search_by_coordinates`
```python
search_by_coordinates(target_ra, target_dec, box_size=100, cone_radius=None, show_result_table_in_browser=False,
                      save_result_table=True, result_table_format='ascii', result_table_extension='dat'):
```
Search the catalog by coordinates (box search).

#### <ins>Parameters</ins>
- ra : float  
    Right ascension in decimal degrees.
- dec : float  
    Declination in decimal degrees.
- box_size : int, optional  
    Box search size and/or image size in arcseconds. The default is 100.
- cone_radius : int, optional  
    Cone search radius in arcseconds. If specified, a cone search will be performed (instead of a box search) around the given coordinates within the given radius.
    However, the value of the ``box_size`` parameter still defines the image size of the finder charts and image blinks.
- show_result_table_in_browser : bool, optional  
    Whether to show the result table in your browser (columns can be sorted). The default is False.
- save_result_table : bool, optional  
    Whether to save the result table to the directory specified in the constructor ``unTimelyCatalogExplorer(directory=)``. The default is True.
- result_table_format : str, optional  
    Result table output format. The default is 'ascii'.
- result_table_extension : str, optional  
    Result table file extension. The default is 'dat'.
-  multi_processing : bool, optional
     Whether to allow multi-processing when downloading and scanning unTimely catalog files (faster but higher CPU usage). The default is False.
     ``unTimelyCatalogExplorer`` methods must be called within following ``if`` statement: ``if __name__ == '__main__':``

#### <ins>Returns</ins>
Astropy table containing the catalog entries located within a field of view of the specified size at the given coordinates.

---

### <kbd>method</kbd> `create_finder_charts`
```python
create_finder_charts(overlays=True, overlay_color='green', overlay_labels=False, overlay_label_color='red',
                     image_contrast=3, open_file=False, file_format='pdf'):
```
Create finder charts for W1 and W2 at each epoch with overplotted catalog positions (overlays)

#### <ins>Parameters</ins>
- overlays : bool, optional  
    Whether to plot W1 and W2 catalog positions on the finder charts (overlays). The default is True.
- overlay_color : str, optional  
    Overlay color. The default is 'green'.
- overlay_labels : bool, optional  
    Whether to plot overlay labels (corresponding to source labels in the result table) on the finder charts. The default is False.
- overlay_label_color : str, optional  
    Label color. The default is 'red'.
- image_contrast : int, optional  
    Contrast of W1 and W2 images. The default is 3.
- open_file : bool, optional  
    Whether to open the saved finder charts automatically. The default is False.
- file_format : str, optional  
    Output file format: pdf, png, eps, etc.. The default is 'pdf'.

#### <ins>Raises</ins>
Exception if method ``search_by_coordinates`` has not been called first.

#### <ins>Returns</ins>
None.

---

### <kbd>method</kbd> `create_light_curves`
```python
create_light_curves(photometry_radius=5, yticks=None, open_file=None, file_format=None, overplot_l1b_phot=False, bin_l1b_phot=False,
                    legend_location='best'):
```
Create light curves using W1 and W2 photometry of all available epochs.

#### <ins>Parameters</ins>
- photometry_radius : float, optional  
    Radius to search for the photometry used to create the light curves. The default is 5.
- yticks : tuple, optional  
    Tuple containing y axis tick values. The default is pyplot's automatic tick allocation.
- open_file : bool, optional  
    Whether to open the saved light curves automatically. The default is None (value given by method ``create_finder_charts`` will be used).
- file_format : str, optional  
    Output file format: pdf, png, eps, etc.. The default is None (value given by method ``create_finder_charts`` will be used).
- overplot_l1b_phot : bool, optional  
    Whether to overplot L1b photometry. The default is False.
- bin_l1b_phot : bool, optional  
    Whether to bin L1b photometry by sky pass and plot the median magnitude. The default is False.
- variability_threshold: float, optional  
    The source is considered as variable if max_magnitude - mean_magnitude >= variability_threshold. The default is 0.1.
- legend_location : str, optional  
    Matplotlib legend location string ('upper left', 'upper right', 'lower left', 'lower right', etc.). The default is 'best'.
- plot_statistics : bool, optional  
    Whether to plot magnitude statistics below the light curves figure. The default is False.

#### <ins>Raises</ins>
Exception if method ``search_by_coordinates`` has not been called first.

#### <ins>Returns</ins>
None.

---

### <kbd>method</kbd> `create_image_blinks`
```python
create_image_blinks(blink_duration=300, image_zoom=10, image_contrast=None, separate_scan_dir=False, display_blinks=False):
```
Create W1 and W2 image blinks with overplotted catalog positions in GIF format.

#### <ins>Parameters</ins>
- blink_duration : int, optional  
    Duration each image is shown in milliseconds. The default is 200.
- image_zoom : int, optional  
    Scaling factor to be applied on W1 and W2 images. The default is 10.
- image_contrast : int, optional  
    Contrast of W1 and W2 images. The default is None (value given by method ``create_finder_charts`` will be used).
- separate_scan_dir : bool, optional  
    Whether to separate sky scans into forward and backward directions. The default is False.
- display_blinks : bool, optional  
    Whether to display the image blinks in your system's media player. The default is False.

#### <ins>Raises</ins>
Exception if method ``create_finder_charts`` has not been called first.

#### <ins>Returns</ins>
None.

---

## References
<a id="1">[1]</a> Meisner, A. M., Caselden, D., and Schlafly, E. F., "unTimely: a Full-sky, Time-Domain unWISE Catalog", 2022. [![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/2209.14327)
