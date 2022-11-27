# unTimely_Catalog_explorer

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

A search and visualization tool for the unTimely Catalog[[1]](#1), a full-sky, time-domain unWISE catalog.

The tool allows to:
- search the catalog by coordinates (box or cone search) ```search_by_coordinates(ra, dec)```,
- create finder charts for W1 and W2 at each epoch with overplotted catalog positions (overlays) ```create_finder_charts()```,
- create light curves using W1 and W2 photometry of all available epochs ```create_light_curves()```,
- create W1 and W2 image blinks with overplotted catalog positions in GIF format ```create_image_blinks()```.

## Module dependencies:

The Python Standard Library, NumPy, Matplotlib, Astropy and Pillow (PIL Fork)

## Usage example

Most of the parameters can be omitted as they have default values (see [API documentation](#apidoc) for more details).

An instance of the unTimelyCatalogExplorer class has to be created first ```ucx = unTimelyCatalogExplorer()```.

The ```search_by_coordinates``` method must always be called before all others.

The ```create_image_blinks``` method depends on the results of the ```create_finder_charts``` method.

```python
from unTimely_Catalog_tools import unTimelyCatalogExplorer
import tempfile

ucx = unTimelyCatalogExplorer(directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
                              catalog_base_url='http://unwise.me/data/neo7/untimely-catalog/',
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

## Example output

### Result table:
```
|object_label|target_dist|        x|        y|      flux|          dx|          dy|     dflux|        qf|     rchi2|  fracflux|   fluxlbs|  dfluxlbs|     fwhm|  spread_model|dspread_model|   fluxiso|          xiso|          yiso|         sky|       ra|      dec|coadd_id|band|          unwise_detid|  nm|primary|flags_unwise|flags_info|epoch|forward|   mjdmin|   mjdmax|  mjdmean|       mag|       dmag|
|        char|      float|    float|    float|     float|       float|       float|     float|     float|     float|     float|     float|     float|    float|         float|        float|     float|         float|         float|       float|    float|    float|    char|long|                  char|long|   long|        long|      long| long|   long|    float|    float|    float|     float|      float|
|            |     arcsec|      pix|      pix|      nMgy|         pix|         pix|      nMgy|          |          |          |      nMgy|      nMgy|      pix|              |             |          |              |              |        nMgy|      deg|      deg|        |    |                      |    |       |            |          |     |       |        d|        d|        d|       mag|        mag|
|        null|       null|     null|     null|      null|        null|        null|      null|      null|      null|      null|      null|      null|     null|          null|         null|      null|          null|          null|        null|     null|     null|    null|null|                  null|null|   null|        null|      null| null|   null|     null|     null|     null|      null|       null|
          1.1   41.133976 267.98914 366.07266  1478.2817   0.02711302  0.026857046  23.787693  0.9999999    1.45462  0.9730243  1476.4236  24.258665 2.4941554  -7.581711e-05  0.0006082071  1478.2878  -0.0029756483  -0.0013511414  0.045344237 26.970234  23.65304 0264p242    1 0264p242w1o0003274e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  14.578548 0.019259011 
          1.2   24.311607 270.75854 354.26965  186.80829   0.21058543   0.20848638  22.858791  0.9999999 0.40986335 0.85064834  182.78297  23.324059   2.49984    0.002097249   0.004466529  186.80855   0.0067806835    0.018961437  -0.20955364 26.980085  23.65512 0264p242    1 0264p242w1o0003307e000   11       1            0          0     0       0 55213.777 55216.355 55215.066   16.84501 0.019311216 
          1.3   26.611568 272.67844 349.33646   350.6362  0.112922154  0.111713514  23.004868        1.0  0.7470472  0.9601253  341.51666   23.47193  2.502432  0.00046521425  0.0023680562  350.63937   -0.007589876   0.0136228325  -0.45855537 26.984205 23.656572 0264p242    1 0264p242w1o0003336e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  16.150314 0.019302838 
          1.4   14.474302  277.3332 361.20578   642.6594  0.061838116  0.061267607  23.219679 0.99999994  0.5048074  0.8716953   635.2518  23.687193 2.5000417   0.0020625591  0.0013212693   642.6394    0.013772808    0.041024134  -0.27924097 26.974323 23.660162 0264p242    1 0264p242w1o0003390e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  15.486827 0.019290635 
          1.5   26.739912 279.12183  346.6229  195.84848   0.20177622   0.19973092  23.160677 0.99999994 0.37018228  0.9257498  199.33455   23.62958 2.5000474   0.0062518716   0.004382402  195.82774   0.0010374679   -0.012010553   -1.0438267  26.98649 23.661486 0264p242    1 0264p242w1o0003413e000   10       1            0          0     0       0 55213.777 55216.355 55215.066  16.792604 0.019293973 
          1.6  0.87831515  279.6703 356.38605  2937.9165  0.014060003 0.0139493905  24.824362        1.0  3.6525035 0.97762567  2932.1765   25.30369 2.5033457   0.0012370348 0.00031610663  2938.0347    -0.00686241    0.017446637   -0.6429966  26.97835 23.661934 0264p242    1 0264p242w1o0003427e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  13.831381 0.019203572 
          1.7   31.302246 285.21735 346.58777  349.63837  0.113236874   0.11249279  23.257622        1.0  0.3575519  0.9440686  347.69266  23.729797 2.5018487    0.000954926  0.0024212874   349.5952    0.008859864   -0.010896928   -1.5677048  26.98654 23.666142 0264p242    1 0264p242w1o0003492e000   10       1            0          0     0       0 55213.777 55216.355 55215.066  16.153444 0.019288493 
          1.8   56.389893 291.25018 373.05084   165.4573   0.23699896   0.23463635   22.81245        1.0 0.71865296 0.99898523  169.73457  23.279415 2.4965863   -0.008771837     0.0052296  165.48672     0.19331563   -0.053440556  -0.75371563  26.96449 23.670828 0264p242    1 0264p242w1o0003565e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  16.979856 0.019313887 
          1.9   48.941624 295.82776 363.08554  132.05058    0.2975982   0.29347035   22.78306        1.0  0.6937911 0.98462033  142.57121  23.250383 2.4956353  -0.0073519945   0.006480334  132.13405    -0.08966763    -0.56501573   -1.0877496 26.972816 23.674294 0264p242    1 0264p242w1o0003624e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  17.231546 0.019315584 
         1.10   55.028603 263.89313 343.63528  106.68029   0.36999014    0.3683969  22.931067        1.0 0.59704506 0.91822433  99.514885  23.393208 2.4997027  -0.0054461956   0.007932083  106.65743      0.2900153    -0.20782283  -0.23481975  26.98893 23.649845 0264p242    1 0264p242w1o0026338e000   11       1            0          0     0       0 55213.777 55216.355 55215.066  17.471283 0.019307062 
```
![Full result table](Example%20output/unTimely_Catalog_search%20results_26.978383%2B23.661691.dat)

### Finder charts:
![Finder charts](Example%20output/unTimely_Catalog_finder_charts_26.978383%2B23.661691.png)

### Light curves:
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691.png)
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691_median_L1b.png)

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
    Directory where the finder charts should be saved. The default is tempfile.gettempdir().
- cache : bool, optional  
    Whether to cache the downloaded files. The default is True.
- show_progress : bool, optional  
    Whether to show the file download progress. The default is True.
- timeout : int, optional  
    Timeout for remote requests in seconds. The default is 300.
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
    Whether to plot catalog entry labels on the finder charts (tied to overlay circles). The default is False.
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
create_light_curves(photometry_radius=5, yticks=None, open_file=None, file_format=None, overplot_l1b_phot=False, bin_l1b_phot=False):
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
