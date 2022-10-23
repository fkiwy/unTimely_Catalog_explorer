# unTimely_Catalog_explorer

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

A search and visualization tool for the unTimely Catalog[[1]](#1), a full-sky, time-domain unWISE catalog.

The tool allows to:
- search the catalog by coordinates (box search) ```search_by_coordinates(ra, dec)```,
- create finder charts for W1 and W2 at each epoch with overplotted catalog positions (overlays) ```create_finder_charts()```,
- create light curves using W1 and W2 photometry of all available epochs ```create_ligh_curves()```,
- create W1 and W2 image blinks with overplotted catalog positions in GIF format ```create_image_blinks()```.

### Module dependencies:
The Python Standard Library, NumPy, Matplotlib, Astropy and Pillow (PIL Fork)

### Usage example

Most of the parameters can be omitted as they have default values (see [docstrings](#docstrings) for more details).

An instance of the unTimelyCatalogExplorer class has to be created first ```ucx = unTimelyCatalogExplorer()```.

The ```search_by_coordinates``` method must always be called before all others.

The ```create_image_blinks``` method depends on the results of the ```create_finder_charts``` method.

```python
from unTimely_Catalog_tools import unTimelyCatalogExplorer
import tempfile

ucx = unTimelyCatalogExplorer(directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
                              catalog_base_url='http://unwise.me/data/neo7/untimely-catalog/',
                              catalog_index_file='untimely_index-neo7.fits')

result_table = ucx.search_by_coordinates(26.9783833, 23.6616914, box_size=100, show_result_table_in_browser=False, save_result_table=True,
                                         result_table_format='ascii.ipac', result_table_extension='dat')

# Do whatever you want with the result table data here
print(result_table.info)

ucx.create_finder_charts(overlays=True, overlay_color='green', overlay_labels=False, overlay_label_color='red',
                         image_contrast=5, open_file=False, file_format='pdf')

ucx.create_ligh_curves(photometry_radius=2, yticks=None, open_file=False, file_format='png')

ucx.create_image_blinks(blink_duration=300, image_zoom=10, image_contrast=5, scan_dir_mode=ucx.ALTERNATE_SCAN, display_blinks=False)
```

### Example output
#### Result table:
![Catalog search results](Example%20output/unTimely_Catalog_search%20results_26.978383%2B23.661691.dat)
#### Finder charts:
![Finder charts](Example%20output/unTimely_Catalog_finder_charts_26.978383%2B23.661691.png)
#### Light curves:
![Light curves](Example%20output/unTimely_Catalog_light_curves_26.978383%2B23.661691.png)
### Image blinks:
![Image blinks - variable](Example%20output/Animated_time_series_w1_26.978383%2B23.661691.gif) | ![Image blinks - color](Example%20output/Animated_time_series_26.978383%2B23.661691.gif)
![Image blinks - high PM](Example%20output/Animated_time_series_w2_133.79476-7.245146.gif)

## <kbd>class</kbd> <a id="docstrings">`unTimelyCatalogExplorer`</a>

### <kbd>constructor</kbd> `__init__`
```python
__init__(directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
         catalog_base_url='https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/',
         catalog_index_file='untimely_index-neo7.fits'):
```
Creates an unTimelyCatalogExplorer instance with the given parameters

Parameters
----------
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

Returns
-------
An unTimelyCatalogExplorer instance.

---

### <kbd>method</kbd> `search_by_coordinates`
```python
search_by_coordinates(target_ra, target_dec, box_size=100, show_result_table_in_browser=False, save_result_table=True,
                      result_table_format='ascii', result_table_extension='dat'):
```
Search the catalog by coordinates (box search).

Parameters
----------
- ra : float  
    Right ascension in decimal degrees.
- dec : float  
    Declination in decimal degrees.
- box_size : int, optional  
    Image size in arcseconds. The default is 100.
- show_result_table_in_browser : bool, optional  
    Whether to show the result table in your browser (columns can be sorted). The default is False.
- save_result_table : bool, optional  
    Whether to save the result table to the directory specified in the constructor ``unTimelyCatalogExplorer(directory=)``. The default is True.
- result_table_format : str, optional  
    Result table output format. The default is 'ascii'.
- result_table_extension : str, optional  
    Result table file extension. The default is 'dat'.

Returns
-------
Astropy table containing the catalog entries located within a field of view of the specified size at the given coordinates.

---

### <kbd>method</kbd> `create_finder_charts`
```python
create_finder_charts(overlays=True, overlay_color='green', overlay_labels=False, overlay_label_color='red',
                     image_contrast=3, open_file=False, file_format='pdf'):
```
Create finder charts for W1 and W2 at each epoch with overplotted catalog positions (overlays)

Parameters
----------
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

Raises
------
Exception if method ``search_by_coordinates`` has not been called first.

Returns
-------
None.

---

### <kbd>method</kbd> `create_ligh_curves`
```python
create_ligh_curves(photometry_radius=5, yticks=None, open_file=None, file_format=None):
```
Create light curves using W1 and W2 photometry of all available epochs.

Parameters
----------
- photometry_radius : float, optional  
    Radius to search for the photometry used to create the light curves. The default is 5.
- yticks : tuple, optional  
    Tuple containing y axis tick values. The default is pyplot's automatic tick allocation.
- open_file : bool, optional  
    Whether to open the saved light curves automatically. The default is None (value given by method ``create_finder_charts`` will be used).
- file_format : bool, optional  
    Output file format: pdf, png, eps, etc.. The default is None (value given by method ``create_finder_charts`` will be used).

Raises
------
Exception if method ``search_by_coordinates`` has not been called first.

Returns
-------
None.

---

### <kbd>method</kbd> `create_image_blinks`
```python
create_image_blinks(blink_duration=300, image_zoom=10, image_contrast=None, scan_dir_mode=ALTERNATE_SCAN, display_blinks=False):
```
Create W1 and W2 image blinks with overplotted catalog positions in GIF format.

Parameters
----------
- blink_duration : int, optional  
    Duration each image is shown in milliseconds. The default is 200.
- image_zoom : int, optional  
    Scaling factor to be applied on W1 and W2 images. The default is 10.
- image_contrast : int, optional  
    Contrast of W1 and W2 images. The default is None (value given by method ``create_finder_charts`` will be used).
- scan_dir_mode : int, optional  
    Order in which the image epochs are displayed. The default is ALTERNATE_SCAN.
    - ALTERNATE_SCAN : epoch0asc, epoch0desc, epoch1asc, epoch1desc, ...
    - SEPARATE_SCAN : epoch0asc, epoch1asc, ... epoch0desc, epoch1desc, ...
    - MERGE_SCAN : epoch0asc+epoch0desc, epoch1asc+epoch1desc, ...
- display_blinks : bool, optional  
    Whether to display the image blinks in your system's media player. The default is False.

Raises
------
Exception if method ``create_finder_charts`` has not been called first.

Returns
-------
None.

---

### References
<a id="1">[1]</a> Meisner, A. M., Caselden, D., and Schlafly, E. F., "unTimely: a Full-sky, Time-Domain unWISE Catalog", 2022. [![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/2209.14327)
