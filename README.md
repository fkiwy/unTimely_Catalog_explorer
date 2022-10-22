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

Most of the parameters can be omitted as they have default values (see docstrings for more details).

An instance of the unTimelyCatalogExplorer class has to be created first ```ucx = unTimelyCatalogExplorer()```.

The ```search_by_coordinates``` method must always be called before all others.

The ```create_image_blinks``` method depends on the results of the ```create_finder_charts``` method.

```
from unTimely_Catalog_tools import unTimelyCatalogExplorer
import tempfile

ucx = unTimelyCatalogExplorer(directory=tempfile.gettempdir(), cache=True, show_progress=True, timeout=300,
                              catalog_base_url='http://unwise.me/data/neo7/untimely-catalog/',
                              catalog_index_file='untimely_index-neo7.fits')

result_table = ucx.search_by_coordinates(0.12345, 0.12345, box_size=100, show_result_table_in_browser=False, save_result_table=True,
                                         result_table_format='ascii.ipac', result_table_extension='dat')

# Do whatever you want with the result table data here
print(result_table.info)

ucx.create_finder_charts(overlays=True, overlay_color='green', overlay_labels=False, overlay_label_color='red',
                         neowise_contrast=5, open_file=False, file_format='pdf')

ucx.create_ligh_curves(photometry_radius=2, open_file=False, file_format='png')

ucx.create_image_blinks(blink_duration=300, display_blinks=False, scan_dir_mode=ucx.SEPARATE_SCAN, neowise_contrast=None)
```

### References
<a id="1">[1]</a> Meisner, A. M., Caselden, D., and Schlafly, E. F., "unTimely: a Full-sky, Time-Domain unWISE Catalog", 2022. [![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/2209.14327)
