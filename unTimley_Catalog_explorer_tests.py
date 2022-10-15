import unTimley_Catalog_explorer as ucx

directory = 'C:/Users/wcq637/Documents/Private/BYW/unTimely'

#box_size = 500
#target_ra = 209.2891781
#target_dec = 55.7474398

#box_size = 200
#target_ra = 133.7947597
#target_dec = -7.2451456

#target_ra = 360
#target_dec = -89

#target_ra = 16.9701886
#target_dec = 0.6992208

target_ra = 16.9701886
target_dec = 0.6992208

box_size = 100

ucx.search_by_coordinates(target_ra, target_dec, box_size=box_size, finder_charts=True, overlays=True, overlay_color='green', overlay_labels=True,
                          overlay_label_color='red', neowise_contrast=3, show_result_table_in_browser=True, save_result_table=True, result_table_format='ascii',
                          result_table_extension='dat', open_finder_charts=True, finder_charts_format='pdf', animated_gif=True, scan_dir_mode=ucx.MERGE_SCAN,
                          directory=directory, cache=True, show_progress=True, timeout=300)
