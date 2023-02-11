from unTimely_Catalog_tools import unTimelyCatalogExplorer

target_ra = 24.244145
target_dec = 47.8580834

ucx = unTimelyCatalogExplorer(directory='output', cache=False, show_progress=True, timeout=300, suppress_console_output=False,
                              catalog_base_url='http://unwise.me/data/neo7/untimely-catalog/',
                              catalog_index_file='untimely_index-neo7.fits')

result_table = ucx.search_by_coordinates(target_ra, target_dec, box_size=100, cone_radius=5, show_result_table_in_browser=False,
                                         save_result_table=True, result_table_format='ipac', result_table_extension='txt')

assert len(result_table) == 34

ucx.create_finder_charts(overlays=True, overlay_color='green', overlay_labels=True, overlay_label_color='red',
                         image_contrast=5, open_file=False, file_format='pdf')

ucx.create_light_curves(photometry_radius=2, yticks=None, open_file=False, file_format='pdf', overplot_l1b_phot=True, bin_l1b_phot=True)

ucx.create_image_blinks(blink_duration=300, image_zoom=10, image_contrast=5, separate_scan_dir=True, display_blinks=False)
