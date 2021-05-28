"""
Tests for the framework.py module
"""
import os
import numpy as np
import pytest
from mapgwm.framework import setup_model_layers, plot_cross_sections, plot_zone_maps


@pytest.fixture(scope='module')
def model_layers_zones(test_data_path, test_output_folder, delta_inset_model_grid):

    # mean dem values for each 1 km gridcell, elevation units of feet
    dem_means_raster = os.path.join(test_data_path, 'rasters', 'dem_mean_elevs_1000.tif')

    # AEM electrical resistivity-based facies classes from tempest and resolve surveys
    # (All facies classes version)
    facies_classes_netcdf = os.path.join(test_data_path, 'netcdf', 'RSTM_res_fac_depth_15m.nc')

    # Original MERAS framework (1 mi resolution), elevation units of feet
    framework_rasters = [
        os.path.join(test_data_path, 'rasters', 'vkbg_surf.tif'),  # Vicksburg-Jackson Group (top)
        os.path.join(test_data_path, 'rasters', 'ucaq_surf.tif'),  # Upper Claiborne aquifer (top)
        os.path.join(test_data_path, 'rasters', 'mccu_surf.tif'),  # Middle Claiborne confining unit (t
        os.path.join(test_data_path, 'rasters', 'mcaq_surf.tif'),  # Middle Claiborne aquifer (top)
        os.path.join(test_data_path, 'rasters', 'lccu_surf.tif'),  # Lower Claiborne confining unit (to
        os.path.join(test_data_path, 'rasters', 'lcaq_surf.tif'),  # Lower Claiborne aquifer (top)
        os.path.join(test_data_path, 'rasters', 'mwaq_surf.tif'),  # Middle Wilcox aquifer (top)
        os.path.join(test_data_path, 'rasters', 'lwaq_surf.tif'),  # Lower Wilcox aquifer (top)
        os.path.join(test_data_path, 'rasters', 'mdwy_surf.tif'),  # Midway confining unit (top)
    ]
    framework_unit_names = [
        'Undifferentiated sediments\nabove the Vicksburg',
        'Vicksburg-Jackson Group',
        'Upper Claiborne aquifer',
        'Middle Claiborne confining unit',
        'Middle Claiborne aquifer',
        'Lower Claiborne confining unit',
        'Lower Claiborne aquifer',
        'Middle Wilcox aquifer',
        'Lower Wilcox aquifer'
    ]

    output_folder = os.path.join(test_output_folder, 'framework')
    layers, zone_array = setup_model_layers(dem_means_raster,
                                            facies_classes_netcdf,
                                            framework_rasters,
                                            delta_inset_model_grid,
                                            facies_class_variable='fac_a',
                                            facies_zedge_variable='zb',
                                            dem_elevation_units='feet',
                                            framework_raster_elevation_units='feet',
                                            model_length_units='meters', output_folder=output_folder,
                                            framework_unit_names=framework_unit_names)
    return layers, zone_array


def test_setup_model_layers(model_layers_zones, test_output_folder):
    output_folder = os.path.join(test_output_folder, 'framework')
    for folder in 'figures', 'botm_array', 'zones', os.path.join('zones', 'rasters'):
        flist = os.listdir(os.path.join(output_folder, folder))
        for f in flist:
            fullpath = os.path.join(output_folder, folder, f)
            if os.path.isfile(fullpath):
                assert os.path.getsize(fullpath) > 0


def test_plot_cross_sections(model_layers_zones, test_output_folder,
                             test_data_path, delta_inset_model_grid):
    framework_unit_names = [
        'Undifferentiated sediments\nabove the Vicksburg',
        'Vicksburg-Jackson Group',
        'Upper Claiborne aquifer',
        'Middle Claiborne confining unit',
        'Middle Claiborne aquifer',
        'Lower Claiborne confining unit',
        'Lower Claiborne aquifer',
        'Middle Wilcox aquifer',
        'Lower Wilcox aquifer'
    ]
    layers, zone_array = model_layers_zones
    voxel_zones = np.unique(zone_array.ravel().astype(int))[:-len(framework_unit_names)]
    framework_zone_numbers = np.unique(zone_array.ravel().astype(int))[-len(framework_unit_names):]
    framework_unit_labels = dict(zip(framework_zone_numbers.astype(int), framework_unit_names))

    out_pdf = os.path.join(test_output_folder, 'framework', 'figures', 'x_sections.pdf')
    production_zone_top_raster = test_data_path / 'iwum/est_prod_zone_top.tif'
    production_zone_botm_raster = test_data_path / 'iwum/est_prod_zone_botm.tif'
    plot_cross_sections(layers, out_pdf, property_data=zone_array,
                        voxel_start_layer=0, voxel_zones=voxel_zones,
                        cmap='copper',
                        voxel_cmap='viridis', unit_labels=framework_unit_labels,
                        add_raster_surfaces={'est. production zone top': production_zone_top_raster,
                                             'est. production zone bottom': production_zone_botm_raster
                                             },
                        modelgrid=delta_inset_model_grid)


def test_plot_zone_maps(model_layers_zones, test_output_folder):
    layers, zone_array = model_layers_zones
    out_pdf = os.path.join(test_output_folder, 'framework', 'figures', 'maps.pdf')

    plot_zone_maps(layers, out_pdf, zones=zone_array,
                   voxel_zones=np.arange(1, 21))
    j=2
