#!/bin/bash

set -e # quit on errors
set -u # quit on undefined variables

# Copy input data to dataDir, if not already present (or if source is newer than destination)
cp -u /lustre/f1/pdata/fms/module_data/riga/sst_data.cpio                                                 $dataDir/INPUT/sst_data.cpio
cp -u /lustre/f1/pdata/fms/module_data/riga/navy_topography.data.nc                                       $dataDir/INPUT/navy_topography.data.nc
cp -u /lustre/f1/pdata/esm/ipcc_ar5/datasets/common/LM3/unpack/surface_reflectance.20090108.nc            $dataDir/INPUT/soil_brdf.nc
cp -u /lustre/f1/pdata/fms/esm/LM3/netcdf/any_grid/cover_frac_m45.20080723.nc                             $dataDir/INPUT/cover_type.nc
cp -u /lustre/f1/pdata/fms/esm/LM3/netcdf/any_grid/ground_frac_m45.20080723.nc                            $dataDir/INPUT/ground_type.nc
cp -u /lustre/f1/unswept/Sergey.Malyshev/DATA/co2/co2_rcp45_gblannualdata_1800-2500_2x2.nc                $dataDir/INPUT/co2_gblannualdata.nc
cp -u /lustre/f1/unswept/Sergey.Malyshev/DATA/CM2.1U_Control-1860_D4/biodata/2060-2109.biodata.nc         $dataDir/INPUT/biodata.nc
cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/lu.1700-2015.nc                                             $dataDir/INPUT/landuse.nc
cp -u /lustre/f1/pdata/esm/ipcc_ar5/datasets/common/LM3/unpack/surface_reflectance.20090108.nc            $dataDir/INPUT/soil_brdf.nc
cp -u /lustre/f1/pdata/esm/ipcc_ar5/datasets/common/LM3/unpack/geohydrology_table.20090108.nc             $dataDir/INPUT/geohydrology_table.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/atmos_hgrid.nc                $dataDir/INPUT/atmos_hgrid.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/atmos_mosaic.nc               $dataDir/INPUT/atmos_mosaic.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/atmos_mosaicXland_mosaic.nc   $dataDir/INPUT/atmos_mosaicXland_mosaic.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/atmos_mosaicXocean_mosaic.nc  $dataDir/INPUT/atmos_mosaicXocean_mosaic.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/land_hgrid.nc                 $dataDir/INPUT/land_hgrid.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/land_mosaic.nc                $dataDir/INPUT/land_mosaic.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/land_mosaicXocean_mosaic.nc   $dataDir/INPUT/land_mosaicXocean_mosaic.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/ocean_hgrid.nc                $dataDir/INPUT/ocean_hgrid.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/ocean_mosaic.nc               $dataDir/INPUT/ocean_mosaic.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/ocean_vgrid.nc                $dataDir/INPUT/ocean_vgrid.nc
cp -u /lustre/f1/pdata/fms/esm/gridspec/ESM2M_mosaic_feb202009.river_regrid/topog.nc                      $dataDir/INPUT/topog.nc
cp -u /ncrc/home1/Sergey.Malyshev/devel/lad2/tikal-sheffield/input/hillslope6.nc                          $dataDir/INPUT/hillslope.nc
cp -u /ncrc/home1/Sergey.Malyshev/devel/lad2/tikal-sheffield/input/coldstartsoildata_lm3zw_hillslope_Trans48-77_20130716.nc    $dataDir/INPUT/soil_wtt.nc
cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/LISOTD_LRMTS_V2.3.2013.CLIMO.20140208.nc                    $dataDir/INPUT/lightning.nc
cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/gdp_1500-2500.20141126.nc                                   $dataDir/INPUT/GDP.nc

cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/popD_HYDE100to2010.latest_highres.nc                        $dataDir/INPUT/population.nc
cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/FcFp_latest.nc                                              $dataDir/INPUT/Fk.nc

# Copy data specific to geohyrology_to_use='hill'. Don't use -u to ensure correct file always overwrites existing file.
cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/lm3_data_20150322/geohydrology.nc                           $dataDir/INPUT/geohydrology.nc
cp -u /lustre/f1/unswept/Sam.Rabin/fire_input/lm3_data_20150322/geohydrology_table_2a2n.nc                $dataDir/INPUT/geohydrology_table_2a2n.nc
