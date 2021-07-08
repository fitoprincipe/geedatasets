# coding=utf-8
""" Google Earth Engine Sentinel Collections """
from .visualization import *
from .datasets import OpticalSatellite, ImageCollection
from .bands import OpticalBand, BitBand, ClassificationBand, ExpressionBand, \
    RangeBand, Band, Precisions
from .masks import Mask
from . import helpers
import geetools
from . import register
import ee


class Modis(OpticalSatellite, ImageCollection):
    """ Modis Base Class"""
    pass


@register
class MOD09GQ(Modis):
    """ Terra Surface Reflectance Daily Global 250m """
    id = 'MODIS/006/MOD09GQ'
    short_name = 'MOD_SR_250'
    process = 'SR'
    num_observations = Band('num_observations', 'num_observations',
                            Precisions.int8, 250, 'count', 1,
                            "Number of observations per 250m pixel")
    red = OpticalBand('sur_refl_b01', 'red', wavelength=(0.62, 0.67),
                               precision=Precisions.int16, resolution=250,
                               scale=0.0001)
    nir = OpticalBand('sur_refl_b02', 'nir', wavelength=(0.841, 0.876),
                               precision=Precisions.int16, resolution=250,
                               scale=0.0001)
    QC_250m = BitBand('QC_250m', 'QA', bits={
                          '4-7': {0: 'b1_highest_quality'},
                          '8-11': {0: 'b2_highest_quality'}
                      }, resolution=250, precision=Precisions.uint16,
                      positives=['b1_highest_quality', 'b2_highest_quality']
                      )
    obscov = RangeBand('obscov', 'obscov', 0, 100, Precisions.int8, scale=0.01,
                       resolution=250, description="Observation coverage percent")
    iobs_res = RangeBand('iobs_res', 'iobs_res', 0, 254, Precisions.uint8,
                         resolution=250, description="Observation number")
    orbit_pnt = RangeBand('orbit_pnt', 'orbit_pnt', 0, 15, Precisions.uint8,
                          resolution=250, description='Orbit pointer')
    granule_pnt = RangeBand('granule_pnt', 'granule_pnt', 0, 254, Precisions.uint8,
                            resolution=250, description='Granule pointer')

    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')

    bands = (num_observations, red, nir, QC_250m, obscov, iobs_res,
             orbit_pnt, granule_pnt)

    extra_bands = (ndvi,)

    masks = (Mask.fromBand('QA', QC_250m),)


@register
class MOD11A1(Modis):
    """ Terra Land Surface Temperature and Emissivity Daily Global 1km """
    id = 'MODIS/006/MOD11A1'
    short_name = 'MOD_TEMP_1000'

    LST_Day_1km = RangeBand('LST_Day_1km', 'LST_day', 7500, 65535,
                            Precisions.uint16,  resolution=1000,
                            units='Kelvin', scale=0.02,
                            description="Daytime Land Surface Temperature")
    QC_Day = BitBand('QC_Day', 'QA', precision=Precisions.uint8,
                     resolution=1000, bits={
                          '0-1': {0: 'good_quality'}
                      }, positives=('good_quality',))
    Day_view_time = RangeBand('Day_view_time', 'hours_day', 0, 240,
                              precision=Precisions.uint8, resolution=1000,
                              units='hours', scale=0.1,
                              description="Local time of day observation")
    Day_view_angle = RangeBand('Day_view_angle', 'angle_day', 0, 130,
                               precision=Precisions.uint8,
                               resolution=1000, units='degrees',
                               description="View zenith angle of day observation")
    LST_Night_1km = RangeBand('LST_Night_1km', 'LST_night', 7500, 65535,
                              precision=Precisions.uint16,
                              resolution=1000, units='Kelvin', scale=0.02,
                              description="Nighttime Land Surface Temperature")
    QC_Night = BitBand('QC_Night', 'QA', precision=Precisions.uint8,
                       resolution=1000, bits={
                         '0-1': {0: 'good_quality'}
                       }, positives=('good_quality',))
    Night_view_time = RangeBand('Night_view_time', 'hours_night', 0, 240,
                                Precisions.uint8,
                                resolution=1000, units='hours', scale=0.1,
                                description="Local time of day observation")
    Night_view_angle = RangeBand('Night_view_angle', 'angle_night', 0, 130,
                                 Precisions.uint8,
                                 resolution=1000, units='degrees',
                                 description="View zenith angle of day observation")
    Emis_31 = OpticalBand('Emis_31', 'Emis_31', precision=Precisions.uint8,
                          scale=0.002, resolution=1000,
                          description="Band 31 emissivity")
    Emis_32 = OpticalBand('Emis_32', 'Emis_32', precision=Precisions.uint8,
                          scale=0.002, resolution=1000,
                          description="Band 32 emissivity")
    Clear_day_cov = Band('Clear_day_cov', 'Clear_day_cov',
                         precision=Precisions.uint16,
                         scale=0.0005, resolution=1000,
                         description="Day clear-sky coverage")
    Clear_night_cov = Band('Clear_night_cov', 'Clear_night_cov',
                           precision=Precisions.uint16,
                           scale=0.0005, resolution=1000,
                           description="Night clear-sky coverage")
    LST_Day_deg = ExpressionBand('LST_Day_deg', 'LST_Day_deg',
                                 "(LST_day*scale)-273.15",
                                 [LST_Day_1km], dict(scale=LST_Day_1km.scale))
    LST_Night_deg = ExpressionBand('LST_Night_deg', 'LST_Night_deg',
                                 "(LST_night*scale)-273.15",
                                 [LST_Night_1km], dict(scale=LST_Night_1km.scale))

    bands = (LST_Day_1km, QC_Day, Day_view_time, Day_view_angle, LST_Night_1km,
             QC_Night, Night_view_time, Night_view_angle, Emis_31, Emis_32,
             Clear_day_cov, Clear_night_cov)

    extra_bands = (LST_Day_deg, LST_Night_deg)

    masks = (Mask.fromBand('QC_Day', QC_Day),
             Mask.fromBand('QC_night', QC_Night))

    visualizers = Visualizers(
        LST_Day = Visualization('LST_Day', (LST_Day_1km,), 13000, 16500,
                  palette=[
                  '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
                  '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
                  '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
                  'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
                  'ff0000', 'de0101', 'c21301', 'a71001', '911003']),
        LST_Night = Visualization('LST_Night', (LST_Night_1km,), 13000, 16500,
                    palette=[
                    '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
                    '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
                    '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
                    'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
                    'ff0000', 'de0101', 'c21301', 'a71001', '911003']),
        LST_Day_deg = Visualization(
            'LST_Day_deg', (LST_Day_deg,), -13.15, 56.85, palette=[
                '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
                '0602ff', '235cb1', '307ef3', '269db1','30c8e2', '32d3ef',
                '3be285', '3ff38f', '86e26f', '3ae237','b5e22e', 'd6e21f',
                'fff705', 'ffd611', 'ffb613', 'ff8b13','ff6e08', 'ff500d',
                'ff0000', 'de0101', 'c21301', 'a71001','911003']),
        LST_Night_deg = Visualization(
            'LST_Night_deg', (LST_Night_deg,), -13.15, 56.85, palette=[
                '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
                '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
                '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
                'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
                'ff0000', 'de0101', 'c21301', 'a71001', '911003']),
    )


@register
class MOD14A1(Modis):
    id = 'MODIS/006/MOD14A1'
    short_name = 'MOD_FIRE_TERRA'
    start_date = '2000-02-18'
    end_date = None

    firemask = BitBand('FireMask', 'FireMask',
                       bits={
                           '0-3': {
                               4: 'cloud',
                               5: 'non_fire',
                               7: 'fire_low',
                               8: 'fire_mid',
                               9: 'fire_high'
                           },
                       }, positives=['fire_low', 'fire_mid', 'fire_high'])

    maxfrp = RangeBand('MaxFRP', 'MaxFRP', 0, 180000, Precisions.int32, 1000,
                       'Megawatts', 0.1, description='Maximum fire radiative power')

    sample = RangeBand('sample', 'sample', 0, 1353, Precisions.uint16, 1000,
                       description='Position of fire pixel within scan')

    qa = BitBand('QA', 'QA', bits={
        '0-1': {0: 'water', 1: 'coast', 2: 'land', 3: 'missing'},
        '2': {0: 'night', 1: 'day'}
    })

    bands = (firemask, maxfrp, sample, qa)

    masks = (Mask.fromBand('FireMask', firemask), Mask.fromBand('QA', qa))
