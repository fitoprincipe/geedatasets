# coding=utf-8
""" Google Earth Engine Landsat Collections """
import ee

from .visualization import *
from .datasets import OpticalSatellite, ImageCollection
from .bands import OpticalBand, BitBand, ClassificationBand, ExpressionBand,\
                   Precisions
from .helpers import TODAY
from functools import partial
from . import register
from .masks import Mask
from .helpers import getCommonBands
import geetools
import datetime


START = {4: '1982-07-16', 5: '1984-01-01', 7: '1999-01-01',
         8: '2013-04-11', 9: '2021-10-31'}
END = {4: '1993-12-14', 5: '2012-05-05', 7: TODAY, 8: TODAY, 9: TODAY}

IDS = [
    'LANDSAT/LT04/C02/T1_L2', 'LANDSAT/LT04/C02/T2_L2',
    'LANDSAT/LT04/C02/T1_TOA', 'LANDSAT/LT04/C02/T2_TOA',
    'LANDSAT/LT05/C02/T1_L2', 'LANDSAT/LT05/C02/T2_L2',
    'LANDSAT/LT05/C02/T1_TOA', 'LANDSAT/LT05/C02/T2_TOA',
    'LANDSAT/LE07/C02/T1_L2', 'LANDSAT/LE07/C02/T2_L2',
    'LANDSAT/LE07/C02/T1_TOA', 'LANDSAT/LE07/C02/T2_TOA',
    'LANDSAT/LC08/C02/T1_L2', 'LANDSAT/LC08/C02/T2_L2',
    'LANDSAT/LC08/C02/T1_TOA', 'LANDSAT/LC08/C02/T2_TOA',
    'LANDSAT/LC09/C02/T1_L2', 'LANDSAT/LC09/C02/T2_L2',
    'LANDSAT/LC09/C02/T1_TOA', 'LANDSAT/LC09/C02/T2_TOA',
]
SLC_OFF = '2003-05-31'


class Landsat(OpticalSatellite, ImageCollection):
    """  Landsat Collection """
    INFO = """ID: {id}
Short Name: {short_name}
Spacecraft: {spacecraft}
Number: {number}
Sensor: {sensor}
Process: {process}
Start Date: {start_date}
End Date: {end_date}
Cloud Cover: {cloud_cover}
Tier: {tier}
Bands: {bands}
Masks: {masks}
visualizers: {visualizers}
"""
    spacecraft = 'LANDSAT'
    cloud_cover = 'CLOUD_COVER'
    number = None
    tier = None
    process = None
    sensor = None

    qa_pixel = BitBand(
        'QA_PIXEL', 'qa_pixel', precision=Precisions.uint16, resolution=30,
        bits={
            '0': {1: 'fill'},
            '1': {1: 'dilated_cloud'},
            '3': {1: 'cloud'},
            '4': {1: 'shadow'},
            '5': {1: 'snow'},
            '6': {1: 'clear'},
            '7': {1: 'water'},
            '8-9': {1: 'cloud_confidence_low', 2: 'cloud_confidence_medium',
                    3: 'cloud_confidence_high'},
            '10-11': {1: 'shadow_confidence_low',
                      2: 'shadow_confidence_medium',
                      3: 'shadow_confidence_high'},
            '12-13': {1: 'snow_confidence_low', 2: 'snow_confidence_medium',
                      3: 'snow_confidence_high'},
            '14-15': {1: 'cirrus_confidence_low',
                      2: 'cirrus_confidence_medium',
                      3: 'cirrus_confidence_high'},
        },
        positives=['clear'],
        negatives=['water', 'shadow', 'snow', 'cloud', 'cloud_confidence_high',
                   'cirrus_confidence_high']
    )

    def __init__(self, **kwargs):
        super(Landsat, self).__init__(**kwargs)


class Tier1:
    tier = 1


class Tier2:
    tier = 2


class SR:
    process = 'SR'
    atran = OpticalBand('ST_ATRAN', 'transmittance', units='Kelvin',
                        precision=Precisions.int16, scale=0.0001,
                        resolution=30, description='Atmospheric Transmittance')
    cdist = OpticalBand('ST_CDIST', 'cloud_distance', units='km',
                        precision=Precisions.int16, scale=0.01,
                        resolution=30, description='Distance to cloud mask')
    drad = OpticalBand('ST_DRAD', 'downwelled', units='W/(m^2*sr*um)/DN',
                       precision=Precisions.int16, scale=0.001,
                       resolution=30, description="Downwelled Radiance")
    emis = OpticalBand('ST_EMIS', 'emissivity', scale=0.0001,
                       precision=Precisions.int16, resolution=30,
                       description='Emissivity of Band 10 estimated from ASTER GED')
    emsd = OpticalBand('ST_EMSD', 'emissivity_std', scale=0.0001,
                       precision=Precisions.int16, resolution=30,
                       description='Emissivity standard deviation')
    stqa = OpticalBand('ST_QA', 'st_qa', units='K', scale=0.01,
                       resolution=30, precision=Precisions.int16,
                       description='Uncertainty of the Surface Temperature band')
    trad = OpticalBand('ST_TRAD', 'thermal_radiance', units='W/(m^2*sr*um)/DN',
                       precision=Precisions.int16, scale=0.001, resolution=30,
                       description='Thermal band converted to radiance')
    urad = OpticalBand('ST_URAD', 'upwelled', units='W/(m^2*sr*um)/DN',
                       precision=Precisions.int16, scale=0.001, resolution=30,
                       description='Upwelled radiance')
    opacity = OpticalBand('SR_ATMOS_OPACITY', 'opacity', scale=0.001,
                          precision=Precisions.int16, resolution=30,
                          description='Atmospheric Opacity, clear: value < 0.1,'
                                      'averge: 0.1 < value < 0.3, '
                                      'haze: value > 0.3')
    sr_cloud_qa = BitBand(
        name='SR_CLOUD_QA', alias='sr_cloud_qa',
        precision=Precisions.uint8, resolution=30,
        bits={
            '0': {1: 'ddv'},
            '1': {1: 'cloud'},
            '2': {1: 'shadow'},
            '3': {1: 'adjacent'},
            '4': {1: 'snow'},
            '5': {1: 'water'}
        },
        negatives=['ddv', 'cloud', 'shadow', 'snow']
    )    


class TOA:
    process = 'TOA'
    _angles = dict(precision=Precisions.int16, resolution=30, units='degrees')
    saa = OpticalBand('SAA', 'solar_azimuth',
                      description='Solar Azimuth Angle', **_angles)
    sza = OpticalBand('SZA', 'solar_zenith',
                      description='Solar Zenith Angle', **_angles)
    vaa = OpticalBand('VAA', 'view_azimuth',
                      description='View Azimuth Angle', **_angles)
    vza = OpticalBand('VZA', 'view_zenith',
                      description='View Zenith Angle', **_angles)
    pan = OpticalBand('B8', 'pan', wavelength=(0.52, 0.9),
                      precision=Precisions.float, resolution=15,
                      description='Panchromatic')


class RAW:
    process = 'RAW'
    saa = TOA.saa
    sza = TOA.sza
    vaa = TOA.vaa
    vza = TOA.vza
    pan = OpticalBand('B8', 'pan', wavelength=(0.52, 0.9),
                      precision=Precisions.uint8, resolution=15,
                      description='Panchromatic')


class MSS:
    green = OpticalBand(
        'B4', 'green', wavelength=(0.5, 0.6), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    red = OpticalBand(
        'B5', 'red', wavelength=(0.6, 0.7), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    nir = OpticalBand(
        'B6', 'nir', wavelength=(0.7, 0.8), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    swir = OpticalBand(
        'B7', 'swir', wavelength=(0.8, 1.1), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    qa_radsat = BitBand(
        'QA_RADSAT', 'qa_radsat', precision=Precisions.uint16, resolution=30,
        bits={
            0: {1: 'B1_saturated'},
            1: {1: 'B2_saturated'},
            2: {1: 'B3_saturated'},
            3: {1: 'B4_saturated'},
            4: {1: 'B5_saturated'},
            5: {1: 'B6_saturated'},
            6: {1: 'B7_saturated'},
            9: {1: 'dropped'}
        }
    )


class MSS45:
    green = OpticalBand(
        'B1', 'green', wavelength=(0.5, 0.6), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    red = OpticalBand(
        'B2', 'red', wavelength=(0.6, 0.7), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    nir = OpticalBand(
        'B3', 'nir', wavelength=(0.7, 0.8), precision=Precisions.uint8,
        resolution=80, scale=1/255)
    swir = OpticalBand(
        'B4', 'swir', wavelength=(0.8, 1.1), precision=Precisions.uint8,
        resolution=80, scale=1/255)


class TM:
    qa_radsat = BitBand(
        'QA_RADSAT', 'qa_radsat', precision=Precisions.uint16, resolution=30,
        bits={
            0: {1: 'B1_saturated'},
            1: {1: 'B2_saturated'},
            2: {1: 'B3_saturated'},
            3: {1: 'B4_saturated'},
            4: {1: 'B5_saturated'},
            5: {1: 'B6_saturated'},
            6: {1: 'B7_saturated'},
            9: {1: 'dropped'}
        }
    )


class ETM:
    qa_radsat = BitBand(
        'QA_RADSAT', 'qa_radsat', precision=Precisions.uint16, resolution=30,
        bits={
            0: {1: 'B1_saturated'},
            1: {1: 'B2_saturated'},
            2: {1: 'B3_saturated'},
            3: {1: 'B4_saturated'},
            4: {1: 'B5_saturated'},
            5: {1: 'B6L_saturated'},
            6: {1: 'B7_saturated'},
            8: {1: 'B6H_saturated'},
            9: {1: 'dropped'}
        }
    )


class OLI:
    qa_radsat = BitBand(
        'QA_RADSAT', 'qa_radsat', precision=Precisions.uint16, resolution=30,
        bits={
            0: {1: 'B1_saturated'},
            1: {1: 'B2_saturated'},
            2: {1: 'B3_saturated'},
            3: {1: 'B4_saturated'},
            4: {1: 'B5_saturated'},
            5: {1: 'B6_saturated'},
            6: {1: 'B7_saturated'},
            8: {1: 'B9_saturated'},
            11: {1: 'terrain_occlusion'}
        }
    )


class TM_SR:
    sensor = 'TM'
    _optical = dict(resolution=30, precision=Precisions.uint16,
                    scale=2.75e-05, offset=-0.2)
    blue = OpticalBand('SR_B1', 'blue', wavelength=(0.45, 0.52), **_optical)
    green = OpticalBand('SR_B2', 'green', wavelength=(0.52, 0.6), **_optical)
    red = OpticalBand('SR_B3', 'red', wavelength=(0.63, 0.69), **_optical)
    nir = OpticalBand('SR_B4', 'nir', wavelength=(0.76, 0.9), **_optical)
    swir = OpticalBand('SR_B5', 'swir', wavelength=(1.55, 1.75), **_optical)
    thermal = OpticalBand('ST_B6', 'thermal', units='Kelvin',
                          precision=Precisions.uint16, resolution=30,
                          wavelength=(10.4, 12.5), scale=0.00341802,
                          offset=149)
    swir2 = OpticalBand('SR_B7', 'swir2', wavelength=(2.08, 2.35), **_optical)


class ETM_SR:
    sensor = 'ETM+'
    blue = TM_SR.blue
    green = TM_SR.green
    red = TM_SR.red
    nir = TM_SR.nir
    swir = TM_SR.swir
    swir2 = TM_SR.swir2
    thermal = TM_SR.thermal


class OLI_SR:
    sensor = 'OLI'
    _optical = dict(resolution=30, precision=Precisions.uint16,
                    scale=2.75e-05, offset=-0.2)
    ublue = OpticalBand('SR_B1', 'ublue', wavelength=(0.435, 0.451), **_optical)
    blue = OpticalBand('SR_B2', 'blue', wavelength=(0.452, 0.512), **_optical)
    green = OpticalBand('SR_B3', 'green', wavelength=(0.533, 0.59), **_optical)
    red = OpticalBand('SR_B4', 'red', wavelength=(0.636, 0.673), **_optical)
    nir = OpticalBand('SR_B5', 'nir', wavelength=(0.851, 0.879), **_optical)
    swir = OpticalBand('SR_B6', 'swir', wavelength=(1.566, 1.651), **_optical)
    swir2 = OpticalBand('SR_B7', 'swir2', wavelength=(2.107, 2.294), **_optical)
    thermal = OpticalBand('ST_B10', 'thermal', units='Kelvin',
                          precision=Precisions.uint16, resolution=30,
                          wavelength=(10.6, 11.19), scale=0.00341802,
                          offset=149)
    sr_qa_aerosol = BitBand(
        'SR_QA_AEROSOL',
        'sr_qa_aerosol', precision=Precisions.uint8, resolution=30,
        bits={
            '0': {1: 'fill'},
            '1': {1: 'valid'},
            '2': {1: 'water'},
            '5': {1: 'interpolated'},
            '6-7': {0: 'climatology',
                    1: 'low', 2: 'medium',
                    3: 'high'}
        }
    )


class TM_TOA:
    sensor = 'TM'
    _optical = dict(resolution=30, precision=Precisions.float,
                    scale=1, offset=0)
    blue = OpticalBand('B1', 'blue', wavelength=(0.45, 0.52), **_optical)
    green = OpticalBand('B2', 'green', wavelength=(0.52, 0.6), **_optical)
    red = OpticalBand('B3', 'red', wavelength=(0.63, 0.69), **_optical)
    nir = OpticalBand('B4', 'nir', wavelength=(0.76, 0.9), **_optical)
    swir = OpticalBand('B5', 'swir', wavelength=(1.55, 1.75), **_optical)
    swir2 = OpticalBand('B7', 'swir2', wavelength=(2.08, 2.35), **_optical)
    thermal = OpticalBand('B6', 'thermal', units='Kelvin',
                          precision=Precisions.uint16, resolution=30,
                          wavelength=(10.4, 12.5), scale=0.00341802,
                          offset=149)


class ETM_TOA:
    sensor = 'ETM+'
    blue = TM_TOA.blue
    green = TM_TOA.green
    red = TM_TOA.red
    nir = TM_TOA.nir
    swir = TM_TOA.swir
    swir2 = TM_TOA.swir2


class OLI_TOA:
    sensor = 'OLI'
    _optical = dict(resolution=30, precision=Precisions.float)
    ublue = OpticalBand('B1', 'ublue', wavelength=(0.43, 0.45), **_optical)
    blue = OpticalBand('B2', 'blue', wavelength=(0.45, 0.51), **_optical)
    green = OpticalBand('B3', 'green', wavelength=(0.53, 0.59), **_optical)
    red = OpticalBand('B4', 'red', wavelength=(0.64, 0.67), **_optical)
    nir = OpticalBand('B5', 'nir', wavelength=(0.85, 0.88), **_optical)
    swir = OpticalBand('B6', 'swir', wavelength=(1.566, 1.651), **_optical)
    swir2 = OpticalBand('B7', 'swir2', wavelength=(2.107, 2.294), **_optical)
    cirrus = OpticalBand('B9', 'cirrus', wavelength=(1.36, 1.38),
                         resolution=15, precision=Precisions.float)
    thermal = OpticalBand('B10', 'thermal', units='Kelvin',
                          precision=Precisions.float, resolution=30,
                          wavelength=(10.6, 11.19))
    thermal2 = OpticalBand('B11', 'thermal2', units='Kelvin',
                          precision=Precisions.float, resolution=30,
                          wavelength=(11.5, 12.51))


class TM_RAW:
    sensor = 'TM'
    _optical = dict(resolution=30, precision=Precisions.uint8,
                    scale=1/255, offset=0)
    blue = OpticalBand('B1', 'blue', wavelength=(0.45, 0.52), **_optical)
    green = OpticalBand('B2', 'green', wavelength=(0.52, 0.6), **_optical)
    red = OpticalBand('B3', 'red', wavelength=(0.63, 0.69), **_optical)
    nir = OpticalBand('B4', 'nir', wavelength=(0.76, 0.9), **_optical)
    swir = OpticalBand('B5', 'swir', wavelength=(1.55, 1.75), **_optical)
    swir2 = OpticalBand('B7', 'swir2', wavelength=(2.08, 2.35), **_optical)
    thermal = OpticalBand('B6', 'thermal', units='Kelvin',
                          precision=Precisions.uint16, resolution=30,
                          wavelength=(10.4, 12.5), scale=0.00341802,
                          offset=149)


class ETM_RAW:
    sensor = 'ETM+'
    blue = TM_RAW.blue
    green = TM_RAW.green
    red = TM_RAW.red
    nir = TM_RAW.nir
    swir = TM_RAW.swir
    swir2 = TM_RAW.swir2


class OLI_RAW:
    sensor = 'OLI'
    _optical = dict(resolution=30, precision=Precisions.uint16, scale=1/65535)
    ublue = OpticalBand('B1', 'ublue', wavelength=(0.43, 0.45), **_optical)
    blue = OpticalBand('B2', 'blue', wavelength=(0.45, 0.51), **_optical)
    green = OpticalBand('B3', 'green', wavelength=(0.53, 0.59), **_optical)
    red = OpticalBand('B4', 'red', wavelength=(0.64, 0.67), **_optical)
    nir = OpticalBand('B5', 'nir', wavelength=(0.85, 0.88), **_optical)
    swir = OpticalBand('B6', 'swir', wavelength=(1.566, 1.651), **_optical)
    swir2 = OpticalBand('B7', 'swir2', wavelength=(2.107, 2.294), **_optical)
    cirrus = OpticalBand('B9', 'cirrus', wavelength=(1.36, 1.38),
                         resolution=15, precision=Precisions.uint16,
                         scale=1/65535)
    thermal = OpticalBand('B10', 'thermal', units='Kelvin',
                          precision=Precisions.uint16, resolution=30,
                          wavelength=(10.6, 11.19), scale=1/65535)
    thermal2 = OpticalBand('B11', 'thermal2', units='Kelvin',
                           precision=Precisions.uint16, resolution=30,
                           wavelength=(11.5, 12.51), scale=1/65535)


@register
class Landsat1RAW(Landsat, Tier1, MSS):
    id = 'LANDSAT/LM01/C02/T1'
    short_name = 'L1RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    bands = (MSS.green, MSS.red, MSS.nir, MSS.swir, Landsat.qa_pixel,
             MSS.qa_radsat)
    visualizers = Visualizers(
        FalseColor=Visualization.falseColor([MSS.nir, MSS.red, MSS.green])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [MSS.nir, MSS.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [MSS.nir, MSS.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)


@register
class Landsat1RAWT2(Tier2, Landsat1RAW):
    id = 'LANDSAT/LM01/C02/T2'
    short_name = 'L1RAWT2'

    def __init__(self, **kwargs):
        super(Landsat1RAWT2, self).__init__(**kwargs)


@register
class Landsat2RAW(Landsat, Tier1, MSS):
    id = 'LANDSAT/LM02/C02/T1'
    short_name = 'L2RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    bands = (MSS.green, MSS.red, MSS.nir, MSS.swir, Landsat.qa_pixel,
             MSS.qa_radsat)
    visualizers = Visualizers(
        FalseColor=Visualization.falseColor([MSS.nir, MSS.red, MSS.green])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [MSS.nir, MSS.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [MSS.nir, MSS.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)


@register
class Landsat2RAWT2(Tier2, Landsat2RAW):
    id = 'LANDSAT/LM02/C02/T2'
    short_name = 'L2RAWT2'

    def __init__(self, **kwargs):
        super(Landsat2RAWT2, self).__init__(**kwargs)


@register
class Landsat3RAW(Landsat, Tier1, MSS):
    id = 'LANDSAT/LM03/C02/T1'
    short_name = 'L3RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    bands = (MSS.green, MSS.red, MSS.nir, MSS.swir, Landsat.qa_pixel,
             MSS.qa_radsat)
    visualizers = Visualizers(
        FalseColor=Visualization.falseColor([MSS.nir, MSS.red, MSS.green])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [MSS.nir, MSS.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [MSS.nir, MSS.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)


@register
class Landsat3RAWT2(Tier2, Landsat3RAW):
    id = 'LANDSAT/LM03/C02/T2'
    short_name = 'L3RAWT2'

    def __init__(self, **kwargs):
        super(Landsat3RAWT2, self).__init__(**kwargs)


@register
class Landsat4MSS(Landsat, Tier1, MSS45):
    id = 'LANDSAT/LM04/C02/T1'
    short_name = 'L4MSS'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    bands = (MSS45.green, MSS45.red, MSS45.nir, MSS45.swir, Landsat.qa_pixel,
             MSS.qa_radsat)
    visualizers = Visualizers(
        FalseColor=Visualization.falseColor([MSS45.nir, MSS45.red, MSS45.green])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [MSS45.nir, MSS45.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [MSS45.nir, MSS45.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)


@register
class Landsat4MSST2(Tier2, Landsat4MSS):
    id = 'LANDSAT/LM04/C02/T2'
    short_name = 'L4MSST2'

    def __init__(self, **kwargs):
        super(Landsat4MSST2, self).__init__(**kwargs)


@register
class Landsat5MSS(Landsat, Tier1, MSS45):
    id = 'LANDSAT/LM05/C02/T1'
    short_name = 'L5MSS'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    bands = (MSS45.green, MSS45.red, MSS45.nir, MSS45.swir, Landsat.qa_pixel,
             MSS.qa_radsat)
    visualizers = Visualizers(
        FalseColor=Visualization.falseColor([MSS45.nir, MSS45.red, MSS45.green])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [MSS45.nir, MSS45.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [MSS45.nir, MSS45.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)


@register
class Landsat5MSST2(Tier2, Landsat5MSS):
    id = 'LANDSAT/LM05/C02/T2'
    short_name = 'L5MSST2'

    def __init__(self, **kwargs):
        super(Landsat5MSST2, self).__init__(**kwargs)


@register
class Landsat4SR(Landsat, Tier1, TM_SR, SR):
    id = 'LANDSAT/LT04/C02/T1_L2'
    short_name = 'L4SR'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),
             Mask.fromBand('sr_cloud_qa', SR.sr_cloud_qa))

    bands = (TM_SR.blue, TM_SR.green, TM_SR.red, TM_SR.nir, TM_SR.swir,
             TM_SR.swir2, TM_SR.thermal, SR.atran, SR.cdist,
             SR.drad, SR.emis, SR.emsd, SR.stqa, SR.trad, SR.urad, SR.opacity,
             SR.sr_cloud_qa, Landsat.qa_pixel, TM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([TM_SR.red, TM_SR.green, TM_SR.blue]),
        FalseColor = Visualization.falseColor([TM_SR.nir, TM_SR.red, TM_SR.green]),
        NSR = Visualization.NSR([TM_SR.nir, TM_SR.swir, TM_SR.red]),
        AtmosphericOpacity = SR.opacity
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [TM_SR.nir, TM_SR.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [TM_SR.nir, TM_SR.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat4SR, self).__init__(**kwargs)


@register
class Landsat4SRT2(Tier2, Landsat4SR):
    id = 'LANDSAT/LT04/C02/T2_L2'
    short_name = 'L4SRT2'

    def __init__(self, **kwargs):
        super(Landsat4SRT2, self).__init__(**kwargs)


@register
class Landsat5SR(Landsat, Tier1, TM_SR, SR):
    start_date = '1984-03-16'
    end_date = '2012-05-05'
    id = 'LANDSAT/LT05/C02/T1_L2'
    short_name = 'L5SR'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),
             Mask.fromBand('sr_cloud_qa', SR.sr_cloud_qa))

    bands = (TM_SR.blue, TM_SR.green, TM_SR.red, TM_SR.nir, TM_SR.swir,
             TM_SR.swir2, TM_SR.thermal, SR.atran, SR.cdist,
             SR.drad, SR.emis, SR.emsd, SR.stqa, SR.trad, SR.urad, SR.opacity,
             SR.sr_cloud_qa, Landsat.qa_pixel, TM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([TM_SR.red, TM_SR.green, TM_SR.blue]),
        FalseColor = Visualization.falseColor([TM_SR.nir, TM_SR.red, TM_SR.green]),
        NSR = Visualization.NSR([TM_SR.nir, TM_SR.swir, TM_SR.red]),
        AtmosphericOpacity = SR.opacity
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [TM_SR.nir, TM_SR.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [TM_SR.nir, TM_SR.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat5SR, self).__init__(**kwargs)


@register
class Landsat5SRT2(Tier2, Landsat5SR):
    id = 'LANDSAT/LT05/C02/T2_L2'
    short_name = 'L5SRT2'

    def __init__(self, **kwargs):
        super(Landsat5SRT2, self).__init__(**kwargs)


@register
class Landsat7SR(Landsat, Tier1, ETM_SR, SR):
    start_date = '1999-05-28'
    end_date = '2022-04-06'
    id = 'LANDSAT/LE07/C02/T1_L2'
    short_name = 'L7SR'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),
             Mask.fromBand('sr_cloud_qa', SR.sr_cloud_qa))

    bands = (ETM_SR.blue, ETM_SR.green, ETM_SR.red, ETM_SR.nir, ETM_SR.swir,
             ETM_SR.swir2, ETM_SR.thermal, SR.atran, SR.cdist,
             SR.drad, SR.emis, SR.emsd, SR.stqa, SR.trad, SR.urad, SR.opacity,
             SR.sr_cloud_qa, Landsat.qa_pixel, ETM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([ETM_SR.red, ETM_SR.green, ETM_SR.blue]),
        FalseColor = Visualization.falseColor([ETM_SR.nir, ETM_SR.red, ETM_SR.green]),
        NSR = Visualization.NSR([ETM_SR.nir, ETM_SR.swir, ETM_SR.red]),
        AtmosphericOpacity = SR.opacity
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [ETM_SR.nir, ETM_SR.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [ETM_SR.nir, ETM_SR.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat7SR, self).__init__(**kwargs)


@register
class Landsat7SRT2(Tier2, Landsat7SR):
    id = 'LANDSAT/LE07/C02/T2_L2'
    short_name = 'L7SRT2'

    def __init__(self, **kwargs):
        super(Landsat7SRT2, self).__init__(**kwargs)


@register
class Landsat8SR(Landsat, Tier1, OLI_SR, SR):
    start_date = '2013-03-18'
    end_date = None
    id = 'LANDSAT/LC08/C02/T1_L2'
    short_name = 'L8SR'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)

    bands = (OLI_SR.ublue, OLI_SR.blue, OLI_SR.green, OLI_SR.red, OLI_SR.nir,
             OLI_SR.swir, OLI_SR.swir2, OLI_SR.thermal, SR.atran, SR.cdist,
             SR.drad, SR.emis, SR.emsd, SR.stqa, SR.trad, SR.urad,
             Landsat.qa_pixel, OLI_SR.sr_qa_aerosol, OLI.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([OLI_SR.red, OLI_SR.green, OLI_SR.blue]),
        FalseColor = Visualization.falseColor([OLI_SR.nir, OLI_SR.red, OLI_SR.green]),
        NSR = Visualization.NSR([OLI_SR.nir, OLI_SR.swir, OLI_SR.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [OLI_SR.nir, OLI_SR.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [OLI_SR.nir, OLI_SR.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat8SR, self).__init__(**kwargs)


@register
class Landsat8SRT2(Tier2, Landsat8SR):
    id = 'LANDSAT/LC08/C02/T2_L2'
    short_name = 'L8SRT2'

    def __init__(self, **kwargs):
        super(Landsat8SRT2, self).__init__(**kwargs)


@register
class Landsat9SR(Landsat, Tier1, OLI_SR, SR):
    start_date = '2021-10-31'
    end_date = None
    id = 'LANDSAT/LC09/C02/T1_L2'
    short_name = 'L9SR'

    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)

    bands = (OLI_SR.ublue, OLI_SR.blue, OLI_SR.green, OLI_SR.red, OLI_SR.nir,
             OLI_SR.swir, OLI_SR.swir2, OLI_SR.thermal, SR.atran, SR.cdist,
             SR.drad, SR.emis, SR.emsd, SR.stqa, SR.trad, SR.urad,
             Landsat.qa_pixel, OLI_SR.sr_qa_aerosol, OLI.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([OLI_SR.red, OLI_SR.green, OLI_SR.blue]),
        FalseColor = Visualization.falseColor([OLI_SR.nir, OLI_SR.red, OLI_SR.green]),
        NSR = Visualization.NSR([OLI_SR.nir, OLI_SR.swir, OLI_SR.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [OLI_SR.nir, OLI_SR.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [OLI_SR.nir, OLI_SR.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat9SR, self).__init__(**kwargs)


@register
class Landsat9SRT2(Tier2, Landsat9SR):
    id = 'LANDSAT/LC09/C02/T2_L2'
    short_name = 'L9SRT2'

    def __init__(self, **kwargs):
        super(Landsat9SRT2, self).__init__(**kwargs)


@register
class Landsat5TOA(Landsat, Tier1, TM_TOA, TOA):
    start_date = '1984-04-19'
    end_date = '2011-11-08'
    id = 'LANDSAT/LT05/C02/T1_TOA'
    short_name = 'L5TOA'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)

    bands = (TM_TOA.blue, TM_TOA.green, TM_TOA.red, TM_TOA.nir, TM_TOA.swir,
             TM_TOA.swir2, TM_TOA.thermal, Landsat.qa_pixel, TOA.saa, TOA.sza,
             TOA.vaa, TOA.vza, TM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([TM_TOA.red, TM_TOA.green, TM_TOA.blue]),
        FalseColor = Visualization.falseColor([TM_TOA.nir, TM_TOA.red, TM_TOA.green]),
        NSR = Visualization.NSR([TM_TOA.nir, TM_TOA.swir, TM_TOA.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [TM_TOA.nir, TM_TOA.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [TM_TOA.nir, TM_TOA.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat5TOA, self).__init__(**kwargs)


@register
class Landsat5TOAT2(Tier2, Landsat5TOA):
    id = 'LANDSAT/LT05/C02/T2_TOA'
    short_name = 'L5TOAT2'

    def __init__(self, **kwargs):
        super(Landsat5TOAT2, self).__init__(**kwargs)


@register
class Landsat7TOA(Landsat, Tier1, ETM_TOA, TOA):
    start_date = '1999-06-29'
    end_date = '2022-04-06'
    id = 'LANDSAT/LE07/C02/T1_TOA'
    short_name = 'L7TOA'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    thermal_low = OpticalBand(
        'B6_VCID_1', 'thermal_low', wavelength=(10.4, 12.5),
        precision=Precisions.float, resolution=60,
        description='Low-gain Thermal Infrared 1'
    )
    thermal_high = OpticalBand(
        'B6_VCID_2', 'thermal_high', wavelength=(10.4, 12.5),
        precision=Precisions.float, resolution=60,
        description='High-gain Thermal Infrared 1'
    )

    bands = (ETM_TOA.blue, ETM_TOA.green, ETM_TOA.red, ETM_TOA.nir, ETM_TOA.swir,
             ETM_TOA.swir2, thermal_low, thermal_high, TOA.pan,
             Landsat.qa_pixel, TOA.saa, TOA.sza, TOA.vaa, TOA.vza,
             ETM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([ETM_TOA.red, ETM_TOA.green, ETM_TOA.blue]),
        FalseColor = Visualization.falseColor([ETM_TOA.nir, ETM_TOA.red, ETM_TOA.green]),
        NSR = Visualization.NSR([ETM_TOA.nir, ETM_TOA.swir, ETM_TOA.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [ETM_TOA.nir, ETM_TOA.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [ETM_TOA.nir, ETM_TOA.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat7TOA, self).__init__(**kwargs)


@register
class Landsat7TOAT2(Tier2, Landsat7TOA):
    id = 'LANDSAT/LE07/C02/T2_TOA'
    short_name = 'L7TOAT2'

    def __init__(self, **kwargs):
        super(Landsat7TOAT2, self).__init__(**kwargs)


@register
class Landsat8TOA(Landsat, Tier1, OLI_TOA, TOA):
    start_date = '2013-03-18'
    end_date = None
    id = 'LANDSAT/LC08/C02/T1_TOA'
    short_name = 'L8TOA'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    bands = (OLI_TOA.blue, OLI_TOA.green, OLI_TOA.red, OLI_TOA.nir, OLI_TOA.swir,
             OLI_TOA.swir2, OLI_TOA.thermal, OLI_TOA.thermal2, OLI_TOA.cirrus,
             Landsat.qa_pixel, TOA.saa, TOA.sza, TOA.vaa, TOA.vza,
             OLI.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([OLI_TOA.red, OLI_TOA.green, OLI_TOA.blue]),
        FalseColor = Visualization.falseColor([OLI_TOA.nir, OLI_TOA.red, OLI_TOA.green]),
        NSR = Visualization.NSR([OLI_TOA.nir, OLI_TOA.swir, OLI_TOA.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [OLI_TOA.nir, OLI_TOA.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [OLI_TOA.nir, OLI_TOA.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat8TOA, self).__init__(**kwargs)


@register
class Landsat8TOAT2(Tier2, Landsat8TOA):
    id = 'LANDSAT/LC08/C02/T2_TOA'
    short_name = 'L8TOAT2'

    def __init__(self, **kwargs):
        super(Landsat8TOAT2, self).__init__(**kwargs)


@register
class Landsat9TOA(Landsat, Tier1, OLI_TOA, TOA):
    start_date = '2021-10-31'
    end_date = None
    id = 'LANDSAT/LC09/C02/T1_TOA'
    short_name = 'L9TOA'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)

    bands = (OLI_TOA.blue, OLI_TOA.green, OLI_TOA.red, OLI_TOA.nir, OLI_TOA.swir,
             OLI_TOA.swir2, OLI_TOA.thermal, OLI_TOA.thermal2, OLI_TOA.cirrus,
             Landsat.qa_pixel, TOA.saa, TOA.sza, TOA.vaa, TOA.vza, OLI.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([OLI_TOA.red, OLI_TOA.green, OLI_TOA.blue]),
        FalseColor = Visualization.falseColor([OLI_TOA.nir, OLI_TOA.red, OLI_TOA.green]),
        NSR = Visualization.NSR([OLI_TOA.nir, OLI_TOA.swir, OLI_TOA.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [OLI_TOA.nir, OLI_TOA.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [OLI_TOA.nir, OLI_TOA.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat9TOA, self).__init__(**kwargs)


@register
class Landsat9TOAT2(Tier2, Landsat9TOA):
    id = 'LANDSAT/LC09/C02/T2_TOA'
    short_name = 'L9TOAT2'

    def __init__(self, **kwargs):
        super(Landsat9TOAT2, self).__init__(**kwargs)


@register
class Landsat5RAW(Landsat, Tier1, TM_RAW, RAW):
    start_date = '1984-01-01'
    end_date = '2012-05-05'
    id = 'LANDSAT/LT05/C02/T1'
    short_name = 'L5RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)    

    bands = (TM_RAW.blue, TM_RAW.green, TM_RAW.red, TM_RAW.nir, TM_RAW.swir,
             TM_RAW.thermal, TM_RAW.swir2, Landsat.qa_pixel, RAW.saa, RAW.sza,
             RAW.vaa, RAW.vza, TM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([TM_RAW.red, TM_RAW.green, TM_RAW.blue]),
        FalseColor = Visualization.falseColor([TM_RAW.nir, TM_RAW.red, TM_RAW.green]),
        NSR = Visualization.NSR([TM_RAW.nir, TM_RAW.swir, TM_RAW.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [TM_RAW.nir, TM_RAW.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [TM_RAW.nir, TM_RAW.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat5RAW, self).__init__(**kwargs)


@register
class Landsat5RAWT2(Tier2, Landsat5RAW):
    id = 'LANDSAT/LT05/C02/T2'
    short_name = 'L5RAWT2'

    def __init__(self, **kwargs):
        super(Landsat5RAWT2, self).__init__(**kwargs)


@register
class Landsat7RAW(Landsat, Tier1, ETM_RAW, RAW):
    start_date = '1999-01-01'
    end_date = '2022-04-06'
    id = 'LANDSAT/LE07/C02/T1'
    short_name = 'L7RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)
    thermal_low = OpticalBand(
        'B6_VCID_1', 'thermal_low', wavelength=(10.4, 12.5),
        precision=Precisions.uint8, resolution=60,
        description='Low-gain Thermal Infrared 1'
    )
    thermal_high = OpticalBand(
        'B6_VCID_2', 'thermal_high', wavelength=(10.4, 12.5),
        precision=Precisions.uint8, resolution=60,
        description='High-gain Thermal Infrared 1'
    )

    bands = (ETM_RAW.blue, ETM_RAW.green, ETM_RAW.red, ETM_RAW.nir, ETM_RAW.swir,
             ETM_RAW.swir2, thermal_low, thermal_high, RAW.pan,
             Landsat.qa_pixel, RAW.saa, RAW.sza, RAW.vaa, RAW.vza,
             ETM.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([ETM_RAW.red, ETM_RAW.green, ETM_RAW.blue]),
        FalseColor = Visualization.falseColor([ETM_RAW.nir, ETM_RAW.red, ETM_RAW.green]),
        NSR = Visualization.NSR([ETM_RAW.nir, ETM_RAW.swir, ETM_RAW.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [ETM_RAW.nir, ETM_RAW.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [ETM_RAW.nir, ETM_RAW.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat7RAW, self).__init__(**kwargs)


@register
class Landsat7RAWT2(Tier2, Landsat7RAW):
    id = 'LANDSAT/LE07/C02/T2'
    short_name = 'L7RAWT2'

    def __init__(self, **kwargs):
        super(Landsat7RAWT2, self).__init__(**kwargs)


@register
class Landsat8RAW(Landsat, Tier1, OLI_RAW, RAW):
    start_date = '2013-03-18'
    end_date = None
    id = 'LANDSAT/LC08/C02/T1'
    short_name = 'L8RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)    

    bands = (OLI_RAW.ublue, OLI_RAW.blue, OLI_RAW.green, OLI_RAW.red, 
             OLI_RAW.nir, OLI_RAW.swir, OLI_RAW.swir2, OLI_RAW.thermal,
             OLI_RAW.thermal2, OLI_RAW.cirrus, RAW.pan,
             Landsat.qa_pixel, RAW.saa, RAW.sza, RAW.vaa, RAW.vza,
             OLI.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([OLI_RAW.red, OLI_RAW.green, OLI_RAW.blue]),
        FalseColor = Visualization.falseColor([OLI_RAW.nir, OLI_RAW.red, OLI_RAW.green]),
        NSR = Visualization.NSR([OLI_RAW.nir, OLI_RAW.swir, OLI_RAW.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [OLI_RAW.nir, OLI_RAW.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [OLI_RAW.nir, OLI_RAW.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat8RAW, self).__init__(**kwargs)


@register
class Landsat8RAWT2(Tier2, Landsat8RAW):
    id = 'LANDSAT/LC08/C02/T2'
    short_name = 'L8RAWT2'

    def __init__(self, **kwargs):
        super(Landsat8RAWT2, self).__init__(**kwargs)


@register
class Landsat9RAW(Landsat, Tier1, OLI_RAW, RAW):
    start_date = '2021-10-31'
    end_date = None
    id = 'LANDSAT/LC09/C02/T1'
    short_name = 'L9RAW'
    masks = (Mask.fromBand('qa_pixel', Landsat.qa_pixel),)    

    bands = (OLI_RAW.ublue, OLI_RAW.blue, OLI_RAW.green, OLI_RAW.red, 
             OLI_RAW.nir, OLI_RAW.swir, OLI_RAW.swir2, OLI_RAW.thermal,
             OLI_RAW.thermal2, OLI_RAW.cirrus, RAW.pan,
             Landsat.qa_pixel, RAW.saa, RAW.sza, RAW.vaa, RAW.vza,
             OLI.qa_radsat)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([OLI_RAW.red, OLI_RAW.green, OLI_RAW.blue]),
        FalseColor = Visualization.falseColor([OLI_RAW.nir, OLI_RAW.red, OLI_RAW.green]),
        NSR = Visualization.NSR([OLI_RAW.nir, OLI_RAW.swir, OLI_RAW.red])
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [OLI_RAW.nir, OLI_RAW.red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [OLI_RAW.nir, OLI_RAW.swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat9RAW, self).__init__(**kwargs)


@register
class Landsat9RAWT2(Tier2, Landsat9RAW):
    id = 'LANDSAT/LC09/C02/T2'
    short_name = 'L9RAWT2'

    def __init__(self, **kwargs):
        super(Landsat9RAWT2, self).__init__(**kwargs)


def complete_sr(start=None, end=None, site=None, satellites=(9, 8, 7, 5),
                include_slc_off=False, apply_mask='negatives',
                mask_positives=('clear',), mask_negatives=('water', 'shadow',
                'snow', 'cloud', 'cloud_confidence_high',
                'cirrus_confidence_high'), cloud_cover=None,
                landsat_band=True):
    """ Create a collection with all Landsat images. The bands will be always
    renamed. You can set landsat_band to True to keep track where the pixels
    are coming from """
    mask = Mask.fromBand('Landsat', Landsat.qa_pixel)
    if apply_mask:
        apply_mask = 'negatives' if apply_mask is True else apply_mask
        if apply_mask.lower() == 'negatives':
            maskf = partial(mask.apply, negatives=mask_negatives, renamed=False)
        elif apply_mask.lower() == 'positives':
            maskf = partial(mask.apply, positives=mask_positives, renamed=False)
        else:
            maskf = None
    else:
        maskf = None

    satrel = {4: Landsat4SR, 5: Landsat5SR, 7: Landsat7SR, 8: Landsat8SR,
              9: Landsat9SR}
    sats = [satrel[s]() for s in satellites if s in satrel.keys()]
    bands = getCommonBands(*sats, match='alias')
    col = ee.ImageCollection([])
    for lid, s in satrel.items():
        if lid not in satellites:
            continue
        sat = satrel[lid]()
        s = start or sat.start_date
        e = end or sat.end_date or TODAY
        c = sat.collection(site, (s, e))
        if not include_slc_off and sat.short_name in [
            'L7SR', 'L7TOA', 'L7RAW', 'L7SRT2', 'L7TOAT2', 'L7RAWT2']:
            slc = ee.Date(SLC_OFF)
            cond = slc.difference(s, 'day').gt(0)
            c = ee.ImageCollection(
                ee.Algorithms.If(cond, c.filterDate(s, SLC_OFF),
                                 c.filterDate('1970-01-01', SLC_OFF)))
        c = c.filter(ee.Filter.lte(sat.cloud_cover, cloud_cover)) if cloud_cover else c
        c = c.map(lambda i: maskf(i)) if maskf else c
        c = c.map(lambda i: sat.rename(i).select(bands))
        # add landsat id as a band
        if landsat_band:
            c = c.map(lambda i: i.addBands(ee.Image.constant(int(lid)).toUint8().rename('landsat')))
        col = col.merge(c)
    return col
