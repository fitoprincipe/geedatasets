# coding=utf-8
""" Google Earth Engine Sentinel Collections """
from .visualization import *
from .datasets import OpticalSatellite, ImageCollection
from .bands import OpticalBand, BitBand, ClassificationBand, ExpressionBand
from .masks import Mask
import geetools
from . import register


class Sentinel2(OpticalSatellite, ImageCollection):
    """ Sentinel 2 Base Class - Shared between TOA and SR """
    INFO = """ID: {id}
short name: {short_name}
SpaceCraft: Sentinel 2
Process: {process}
Start Date: {start_date}
End Date: {end_date}
Cloud Cover: {cloud_cover}
Bands: {bands}
Masks: {masks}
Visualizers: {visualizers}
"""
    number = 2
    spacecraft = 'SENTINEL2'
    start_date = '2015-06-23'
    end_date = None
    cloud_cover = 'CLOUD_COVERAGE_ASSESSMENT'
    _common = {'scale':0.0001, 'precision': 'uint16'}
    aerosol = OpticalBand('B1', 'aerosol', resolution=60, **_common)
    blue = OpticalBand('B2', 'blue', resolution=10, **_common)
    green = OpticalBand('B3', 'green', resolution=10, **_common)
    red = OpticalBand('B4', 'red', resolution=10, **_common)
    red_edge_1 = OpticalBand('B5', 'red_edge_1', resolution=20, **_common)
    red_edge_2 = OpticalBand('B6', 'red_edge_2', resolution=20, **_common)
    red_edge_3 = OpticalBand('B7', 'red_edge_3', resolution=20, **_common)
    nir = OpticalBand('B8', 'nir', resolution=10, **_common)
    red_edge_4 = OpticalBand('B8A', 'red_edge_4', resolution=20, **_common)
    water_vapor = OpticalBand('B9', 'water_vapor', resolution=60, **_common)
    swir = OpticalBand('B11', 'swir', resolution=20, **_common)
    swir2 = OpticalBand('B12', 'swir2', resolution=20, **_common)
    qa10 = BitBand(name='QA10', alias='qa10', bits=dict(), resolution=10,
                   precision='uint16')
    qa20 = BitBand(name='QA20', alias='qa20', bits=dict(), resolution=20,
                   precision='uint32')
    qa60 = BitBand(
        'QA60', 'qa60',
        precision='uint16',
        resolution=60,
        bits={'10':{1:'cloud'}, '11':{1:'cirrus_cloud'}},
        negatives=['cloud', 'cirrus_cloud']
    )

    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision='float')

    ndre = ExpressionBand('NDRE', 'ndre', '(nir-red_edge_1)/(nir+red_edge_1)',
                          [nir, red_edge_1], precision='float')

    algorithms = {
        'hollstein': geetools.cloud_mask.applyHollstein
    }

    visualizers = (
        Visualization.NSR([nir, swir, red]),
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green])
    )

    def __init__(self, **kwargs):
        super(Sentinel2, self).__init__(**kwargs)


def decode_hollstein(image, renamed=None):
    """ Decode Hollstein classification into an image """
    classes = ('cloud', 'snow', 'shadow', 'water', 'cirrus')
    aerosol = Sentinel2.aerosol.alias if renamed else Sentinel2.aerosol.name
    blue = Sentinel2.blue.alias if renamed else Sentinel2.blue.name
    green = Sentinel2.green.alias if renamed else Sentinel2.green.name
    red_edge_1 = Sentinel2.red_edge_1.alias if renamed else Sentinel2.red_edge_1.name
    red_edge_2 = Sentinel2.red_edge_2.alias if renamed else Sentinel2.red_edge_2.name
    red_edge_3 = Sentinel2.red_edge_3.alias if renamed else Sentinel2.red_edge_3.name
    red_edge_4 = Sentinel2.red_edge_4.alias if renamed else Sentinel2.red_edge_4.name
    water_vapor = Sentinel2.water_vapor.alias if renamed else Sentinel2.water_vapor.name
    cirrus = Sentinel2TOA.cirrus.alias if renamed else Sentinel2TOA.cirrus.name
    swir = Sentinel2.swir.alias if renamed else Sentinel2.swir.name
    decoded = geetools.cloud_mask.hollsteinMask(
        image, classes, aerosol, blue, green, red_edge_1, red_edge_2,
        red_edge_3, red_edge_4, water_vapor, cirrus, swir
    )
    return decoded


@register
class Sentinel2TOA(Sentinel2):
    """ Sentinel 2 TOA """
    id = 'COPERNICUS/S2'
    short_name = 'S2TOA'
    process = 'TOA'
    cirrus = OpticalBand('B10', 'cirrus', precision='uint16', resolution=60)
    bands = (
        Sentinel2.aerosol, Sentinel2.blue, Sentinel2.green, Sentinel2.red, Sentinel2.red_edge_1,
        Sentinel2.red_edge_2, Sentinel2.red_edge_3, Sentinel2.nir, Sentinel2.red_edge_4, Sentinel2.water_vapor,
        cirrus, Sentinel2.swir, Sentinel2.swir2, Sentinel2.qa10, Sentinel2.qa20,
        Sentinel2.qa60
    )

    masks = (Mask.fromBand('QA60', Sentinel2.qa60),
             Mask('Hollstein',
                  negatives= ('cloud', 'snow', 'shadow', 'water', 'cirrus'),
                  decoder=decode_hollstein))

    def __init__(self, **kwargs):
        super(Sentinel2TOA, self).__init__(**kwargs)


def scl_decoder(image):
    saturated = image.eq(1).rename('saturated')
    dark = image.eq(2).rename('dark')
    shadow = image.eq(3).rename('shadow')
    vegetation = image.eq(4).rename('vegetation')
    bare_soil = image.eq(5).rename('bare_soil')
    water = image.eq(6).rename('water')
    unclassified = image.eq(7).rename('unclassified')
    clouds_mid_probability = image.eq(8).rename('clouds_mid_probability')
    clouds_high_probability = image.eq(9).rename('clouds_high_probability')
    cirrus = image.eq(10).rename('cirrus')
    snow = image.eq(11).rename('snow')
    cloud = clouds_mid_probability.Or(clouds_high_probability).rename('cloud')
    return geetools.tools.image.mixBands([
        saturated, dark, shadow, vegetation, bare_soil, water, unclassified,
        clouds_mid_probability, clouds_high_probability, cirrus, snow, cloud
    ])


@register
class Sentinel2SR(Sentinel2):
    """ Sentinel 2 Surface Reflectance """
    id = 'COPERNICUS/S2_SR'
    short_name = 'S2SR'
    start_date = '2017-03-28'
    process = 'SR'
    aot = OpticalBand('AOT', 'aerosol_thickness', precision='uint16',
                      resolution=10)
    wvp = OpticalBand('WVP', 'water_vapor_pressure', precision='uint16',
                      resolution=10)

    scl = ClassificationBand(
        name='SCL',
        alias='scene_classification_map',
        precision='uint8',
        resolution=20,
        classes= dict(
            saturated=1,
            dark=2,
            shadow=3,
            vegetation=4,
            bare_soil=5,
            water=6,
            unclassified=7,
            clouds_mid_probability=8,
            clouds_high_probability=9,
            cirrus=10,
            snow=11
        ),
        decoder=scl_decoder,
        negatives=['saturated', 'dark', 'shadow', 'water', 'clouds_mid_probability',
                   'clouds_high_probability', 'cirrus', 'snow']
    )
    tci_r = OpticalBand(name='TCI_R', alias='true_color_red',
                        precision='uint8', resolution=10)
    tci_g = OpticalBand(name='TCI_G', alias='true_color_blue',
                        precision='uint8', resolution=10)
    tci_b = OpticalBand(name='TCI_B', alias='true_color_green',
                        precision='uint8', resolution=10)

    # Add visualizer
    visualizers = (
        Visualization.NSR([Sentinel2.nir, Sentinel2.swir, Sentinel2.red]),
        Visualization.trueColor([Sentinel2.red, Sentinel2.green, Sentinel2.blue]),
        Visualization.falseColor([Sentinel2.nir, Sentinel2.red, Sentinel2.green]),
        Visualization('SCL', [scl], 1, 11, palette=[
            '#ff0004', '#868686', '#774b0a', '#10d22c', '#ffff52', '#0000ff',
            '#818181', '#c0c0c0', '#f1f1f1', '#bac5eb', '#52fff9'
        ])
    )

    bands = (
        Sentinel2.aerosol, Sentinel2.blue, Sentinel2.green, Sentinel2.red,
        Sentinel2.red_edge_1, Sentinel2.red_edge_2, Sentinel2.red_edge_3,
        Sentinel2.nir, Sentinel2.red_edge_4, Sentinel2.water_vapor,
        Sentinel2.swir, Sentinel2.swir2, aot, wvp, scl, tci_r, tci_g, tci_b,
        Sentinel2.qa10, Sentinel2.qa20, Sentinel2.qa60
    )

    masks = (Mask.fromBand('SCL', scl), Mask.fromBand('QA60', Sentinel2.qa60))

    def __init__(self, **kwargs):
        super(Sentinel2SR, self).__init__(**kwargs)
