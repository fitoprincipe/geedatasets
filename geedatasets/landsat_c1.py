# coding=utf-8
""" Google Earth Engine Landsat Collections """
from .visualization import *
from .datasets import OpticalSatellite, ImageCollection
from .bands import OpticalBand, BitBand, ClassificationBand, ExpressionBand,\
    Precisions
from .helpers import TODAY
from functools import partial
from . import register
from .masks import Mask
import geetools


START = {1: '1972-07-23', 2: '1975-01-22', 3: '1978-03-05',
         4: '1982-07-16', 5: '1984-01-01', 7: '1999-01-01',
         8: '2013-04-11'}
END = {1: '1978-01-07', 2: '1982-02-26', 3: '1983-03-31',
       4: '1993-12-14', 5: '2012-05-05', 7: TODAY, 8: TODAY}

IDS = [
    'LANDSAT/LM01/C01/T2',
    'LANDSAT/LM02/C01/T2',
    'LANDSAT/LM03/C01/T2',
    'LANDSAT/LM04/C01/T2',
    'LANDSAT/LM05/C01/T1',
    'LANDSAT/LM05/C01/T2',
    'LANDSAT/LT04/C01/T1', 'LANDSAT/LT04/C01/T1_TOA', 'LANDSAT/LT04/C01/T1_SR',
    'LANDSAT/LT04/C01/T2', 'LANDSAT/LT04/C01/T2_TOA', 'LANDSAT/LT04/C01/T2_SR',
    'LANDSAT/LT05/C01/T1', 'LANDSAT/LT05/C01/T1_TOA', 'LANDSAT/LT05/C01/T1_SR',
    'LANDSAT/LT05/C01/T2', 'LANDSAT/LT05/C01/T2_TOA', 'LANDSAT/LT05/C01/T2_SR',
    'LANDSAT/LE07/C01/T1', 'LANDSAT/LE07/C01/T1_TOA', 'LANDSAT/LE07/C01/T1_SR',
    'LANDSAT/LE07/C01/T2', 'LANDSAT/LE07/C01/T2_TOA', 'LANDSAT/LE07/C01/T2_SR',
    'LANDSAT/LC08/C01/T1', 'LANDSAT/LC08/C01/T1_TOA', 'LANDSAT/LC08/C01/T1_SR',
    'LANDSAT/LC08/C01/T2', 'LANDSAT/LC08/C01/T2_TOA', 'LANDSAT/LC08/C01/T2_SR',
]


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

    def __init__(self, **kwargs):
        super(Landsat, self).__init__(**kwargs)


class Tier1:
    tier = 1


class Tier2:
    tier = 2


class RAW:
    _extra = dict(precision='uint8', scale=1/255)
    process = 'RAW'


class TOA:
    _extra = dict(precision='float')
    process = 'TOA'


def atm_op_decoder(image):
    scale = 0.001
    image = image.multiply(scale)
    clear = image.lt(0.1).rename('clear')
    hazy= image.gt(0.3).rename('hazy')
    average = image.gte(0.1).And(image.lte(0.3)).rename('average')
    return geetools.tools.image.mixBands([clear, hazy, average])


class SR:
    _extra = dict(scale=0.0001, precision=Precisions.int16)
    process = 'SR'

    atm_op = ClassificationBand(
        name='sr_atmos_opacity',
        alias='atmos_opacity',
        precision=Precisions.int16,
        resolution=30,
        classes= dict(
            clear= 'value<0.1',
            average= '0.1>=value<=0.3',
            hazy= 'value>0.3'
        ),
        decoder=atm_op_decoder,
        positives=['clear', 'average'],
        negatives=['hazy']
    )
    atm_op_vis = Visualization('atmos_opacity', [atm_op],
                               0, 300, ['green', 'red'])

    sr_cloud_qa = BitBand(
        name='sr_cloud_qa',
        alias='cloud_qa',
        precision=Precisions.uint8,
        resolution=30,
        bits={
            '0': {1:'ddv'},
            '1': {1:'cloud'},
            '2': {1:'shadow'},
            '3': {1:'adjacent'},
            '4': {1:'snow'},
            '5': {1:'water'}
        },
        negatives=['ddv', 'cloud', 'shadow', 'adjacent', 'snow', 'water']
    )

    pixel_qa = BitBand(
        name='pixel_qa',
        alias='pixel_qa',
        precision=Precisions.uint16,
        resolution=30,
        bits={'1': {1:'clear'}, '2': {1:'water'},
              '3': {1:'shadow'}, '4': {1:'snow'},
              '5': {1:'cloud'},
              '6-7':{3:'high_confidence_cloud'},
              '8-9':{3:'high_confidence_cirrus'}
              },
        positives=['clear'],
        negatives=['water', 'shadow', 'snow', 'cloud', 'high_confidence_cloud',
                   'high_confidence_cirrus']
    )

    radsat_qa = BitBand(
        name='radsat_qa',
        alias='radsat_qa',
        precision=Precisions.uint8,
        resolution=30,
        bits={
            1: {1:'B1_saturated'},
            2: {1:'B2_saturated'},
            3: {1:'B3_saturated'},
            4: {1:'B4_saturated'},
            5: {1:'B5_saturated'},
            6: {1:'B6_saturated'},
            7: {1:'B7_saturated'},
        },
        negatives=['B1_saturated', 'B2_saturated', 'B3_saturated',
                   'B4_saturated', 'B5_saturated', 'B6_saturated',
                   'B7_saturated']
    )


class MSS:
    sensor = 'MSS'
    green = partial(OpticalBand, alias='green', resolution=60,
                    units='DN', wavelength=(0.5, 0.6))
    red = partial(OpticalBand, alias='red', resolution=60,
                  units='DN', wavelength=(0.6, 0.7))
    nir = partial(OpticalBand, alias='nir', resolution=60,
                  units='DN', wavelength=(0.7, 0.8))
    nir2 = partial(OpticalBand, alias='nir2', resolution=30,
                   units='DN', wavelength=(0.8, 1.1))

    bqa = BitBand(
        name='BQA',
        alias='bqa',
        precision=Precisions.uint16,
        resolution=60,
        bits={'4': {1: 'cloud'}},
        negatives=['cloud']
    )
    masks = (Mask.fromBand('BQA', bqa),)


class TM:
    sensor = 'TM'
    blue = partial(OpticalBand, 'B1', 'blue', resolution=30,
                   wavelength=(0.45, 0.52))
    green = partial(OpticalBand, 'B2', 'green', resolution=30,
                    wavelength=(0.52, 0.6))
    red = partial(OpticalBand, 'B3', 'red', resolution=30,
                  wavelength=(0.63, 0.69))
    nir = partial(OpticalBand, 'B4', 'nir', resolution=30,
                  wavelength=(0.76, 0.9))
    swir = partial(OpticalBand, 'B5', 'swir', resolution=30,
                   wavelength=(1.55, 1.75))
    thermal = partial(OpticalBand, 'B6', 'thermal', units='Kelvin',
                      resolution=30, wavelength=(10.4, 12.5))
    swir2 = partial(OpticalBand, 'B7', 'swir2', resolution=30,
                    wavelength=(2.08, 2.35))
    bqa = BitBand(
        name='BQA', alias='bqa',
        precision=Precisions.uint16, resolution=30,
        bits= {
            '4': {1: 'cloud'},
            '5-6': {3: 'high_confidence_cloud'},
            '7-8': {3: 'shadow'},
            '9-10': {3: 'snow'}
        }
    )

    masks = (Mask.fromBand('BQA', bqa),)


class ETM:
    sensor = 'ETM+'
    blue = TM.blue
    green = TM.green
    red = TM.red
    nir = TM.nir
    swir = TM.swir
    thermal = TM.thermal
    thermal_vcid_1 = partial(OpticalBand, 'B6_VCID_1', 'B6_vcid_1',
                             units='Kelvin', resolution=30,
                             wavelength=(10.4, 12.5))
    thermal_vcid_2 = partial(OpticalBand, 'B6_VCID_2', 'B6_vcid_2',
                             units='Kelvin', resolution=30,
                             wavelength=(10.4, 12.5))
    swir2 = TM.swir2
    bqa = TM.bqa
    masks = (Mask.fromBand('BQA', bqa),)


class OLI:
    sensor = 'OLI'
    aerosol = partial(OpticalBand, 'B1', 'coastal_aerosol',
                      resolution=30, wavelength=(0.43, 0.45))
    blue = partial(OpticalBand, 'B2', 'blue',
                   resolution=30, wavelength=(0.45, 0.51))
    green = partial(OpticalBand, 'B3', 'green',
                    resolution=30, wavelength=(0.53, 0.59))
    red = partial(OpticalBand, 'B4', 'red',
                  resolution=30, wavelength=(0.64, 0.67))
    nir = partial(OpticalBand, 'B5', 'nir',
                  resolution=30, wavelength=(0.85, 0.88))
    swir = partial(OpticalBand, 'B6', 'swir',
                   resolution=30, wavelength=(1.57, 1.65))
    swir2 = partial(OpticalBand, 'B7', 'swir2',
                    resolution=30, wavelength=(2.11, 2.29))
    pan = partial(OpticalBand, 'B8', 'pan',
                  resolution=15, wavelength=(0.52, 0.9))
    cirrus = partial(OpticalBand, 'B9', 'cirrus',
                     resolution=15, wavelength=(1.36, 1.38))
    thermal = partial(OpticalBand, 'B10', 'thermal',
                      resolution=30, scale=0.1, wavelength=(10.60, 11.19))
    thermal2 = partial(OpticalBand, 'B11', 'thermal2',
                       resolution=30, scale=0.1, wavelength=(11.50, 12.51))
    bqa = BitBand(
        name='BQA', alias='bqa',
        precision=Precisions.uint16, resolution=30,
        bits= {
            '4': {1: 'cloud'},
            '5-6': {3: 'high_confidence_cloud'},
            '7-8': {3: 'shadow'},
            '9-10': {3: 'snow'},
            '11-12': {3: 'cirrus'}
        }
    )

    class RAW:
        _extra = dict(precision=Precisions.uint16)

    class SR:
        aerosol = BitBand('sr_aerosol', 'sr_aerosol',
                          precision=Precisions.uint8,
                          bits= {
                              '1': {1: 'aerosol_valid'},
                              '2': {1: 'aerosol_interpolated'},
                              '3': {1: 'water'},
                              '6-7': {0: 'climatology', 1: 'low',
                                      2:'medium', 3:'high'}
                          })

        radsat_qa = BitBand(
            name='radsat_qa',
            alias='radsat_qa',
            precision=Precisions.uint16,
            resolution=30,
            bits={
                1: {1:'B1_saturated'},
                2: {1:'B2_saturated'},
                3: {1:'B3_saturated'},
                4: {1:'B4_saturated'},
                5: {1:'B5_saturated'},
                6: {1:'B6_saturated'},
                7: {1:'B7_saturated'},
                9: {1:'B9_saturated'},
                10: {1:'B10_saturated'},
                11: {1:'B11_saturated'}
            },
            negatives=['B1_saturated', 'B2_saturated', 'B3_saturated',
                       'B4_saturated', 'B5_saturated', 'B6_saturated',
                       'B7_saturated', 'B9_saturated', 'B10_saturated',
                       'B11_saturated']
        )


class Landsat1(MSS, Landsat):
    """ Landsat 1 """
    number = 1
    start_date = '1972-07-26'
    end_date = '1978-01-06'


@register
class Landsat1RAW_C1(Tier1, RAW, Landsat1):
    """ Landsat 1 Tier 1 """
    id = 'LANDSAT/LM01/C01/T1'
    short_name = 'L1RAW_C1'
    green = MSS.green(name='B4', **RAW._extra)
    red = MSS.red(name='B5', **RAW._extra)
    nir = MSS.nir(name='B6', **RAW._extra)
    nir2 = MSS.nir2(name='B7', **RAW._extra)
    bands = (green, red, nir, nir2, MSS.bqa)
    visualizers = Visualizers(
        FalseColor = Visualization.falseColor([nir, red, green]),
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')
    extra_bands = (ndvi,)

    def __init__(self, **kwargs):
        super(Landsat1RAW_C1, self).__init__(**kwargs)


@register
class Landsat1RAWT2_C1(Tier2, Landsat1RAW_C1):
    """ Landsat 1 Tier 2 """
    id = 'LANDSAT/LM01/C01/T2'
    short_name = 'L1RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat1RAWT2_C1, self).__init__(**kwargs)


class Landsat2(MSS, Landsat):
    number = 2
    start_date = '1975-01-31'
    end_date = '1982-02-03'


@register
class Landsat2RAW_C1(Tier1, RAW, Landsat2):
    """ Landsat 2 Tier 1 """
    id = 'LANDSAT/LM02/C01/T1'
    short_name = 'L2RAW_C1'
    bands = (Landsat1RAW_C1.green, Landsat1RAW_C1.red, Landsat1RAW_C1.nir,
             Landsat1RAW_C1.nir2, MSS.bqa)
    visualizers = Landsat1RAW_C1.visualizers
    extra_bands = Landsat1RAW_C1.extra_bands

    def __init__(self, **kwargs):
        super(Landsat2RAW_C1, self).__init__(**kwargs)


@register
class Landsat2RAWT2_C1(Tier2, Landsat2RAW_C1):
    id = 'LANDSAT/LM02/C01/T2'
    short_name = 'L2RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat2RAWT2_C1, self).__init__(**kwargs)


class Landsat3(MSS, Landsat):
    number = 3
    start_date = '1978-06-03'
    end_date = '1983-02-23'


@register
class Landsat3RAW_C1(Tier1, RAW, Landsat3):
    id = 'LANDSAT/LM03/C01/T1'
    short_name = 'L3RAW_C1'
    bands = (Landsat1RAW_C1.green, Landsat1RAW_C1.red, Landsat1RAW_C1.nir, Landsat1RAW_C1.nir2,
             MSS.bqa)
    visualizers = Landsat1RAW_C1.visualizers
    extra_bands = Landsat1RAW_C1.extra_bands

    def __init__(self, **kwargs):
        super(Landsat3RAW_C1, self).__init__(**kwargs)


@register
class Landsat3RAWT2_C1(Tier2, Landsat3RAW_C1):
    id = 'LANDSAT/LM03/C01/T2'
    short_name = 'L3RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat3RAWT2_C1, self).__init__(**kwargs)


class Landsat4MSS(MSS, Landsat):
    number = 4
    start_date = '1982-08-14'
    end_date = '1992-08-28'


@register
class Landsat4MSSRAW_C1(Tier1, RAW, Landsat4MSS):
    """ Landsat 4 MSS """
    id = 'LANDSAT/LM04/C01/T1'
    short_name = 'L4MSSRAW_C1'
    green = MSS.green(name='B1', **RAW._extra)
    red = MSS.red(name='B2', **RAW._extra)
    nir = MSS.nir(name='B3', **RAW._extra)
    nir2 = MSS.nir2(name='B4', **RAW._extra)
    bands = (green, red, nir, nir2, MSS.bqa)
    visualizers = Visualizers(
        FalseColor = Visualization.falseColor([nir, red, green]),
    )
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')
    extra_bands = (ndvi,)

    def __init__(self, **kwargs):
        super(Landsat4MSSRAW_C1, self).__init__(**kwargs)


@register
class Landsat4MSSRAWT2_C1(Tier2, Landsat4MSSRAW_C1):
    """ Landsat 4 MSS """
    id = 'LANDSAT/LM04/C01/T2'
    short_name = 'L4MSSRAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat4MSSRAWT2_C1, self).__init__(**kwargs)


class Landsat4TM(TM, Landsat):
    """ Landsat 4 TM Raw Tier 1 """
    number = 4
    start_date = '1982-08-22'
    end_date = '1993-11-18'


@register
class Landsat4RAW_C1(Tier1, RAW, Landsat4TM):
    id = 'LANDSAT/LT04/C01/T1'
    short_name = 'L4RAW_C1'
    blue = TM.blue(**RAW._extra)
    green = TM.green(**RAW._extra)
    red = TM.red(**RAW._extra)
    nir = TM.nir(**RAW._extra)
    swir = TM.swir(**RAW._extra)
    thermal = TM.thermal(**RAW._extra)
    swir2 = TM.swir2(**RAW._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision='float')
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat4RAW_C1, self).__init__(**kwargs)


@register
class Landsat4RAWT2_C1(Tier2, Landsat4RAW_C1):
    id = 'LANDSAT/LT04/C01/T2'
    short_name = 'L4RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat4RAWT2_C1, self).__init__(**kwargs)


@register
class Landsat4TOA_C1(Tier1, TOA, Landsat4TM):
    id = 'LANDSAT/LT04/C01/T1_TOA'
    short_name = 'L4TOA_C1'
    blue = TM.blue(**TOA._extra)
    green = TM.green(**TOA._extra)
    red = TM.red(**TOA._extra)
    nir = TM.nir(**TOA._extra)
    swir = TM.swir(**TOA._extra)
    thermal = TM.thermal(**TOA._extra)
    swir2 = TM.swir2(**TOA._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision='float')
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat4TOA_C1, self).__init__(**kwargs)


@register
class Landsat4TOAT2_C1(Tier2, Landsat4TOA_C1):
    id = 'LANDSAT/LT04/C01/T2_TOA'
    short_name = 'L4TOAT2_C1'

    def __init__(self, **kwargs):
        super(Landsat4TOAT2_C1, self).__init__(**kwargs)


@register
class Landsat4SR_C1(Tier1, SR, Landsat4TM):
    id = 'LANDSAT/LT04/C01/T1_SR'
    short_name = 'L4SR_C1'
    blue = TM.blue(**SR._extra)
    green = TM.green(**SR._extra)
    red = TM.red(**SR._extra)
    nir = TM.nir(**SR._extra)
    swir = TM.swir(**SR._extra)
    thermal = TM.thermal(**SR._extra)
    swir2 = TM.swir2(**SR._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, SR.atm_op,
             SR.sr_cloud_qa, SR.pixel_qa, SR.radsat_qa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red]),
        AtmosphericOpacity = SR.atm_op_vis
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),
             Mask.fromBand('cloud_qa', SR.sr_cloud_qa))
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)

    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat4SR_C1, self).__init__(**kwargs)


@register
class Landsat4SRT2_C1(Tier2, Landsat4SR_C1):
    id = 'LANDSAT/LT04/C01/T2_SR'
    short_name = 'L4SRT2_C1'

    def __init__(self, **kwargs):
        super(Landsat4SRT2_C1, self).__init__(**kwargs)


class Landsat5MSS(MSS, Landsat):
    number = 5
    start_date = '1984-04-09'
    end_date = '2013-01-01'


@register
class Landsat5MSSRAW_C1(Tier1, RAW, Landsat4MSS):
    """ Landsat 5 MSS """
    id = 'LANDSAT/LM05/C01/T1'
    short_name = 'L5MSSRAW_C1'
    green = MSS.green(name='B1', **RAW._extra)
    red = MSS.red(name='B2', **RAW._extra)
    nir = MSS.nir(name='B3', **RAW._extra)
    nir2 = MSS.nir2(name='B4', **RAW._extra)
    bands = (green, red, nir, nir2, MSS.bqa)
    visualizers = Visualizers(
        FalseColor = Visualization.falseColor([nir, red, green]),
    )
    masks = (Mask.fromBand('BQA', MSS.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision='float')
    extra_bands = (ndvi,)

    def __init__(self, **kwargs):
        super(Landsat5MSSRAW_C1, self).__init__(**kwargs)


@register
class Landsat5MSSRAWT2_C1(Tier2, Landsat5MSSRAW_C1):
    """ Landsat 5 MSS """
    id = 'LANDSAT/LM05/C01/T2'
    short_name = 'L5MSSRAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat5MSSRAWT2_C1, self).__init__(**kwargs)


class Landsat5TM(TM, Landsat):
    """ Landsat 5 TM Raw Tier 1 """
    number = 5
    start_date = '1984-03-16'
    end_date = '2012-05-05'


@register
class Landsat5RAW_C1(Tier1, RAW, Landsat5TM):
    id = 'LANDSAT/LT05/C01/T1'
    short_name = 'L5RAW_C1'
    blue = TM.blue(**RAW._extra)
    green = TM.green(**RAW._extra)
    red = TM.red(**RAW._extra)
    nir = TM.nir(**RAW._extra)
    swir = TM.swir(**RAW._extra)
    thermal = TM.thermal(**RAW._extra)
    swir2 = TM.swir2(**RAW._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat5RAW_C1, self).__init__(**kwargs)


@register
class Landsat5RAWT2_C1(Tier2, Landsat5RAW_C1):
    id = 'LANDSAT/LT05/C01/T2'
    short_name = 'L5RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat5RAWT2_C1, self).__init__(**kwargs)


@register
class Landsat5TOA_C1(Tier1, TOA, Landsat5TM):
    id = 'LANDSAT/LT05/C01/T1_TOA'
    short_name = 'L5TOA_C1'
    blue = TM.blue(**TOA._extra)
    green = TM.green(**TOA._extra)
    red = TM.red(**TOA._extra)
    nir = TM.nir(**TOA._extra)
    swir = TM.swir(**TOA._extra)
    thermal = TM.thermal(**TOA._extra)
    swir2 = TM.swir2(**TOA._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat5TOA_C1, self).__init__(**kwargs)


@register
class Landsat5TOAT2_C1(Tier2, Landsat5TOA_C1):
    id = 'LANDSAT/LT05/C01/T2_TOA'
    short_name = 'L5TOAT2_C1'

    def __init__(self, **kwargs):
        super(Landsat5TOAT2_C1, self).__init__(**kwargs)


@register
class Landsat5SR_C1(Tier1, SR, Landsat5TM):
    id = 'LANDSAT/LT05/C01/T1_SR'
    short_name = 'L5SR_C1'
    blue = TM.blue(**SR._extra)
    green = TM.green(**SR._extra)
    red = TM.red(**SR._extra)
    nir = TM.nir(**SR._extra)
    swir = TM.swir(**SR._extra)
    thermal = TM.thermal(**SR._extra)
    swir2 = TM.swir2(**SR._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, SR.atm_op,
             SR.sr_cloud_qa, SR.pixel_qa, SR.radsat_qa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red]),
        AtmosphericOpacity = SR.atm_op_vis
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),
             Mask.fromBand('cloud_qa', SR.sr_cloud_qa))
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat5SR_C1, self).__init__(**kwargs)


@register
class Landsat5SRT2_C1(Tier2, Landsat5SR_C1):
    id = 'LANDSAT/LT05/C01/T2_SR'
    short_name = 'L5SRT2_C1'

    def __init__(self, **kwargs):
        super(Landsat5SRT2_C1, self).__init__(**kwargs)


class Landsat7ETM(ETM, Landsat):
    """ Landsat 7 ETM Raw Tier 1 """
    number = 7
    start_date = '1999-05-28'
    end_date = None
    scl_off = '2003-05-31'


@register
class Landsat7RAW_C1(Tier1, RAW, Landsat7ETM):
    id = 'LANDSAT/LE07/C01/T1'
    short_name = 'L7RAW_C1'
    blue = ETM.blue(**RAW._extra)
    green = ETM.green(**RAW._extra)
    red = ETM.red(**RAW._extra)
    nir = ETM.nir(**RAW._extra)
    swir = ETM.swir(**RAW._extra)    
    swir2 = ETM.swir2(**RAW._extra)
    thermal1 = ETM.thermal_vcid_1(**RAW._extra)
    thermal2 = ETM.thermal_vcid_2(**RAW._extra)
    bands = (blue, green, red, nir, swir, thermal1, thermal2, swir2, TM.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat7RAW_C1, self).__init__(**kwargs)


@register
class Landsat7RAWT2_C1(Tier2, Landsat7RAW_C1):
    id = 'LANDSAT/LE07/C01/T2'
    short_name = 'L7RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat7RAWT2_C1, self).__init__(**kwargs)


@register
class Landsat7TOA_C1(Tier1, TOA, Landsat7ETM):
    id = 'LANDSAT/LE07/C01/T1_TOA'
    short_name = 'L7TOA_C1'
    blue = ETM.blue(**TOA._extra)
    green = ETM.green(**TOA._extra)
    red = ETM.red(**TOA._extra)
    nir = ETM.nir(**TOA._extra)
    swir = ETM.swir(**TOA._extra)
    swir2 = ETM.swir2(**TOA._extra)
    thermal1 = ETM.thermal_vcid_1(**TOA._extra)
    thermal2 = ETM.thermal_vcid_2(**TOA._extra)
    bands = (blue, green, red, nir, swir, thermal1, thermal2, swir2, TM.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat7TOA_C1, self).__init__(**kwargs)


@register
class Landsat7TOAT2_C1(Tier2, Landsat7TOA_C1):
    id = 'LANDSAT/LE07/C01/T2_TOA'
    short_name = 'L7TOAT2_C1'

    def __init__(self, **kwargs):
        super(Landsat7TOAT2_C1, self).__init__(**kwargs)


@register
class Landsat7SR_C1(Tier1, SR, Landsat7ETM):
    id = 'LANDSAT/LE07/C01/T1_SR'
    short_name = 'L7SR_C1'
    blue = ETM.blue(**SR._extra)
    green = ETM.green(**SR._extra)
    red = ETM.red(**SR._extra)
    nir = ETM.nir(**SR._extra)
    swir = ETM.swir(**SR._extra)
    thermal = ETM.thermal(**SR._extra)
    swir2 = ETM.swir2(**SR._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, SR.atm_op,
             SR.sr_cloud_qa, SR.pixel_qa, SR.radsat_qa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red]),
        AtmosphericOpacity = SR.atm_op_vis
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),
             Mask.fromBand('cloud_qa', SR.sr_cloud_qa))
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat7SR_C1, self).__init__(**kwargs)


@register
class Landsat7SRT2_C1(Tier2, Landsat7SR_C1):
    id = 'LANDSAT/LE07/C01/T2_SR'
    short_name = 'L7SRT2_C1'

    def __init__(self, **kwargs):
        super(Landsat7SRT2_C1, self).__init__(**kwargs)


class Landsat8OLI(OLI, Landsat):
    """ Landsat 8 OLI Raw Tier 1 """
    number = 8
    start_date = '2013-03-18'
    end_date = None


@register
class Landsat8RAW_C1(Tier1, RAW, Landsat8OLI):
    id = 'LANDSAT/LC08/C01/T1'
    short_name = 'L8RAW_C1'
    aerosol = OLI.aerosol(**OLI.RAW._extra)
    blue = OLI.blue(**OLI.RAW._extra)
    green = OLI.green(**OLI.RAW._extra)
    red = OLI.red(**OLI.RAW._extra)
    nir = OLI.nir(**OLI.RAW._extra)
    swir = OLI.swir(**OLI.RAW._extra)
    swir2 = OLI.swir2(**OLI.RAW._extra)
    pan = OLI.pan(**OLI.RAW._extra)
    cirrus = OLI.cirrus(**OLI.RAW._extra)
    thermal1 = OLI.thermal(**OLI.RAW._extra)
    thermal2 = OLI.thermal2(**OLI.RAW._extra)
    bands = (aerosol, blue, green, red, nir, swir, swir2, pan, cirrus,
             thermal1, thermal2, OLI.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', OLI.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat8RAW_C1, self).__init__(**kwargs)


@register
class Landsat8RAWT2_C1(Tier2, Landsat8RAW_C1):
    id = 'LANDSAT/LC08/C01/T2'
    short_name = 'L8RAWT2_C1'

    def __init__(self, **kwargs):
        super(Landsat8RAWT2_C1, self).__init__(**kwargs)


@register
class Landsat8TOA_C1(Tier1, TOA, Landsat8OLI):
    id = 'LANDSAT/LC08/C01/T1_TOA'
    short_name = 'L8TOA_C1'
    aerosol = OLI.aerosol(**TOA._extra)
    blue = OLI.blue(**TOA._extra)
    green = OLI.green(**TOA._extra)
    red = OLI.red(**TOA._extra)
    nir = OLI.nir(**TOA._extra)
    swir = OLI.swir(**TOA._extra)
    swir2 = OLI.swir2(**TOA._extra)
    pan = OLI.pan(**TOA._extra)
    cirrus = OLI.cirrus(**TOA._extra)
    thermal1 = OLI.thermal(**TOA._extra)
    thermal2 = OLI.thermal2(**TOA._extra)
    bands = (aerosol, blue, green, red, nir, swir, swir2, pan, cirrus,
             thermal1, thermal2, OLI.bqa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', OLI.bqa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat8TOA_C1, self).__init__(**kwargs)


@register
class Landsat8TOAT2_C1(Tier2, Landsat8TOA_C1):
    id = 'LANDSAT/LC08/C01/T2_TOA'
    short_name = 'L8TOAT2_C1'

    def __init__(self, **kwargs):
        super(Landsat8TOAT2_C1, self).__init__(**kwargs)


@register
class Landsat8SR_C1(Tier1, SR, Landsat8OLI):
    id = 'LANDSAT/LC08/C01/T1_SR'
    short_name = 'L8SR_C1'
    aerosol = OLI.aerosol(**SR._extra)
    blue = OLI.blue(**SR._extra)
    green = OLI.green(**SR._extra)
    red = OLI.red(**SR._extra)
    nir = OLI.nir(**SR._extra)
    swir = OLI.swir(**SR._extra)
    swir2 = OLI.swir2(**SR._extra)
    thermal1 = OLI.thermal(**SR._extra)
    thermal2 = OLI.thermal2(**SR._extra)
    bands = (aerosol, blue, green, red, nir, swir, swir2, thermal1, thermal2,
             OLI.SR.aerosol, SR.pixel_qa, OLI.SR.radsat_qa)
    visualizers = Visualizers(
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        NSR = Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),)
    ndvi = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                          [nir, red], precision=Precisions.float)
    nbr = ExpressionBand('NBR', 'nbr', '(nir-swir)/(nir+swir)',
                         [nir, swir], precision=Precisions.float)
    extra_bands = (ndvi, nbr)

    def __init__(self, **kwargs):
        super(Landsat8SR_C1, self).__init__(**kwargs)


@register
class Landsat8SRT2_C1(Tier2, Landsat8SR_C1):
    id = 'LANDSAT/LC08/C01/T2_SR'
    short_name = 'L8SRT2_C1'

    def __init__(self, **kwargs):
        super(Landsat8SRT2_C1, self).__init__(**kwargs)