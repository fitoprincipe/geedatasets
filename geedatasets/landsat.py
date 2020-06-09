# coding=utf-8
""" Google Earth Engine Landsat Collections """
from .visualization import *
from .datasets import OpticalSatellite, ImageCollection
from .bands import OpticalBand, BitBand, ClassificationBand
from .helpers import TODAY
from functools import partial
from . import register
from .masks import Mask


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
    _extra = dict(precision='uint8')
    process = 'RAW'


class TOA:
    _extra = dict(precision='float')
    process = 'TOA'


class SR:
    _extra = dict(scale=0.0001, precision='int16')
    process = 'SR'

    atm_op = ClassificationBand(
        name='sr_atmos_opacity',
        alias='atmos_opacity',
        precision='int16',
        resolution=30,
        classes= dict(
            clear='value<0.1',
            average='value>=0.1&&value<=0.3',
            hazy='value>0.3'
        ),
        positives=['clear', 'average'],
        negatives=['hazy']
    )
    sr_cloud_qa = BitBand(
        name='sr_cloud_qa',
        alias='cloud_qa',
        precision='uint8',
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
        precision='uint16',
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
        precision='uint8',
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
    green = partial(OpticalBand, alias='green', resolution=60, units='DN', wavelength=(0.5, 0.6))
    red = partial(OpticalBand, alias='red', resolution=60, units='DN', wavelength=(0.6, 0.7))
    nir = partial(OpticalBand, alias='nir', resolution=60, units='DN', wavelength=(0.7, 0.8))
    nir2 = partial(OpticalBand, alias='nir2', resolution=30, units='DN', wavelength=(0.8, 1.1))
    bqa = BitBand(
        name='BQA',
        alias='bqa',
        precision='uint16',
        resolution=60,
        bits={'4': {1: 'cloud'}},
        negatives=['cloud']
    )
    common_masks = (Mask.fromBand('BQA', bqa),)


class TM:
    sensor = 'TM'
    blue = partial(OpticalBand, 'B1', 'blue', resolution=30, wavelength=(0.45, 0.52))
    green = partial(OpticalBand, 'B2', 'green', resolution=30, wavelength=(0.52, 0.6))
    red = partial(OpticalBand, 'B3', 'red', resolution=30, wavelength=(0.63, 0.69))
    nir = partial(OpticalBand, 'B4', 'nir', resolution=30, wavelength=(0.76, 0.9))
    swir = partial(OpticalBand, 'B5', 'swir', resolution=30, wavelength=(1.55, 1.75))
    thermal = partial(OpticalBand, 'B6', 'thermal', units='Kelvin', resolution=30, wavelength=(10.4, 12.5))
    swir2 = partial(OpticalBand, 'B7', 'swir2', resolution=30, wavelength=(2.08, 2.35))
    bqa = BitBand(
        name='BQA', alias='bqa',
        precision='uint16', resolution=30,
        bits= {
            '4': {1: 'cloud'},
            '5-6': {3: 'high_confidence_cloud'},
            '7-8': {3: 'shadow'},
            '9-10': {3: 'snow'}
        }
    )
    common_masks = (Mask.fromBand('BQA', bqa),)


class ETM:
    sensor = 'ETM+'
    blue = TM.blue
    green = TM.green
    red = TM.red
    nir = TM.nir
    swir = TM.swir
    thermal = TM.thermal
    thermal_vcid_1 = partial(OpticalBand, 'B6_VCID_1', 'B6_vcid_1', units='Kelvin',
                                 resolution=30, wavelength=(10.4, 12.5))
    thermal_vcid_2 = partial(OpticalBand, 'B6_VCID_2', 'B6_vcid_2', units='Kelvin',
                                 resolution=30, wavelength=(10.4, 12.5))
    swir2 = TM.swir2
    bqa = TM.bqa
    common_masks = (Mask.fromBand('BQA', bqa),)


class OLI:
    sensor = 'OLI'
    aerosol = partial(OpticalBand, 'B1', 'coastal_aerosol', resolution=30, wavelength=(0.43, 0.45))
    blue = partial(OpticalBand, 'B2', 'blue', resolution=30, wavelength=(0.45, 0.51))
    green = partial(OpticalBand, 'B3', 'green', resolution=30, wavelength=(0.53, 0.59))
    red = partial(OpticalBand, 'B4', 'red', resolution=30, wavelength=(0.64, 0.67))
    nir = partial(OpticalBand, 'B5', 'nir', resolution=30, wavelength=(0.85, 0.88))
    swir = partial(OpticalBand, 'B6', 'swir', resolution=30, wavelength=(1.57, 1.65))
    swir2 = partial(OpticalBand, 'B7', 'swir2', resolution=30, wavelength=(2.11, 2.29))
    pan = partial(OpticalBand, 'B8', 'pan', resolution=15, wavelength=(0.52, 0.9))
    cirrus = partial(OpticalBand, 'B9', 'cirrus', resolution=15, wavelength=(1.36, 1.38))
    thermal = partial(OpticalBand, 'B10', 'thermal', resolution=30, scale=0.1, wavelength=(10.60, 11.19))
    thermal2 = partial(OpticalBand, 'B11', 'thermal2', resolution=30, scale=0.1, wavelength=(11.50, 12.51))
    bqa = BitBand(
        name='BQA', alias='bqa',
        precision='uint16', resolution=30,
        bits= {
            '4': {1: 'cloud'},
            '5-6': {3: 'high_confidence_cloud'},
            '7-8': {3: 'shadow'},
            '9-10': {3: 'snow'},
            '11-12': {3: 'cirrus'}
        }
    )

    class RAW:
        _extra = dict(precision='uint16')

    class SR:
        aerosol = BitBand('sr_aerosol', 'sr_aerosol',
                          precision='uint8',
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
            precision='uint16',
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
    start_date = '1975-01-22'
    end_date = '1978-01-07'


@register
class Landsat1RAW(Tier1, RAW, Landsat1):
    """ Landsat 1 Tier 1 """
    id = 'LANDSAT/LM01/C01/T1'
    short_name = 'L1RAW'
    green = MSS.green(name='B4', **RAW._extra)
    red = MSS.red(name='B5', **RAW._extra)
    nir = MSS.nir(name='B6', **RAW._extra)
    nir2 = MSS.nir2(name='B7', **RAW._extra)
    bands = (green, red, nir, nir2, MSS.bqa)
    visualizers = (
        Visualization.falseColor([nir, red, green]),
    )

    def __init__(self, **kwargs):
        super(Landsat1RAW, self).__init__(**kwargs)


@register
class Landsat1RAWT2(Tier2, Landsat1RAW):
    """ Landsat 1 Tier 2 """
    id = 'LANDSAT/LM01/C01/T2'
    short_name = 'L1RAWT2'

    def __init__(self, **kwargs):
        super(Landsat1RAWT2, self).__init__(**kwargs)


class Landsat2(MSS, Landsat):
    number = 2
    start_date = '1972-07-23'
    end_date = '1982-02-26'


@register
class Landsat2RAW(Tier1, RAW, Landsat2):
    """ Landsat 2 Tier 1 """
    id = 'LANDSAT/LM02/C01/T1'
    short_name = 'L2RAW'
    bands = (Landsat1RAW.green, Landsat1RAW.red, Landsat1RAW.nir, Landsat1RAW.nir2,
             MSS.bqa)
    visualizers = Landsat1RAW.visualizers

    def __init__(self, **kwargs):
        super(Landsat2RAW, self).__init__(**kwargs)


@register
class Landsat2RAWT2(Tier2, Landsat2RAW):
    id = 'LANDSAT/LM02/C01/T2'
    short_name = 'L2RAWT2'

    def __init__(self, **kwargs):
        super(Landsat2RAWT2, self).__init__(**kwargs)


class Landsat3(MSS, Landsat):
    number = 3
    start_date = '1978-03-05'
    end_date = '1983-03-31'


@register
class Landsat3RAW(Tier1, RAW, Landsat3):
    id = 'LANDSAT/LM03/C01/T1'
    short_name = 'L3RAW'
    bands = (Landsat1RAW.green, Landsat1RAW.red, Landsat1RAW.nir, Landsat1RAW.nir2,
             MSS.bqa)
    visualizers = Landsat1RAW.visualizers

    def __init__(self, **kwargs):
        super(Landsat3RAW, self).__init__(**kwargs)


@register
class Landsat3RAWT2(Tier2, Landsat3RAW):
    id = 'LANDSAT/LM03/C01/T2'
    short_name = 'L3RAWT2'

    def __init__(self, **kwargs):
        super(Landsat3RAWT2, self).__init__(**kwargs)


class Landsat4MSS(MSS, Landsat):
    number = 4
    start_date = '1982-07-16'
    end_date = '1993-12-14'


@register
class Landsat4MSSRAW(Tier1, RAW, Landsat4MSS):
    """ Landsat 4 MSS """
    id = 'LANDSAT/LM04/C01/T1'
    short_name = 'L4MSSRAW'
    green = MSS.green(name='B1', **RAW._extra)
    red = MSS.red(name='B2', **RAW._extra)
    nir = MSS.nir(name='B3', **RAW._extra)
    nir2 = MSS.nir2(name='B4', **RAW._extra)
    bands = (green, red, nir, nir2, MSS.bqa)
    visualizers = (
        Visualization.falseColor([nir, red, green]),
    )

    def __init__(self, **kwargs):
        super(Landsat4MSSRAW, self).__init__(**kwargs)


@register
class Landsat4MSSRAWT2(Tier2, Landsat4MSSRAW):
    """ Landsat 4 MSS """
    id = 'LANDSAT/LM04/C01/T2'
    short_name = 'L4MSSRAWT2'

    def __init__(self, **kwargs):
        super(Landsat4MSSRAWT2, self).__init__(**kwargs)


class Landsat4TM(TM, Landsat):
    """ Landsat 4 TM Raw Tier 1 """
    number = 4
    start_date = '1982-07-16'
    end_date = '1993-12-14'


@register
class Landsat4RAW(Tier1, RAW, Landsat4TM):
    id = 'LANDSAT/LT04/C01/T1'
    short_name = 'L4RAW'
    blue = TM.blue(**RAW._extra)
    green = TM.green(**RAW._extra)
    red = TM.red(**RAW._extra)
    nir = TM.nir(**RAW._extra)
    swir = TM.swir(**RAW._extra)
    thermal = TM.thermal(**RAW._extra)
    swir2 = TM.swir2(**RAW._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)

    def __init__(self, **kwargs):
        super(Landsat4RAW, self).__init__(**kwargs)


@register
class Landsat4RAWT2(Tier2, Landsat4RAW):
    id = 'LANDSAT/LT04/C01/T2'
    short_name = 'L4RAWT2'

    def __init__(self, **kwargs):
        super(Landsat4RAWT2, self).__init__(**kwargs)


@register
class Landsat4TOA(Tier1, TOA, Landsat4TM):
    id = 'LANDSAT/LT04/C01/T1_TOA'
    short_name = 'L4TOA'
    blue = TM.blue(**TOA._extra)
    green = TM.green(**TOA._extra)
    red = TM.red(**TOA._extra)
    nir = TM.nir(**TOA._extra)
    swir = TM.swir(**TOA._extra)
    thermal = TM.thermal(**TOA._extra)
    swir2 = TM.swir2(**TOA._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)

    def __init__(self, **kwargs):
        super(Landsat4TOA, self).__init__(**kwargs)


@register
class Landsat4TOAT2(Tier2, Landsat4TOA):
    id = 'LANDSAT/LT04/C01/T2_TOA'
    short_name = 'L4TOAT2'

    def __init__(self, **kwargs):
        super(Landsat4TOAT2, self).__init__(**kwargs)


@register
class Landsat4SR(Tier1, SR, Landsat4TM):
    id = 'LANDSAT/LT04/C01/T1_SR'
    short_name = 'L4SR'
    blue = TM.blue(**SR._extra)
    green = TM.green(**SR._extra)
    red = TM.red(**SR._extra)
    nir = TM.nir(**SR._extra)
    swir = TM.swir(**SR._extra)
    thermal = TM.thermal(**SR._extra)
    swir2 = TM.swir2(**SR._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, SR.atm_op,
             SR.sr_cloud_qa, SR.pixel_qa, SR.radsat_qa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),
             Mask.fromBand('cloud_qa', SR.sr_cloud_qa))

    def __init__(self, **kwargs):
        super(Landsat4SR, self).__init__(**kwargs)


@register
class Landsat4SRT2(Tier2, Landsat4SR):
    id = 'LANDSAT/LT04/C01/T2_SR'
    short_name = 'L4SRT2'

    def __init__(self, **kwargs):
        super(Landsat4SRT2, self).__init__(**kwargs)


class Landsat5MSS(MSS, Landsat):
    number = 5
    start_date = '1984-01-01'
    end_date = '2012-05-05'


@register
class Landsat5MSSRAW(Tier1, RAW, Landsat4MSS):
    """ Landsat 5 MSS """
    id = 'LANDSAT/LM05/C01/T1'
    short_name = 'L5MSSRAW'
    green = MSS.green(name='B1', **RAW._extra)
    red = MSS.red(name='B2', **RAW._extra)
    nir = MSS.nir(name='B3', **RAW._extra)
    nir2 = MSS.nir2(name='B4', **RAW._extra)
    bands = (green, red, nir, nir2, MSS.bqa)
    visualizers = (
        Visualization.falseColor([nir, red, green]),
    )
    masks = (Mask.fromBand('BQA', MSS.bqa),)

    def __init__(self, **kwargs):
        super(Landsat5MSSRAW, self).__init__(**kwargs)


@register
class Landsat5MSSRAWT2(Tier2, Landsat5MSSRAW):
    """ Landsat 5 MSS """
    id = 'LANDSAT/LM05/C01/T2'
    short_name = 'L5MSSRAWT2'

    def __init__(self, **kwargs):
        super(Landsat5MSSRAWT2, self).__init__(**kwargs)


class Landsat5TM(TM, Landsat):
    """ Landsat 5 TM Raw Tier 1 """
    number = 5
    start_date = '1984-04-16'
    end_date = '2012-05-05'


@register
class Landsat5RAW(Tier1, RAW, Landsat5TM):
    id = 'LANDSAT/LT05/C01/T1'
    short_name = 'L5RAW'
    blue = TM.blue(**RAW._extra)
    green = TM.green(**RAW._extra)
    red = TM.red(**RAW._extra)
    nir = TM.nir(**RAW._extra)
    swir = TM.swir(**RAW._extra)
    thermal = TM.thermal(**RAW._extra)
    swir2 = TM.swir2(**RAW._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)

    def __init__(self, **kwargs):
        super(Landsat5RAW, self).__init__(**kwargs)


@register
class Landsat5RAWT2(Tier2, Landsat5RAW):
    id = 'LANDSAT/LT05/C01/T2'
    short_name = 'L5RAWT2'

    def __init__(self, **kwargs):
        super(Landsat5RAWT2, self).__init__(**kwargs)


@register
class Landsat5TOA(Tier1, TOA, Landsat5TM):
    id = 'LANDSAT/LT05/C01/T1_TOA'
    short_name = 'L5TOA'
    blue = TM.blue(**TOA._extra)
    green = TM.green(**TOA._extra)
    red = TM.red(**TOA._extra)
    nir = TM.nir(**TOA._extra)
    swir = TM.swir(**TOA._extra)
    thermal = TM.thermal(**TOA._extra)
    swir2 = TM.swir2(**TOA._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, TM.bqa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)

    def __init__(self, **kwargs):
        super(Landsat5TOA, self).__init__(**kwargs)


@register
class Landsat5TOAT2(Tier2, Landsat5TOA):
    id = 'LANDSAT/LT05/C01/T2_TOA'
    short_name = 'L5TOAT2'

    def __init__(self, **kwargs):
        super(Landsat5TOAT2, self).__init__(**kwargs)


@register
class Landsat5SR(Tier1, SR, Landsat5TM):
    id = 'LANDSAT/LT05/C01/T1_SR'
    short_name = 'L5SR'
    blue = TM.blue(**SR._extra)
    green = TM.green(**SR._extra)
    red = TM.red(**SR._extra)
    nir = TM.nir(**SR._extra)
    swir = TM.swir(**SR._extra)
    thermal = TM.thermal(**SR._extra)
    swir2 = TM.swir2(**SR._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, SR.atm_op,
             SR.sr_cloud_qa, SR.pixel_qa, SR.radsat_qa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),
             Mask.fromBand('cloud_qa', SR.sr_cloud_qa))

    def __init__(self, **kwargs):
        super(Landsat5SR, self).__init__(**kwargs)


@register
class Landsat5SRT2(Tier2, Landsat5SR):
    id = 'LANDSAT/LT05/C01/T2_SR'
    short_name = 'L5SRT2'

    def __init__(self, **kwargs):
        super(Landsat5SRT2, self).__init__(**kwargs)


class Landsat7ETM(ETM, Landsat):
    """ Landsat 7 ETM Raw Tier 1 """
    number = 7
    start_date = '1999-05-28'
    end_date = None


@register
class Landsat7RAW(Tier1, RAW, Landsat7ETM):
    id = 'LANDSAT/LE07/C01/T1'
    short_name = 'L7RAW'
    blue = ETM.blue(**RAW._extra)
    green = ETM.green(**RAW._extra)
    red = ETM.red(**RAW._extra)
    nir = ETM.nir(**RAW._extra)
    swir = ETM.swir(**RAW._extra)    
    swir2 = ETM.swir2(**RAW._extra)
    thermal1 = ETM.thermal_vcid_1(**RAW._extra)
    thermal2 = ETM.thermal_vcid_2(**RAW._extra)
    bands = (blue, green, red, nir, swir, thermal1, thermal2, swir2, TM.bqa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)

    def __init__(self, **kwargs):
        super(Landsat7RAW, self).__init__(**kwargs)


@register
class Landsat7RAWT2(Tier2, Landsat7RAW):
    id = 'LANDSAT/LE07/C01/T2'
    short_name = 'L7RAWT2'

    def __init__(self, **kwargs):
        super(Landsat7RAWT2, self).__init__(**kwargs)


@register
class Landsat7TOA(Tier1, TOA, Landsat7ETM):
    id = 'LANDSAT/LE07/C01/T1_TOA'
    short_name = 'L7TOA'
    blue = ETM.blue(**TOA._extra)
    green = ETM.green(**TOA._extra)
    red = ETM.red(**TOA._extra)
    nir = ETM.nir(**TOA._extra)
    swir = ETM.swir(**TOA._extra)
    swir2 = ETM.swir2(**TOA._extra)
    thermal1 = ETM.thermal_vcid_1(**TOA._extra)
    thermal2 = ETM.thermal_vcid_2(**TOA._extra)
    bands = (blue, green, red, nir, swir, thermal1, thermal2, swir2, TM.bqa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', TM.bqa),)

    def __init__(self, **kwargs):
        super(Landsat7TOA, self).__init__(**kwargs)


@register
class Landsat7TOAT2(Tier2, Landsat7TOA):
    id = 'LANDSAT/LE07/C01/T2_TOA'
    short_name = 'L7TOAT2'

    def __init__(self, **kwargs):
        super(Landsat7TOAT2, self).__init__(**kwargs)


@register
class Landsat7SR(Tier1, SR, Landsat7ETM):
    id = 'LANDSAT/LE07/C01/T1_SR'
    short_name = 'L7SR'
    blue = ETM.blue(**SR._extra)
    green = ETM.green(**SR._extra)
    red = ETM.red(**SR._extra)
    nir = ETM.nir(**SR._extra)
    swir = ETM.swir(**SR._extra)
    thermal = ETM.thermal(**SR._extra)
    swir2 = ETM.swir2(**SR._extra)
    bands = (blue, green, red, nir, swir, thermal, swir2, SR.atm_op,
             SR.sr_cloud_qa, SR.pixel_qa, SR.radsat_qa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),
             Mask.fromBand('cloud_qa', SR.sr_cloud_qa))

    def __init__(self, **kwargs):
        super(Landsat7SR, self).__init__(**kwargs)


@register
class Landsat7SRT2(Tier2, Landsat7SR):
    id = 'LANDSAT/LE07/C01/T2_SR'
    short_name = 'L7SRT2'

    def __init__(self, **kwargs):
        super(Landsat7SRT2, self).__init__(**kwargs)


class Landsat8OLI(OLI, Landsat):
    """ Landsat 8 OLI Raw Tier 1 """
    number = 8
    start_date = '2013-03-18'
    end_date = None


@register
class Landsat8RAW(Tier1, RAW, Landsat8OLI):
    id = 'LANDSAT/LC08/C01/T1'
    short_name = 'L8RAW'
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
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', OLI.bqa),)

    def __init__(self, **kwargs):
        super(Landsat8RAW, self).__init__(**kwargs)


@register
class Landsat8RAWT2(Tier2, Landsat8RAW):
    id = 'LANDSAT/LC08/C01/T2'
    short_name = 'L8RAWT2'

    def __init__(self, **kwargs):
        super(Landsat8RAWT2, self).__init__(**kwargs)


@register
class Landsat8TOA(Tier1, TOA, Landsat8OLI):
    id = 'LANDSAT/LC08/C01/T1_TOA'
    short_name = 'L8TOA'
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
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('BQA', OLI.bqa),)

    def __init__(self, **kwargs):
        super(Landsat8TOA, self).__init__(**kwargs)


@register
class Landsat8TOAT2(Tier2, Landsat8TOA):
    id = 'LANDSAT/LC08/C01/T2_TOA'
    short_name = 'L8TOAT2'

    def __init__(self, **kwargs):
        super(Landsat8TOAT2, self).__init__(**kwargs)


@register
class Landsat8SR(Tier1, SR, Landsat8OLI):
    id = 'LANDSAT/LC08/C01/T1_SR'
    short_name = 'L8SR'
    aerosol = OLI.aerosol(**SR._extra)
    blue = OLI.blue(**SR._extra)
    green = OLI.green(**SR._extra)
    red = OLI.red(**SR._extra)
    nir = OLI.nir(**SR._extra)
    swir = OLI.swir(**SR._extra)
    swir2 = OLI.swir2(**SR._extra)
    thermal1 = OLI.thermal(**SR._extra)
    thermal2 = OLI.thermal2(**SR._extra)
    bands = (aerosol, blue, green, red, nir, swir, thermal1, thermal2,
             OLI.SR.aerosol, SR.pixel_qa, OLI.SR.radsat_qa)
    visualizers = (
        Visualization.trueColor([red, green, blue]),
        Visualization.falseColor([nir, red, green]),
        Visualization.NSR([nir, swir, red])
    )
    masks = (Mask.fromBand('pixel_qa', SR.pixel_qa),)

    def __init__(self, **kwargs):
        super(Landsat8SR, self).__init__(**kwargs)


@register
class Landsat8SRT2(Tier2, Landsat8SR):
    id = 'LANDSAT/LC08/C01/T2_SR'
    short_name = 'L8SRT2'

    def __init__(self, **kwargs):
        super(Landsat8SRT2, self).__init__(**kwargs)