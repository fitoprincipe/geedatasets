# coding=utf-8
from datetime import date
import ee
import geetools

TODAY = date.today().isoformat()


def get_mask(image, flag, value):
    """ Get a mask using

    - eq, neq, gt, gte, lt, lte
    """
    mask = ee.Algorithms.If(
        flag.compareTo('eq').eq(0), image.eq(value),
        ee.Algorithms.If(
            flag.compareTo('neq').eq(0), image.neq(value),
            ee.Algorithms.If(
                flag.compareTo('gt').eq(0), image.gt(value),
                ee.Algorithms.If(
                    flag.compareTo('gte').eq(0), image.gte(value),
                    ee.Algorithms.If(
                        flag.compareTo('lt').eq(0), image.lt(value),
                        ee.Algorithms.If(
                            flag.compareTo('lte').eq(0), image.lte(value),
                            image
    ))))))
    return ee.Image(mask)


def allequal(iterable):
    """ Check if all elements inside an iterable are equal """
    first = iterable[0]
    rest = iterable[1:]
    for item in rest:
        if item == first: continue
        else: return False
    return True


def dataTypeInfo(type):
    def make(bits, signed=True):
        nmax = 2**bits
        if signed:
            min = int(-nmax/2)
            max = int(nmax/2-1)
        else:
            min = 0
            max = int(nmax-1)
        return dict(precision='int', min=min, max=max)

    data = {
        'float': dict(precision='float'),
        'double': dict(precision='double'),
        'int8': make(8, True),
        'uint8': make(8, False),
        'int16': make(16, True),
        'uint16': make(16, False),
        'int32': make(32, True),
        'uint32': make(32, False),
        'int64': make(64, True),
    }
    lowtype = type.lower()
    return data.get(lowtype)


def convertPrecision(image, precision):
    """ Convert data type precisions

    DEPRECATED: use ee.Image.cast instead
    """
    TYPES = ee.Dictionary({'float': image.toFloat(),
                           'double': image.toDouble(),
                           'int8': image.toInt8(),
                           'uint8': image.toUint8(),
                           'uint16': image.toUint16(),
                           'int16': image.toInt16(),
                           'uint32': image.toUint32(),
                           'int32': image.toInt32(),
                           'int64': image.toInt64()
                           })
    return ee.Image(TYPES.get(precision.lower()))


def convertPrecisions(image, precision_dict):
    """
    DEPRECATED: use ee.Image.cast instead
    """
    precisions = ee.Dictionary(precision_dict)
    bands = ee.List(precisions.keys())
    def iteration(band, ini):
        ini = ee.Image(ini)
        imgband = ini.select([band])
        precision = ee.String(precisions.get(band))
        newband = convertPrecision(imgband, precision)
        return geetools.tools.image.replace(ini, band, newband)

    return ee.Image(bands.iterate(iteration, image))


def info(collection, renamed=False):
    """ Get the information for the parsed collection

    :param collection: the collection to get the information from
    :type collection: Collection
    :rtype: dict
    """
    data = dict()
    data['spacecraft'] = collection.spacecraft
    data['id'] = collection.name
    data['bands'] = [b.name for b in collection.bands]
    data['band_alias'] = [b.alias for b in collection.bands]
    data['scales'] = collection.scales(renamed=renamed)
    data['resolutions'] = collection.resolutions(renamed=renamed)
    data['start_date'] = collection.start_date
    data['end_date'] = collection.end_date
    data['algorithms'] = collection.algorithms
    if not renamed:
        data['bit_bands'] = [b.name for b in collection.bitBands]
        data['optical_bands'] = [b.name for b in collection.opticalBands]
        data['classification_bands'] = [b.name for b in collection.classificationBands]
    else:
        data['bit_bands'] = [b.alias for b in collection.bitBands]
        data['optical_bands'] = [b.alias for b in collection.opticalBands]
        data['classification_bands'] = [b.alias for b in collection.classificationBands]

    data['cloud_cover'] = collection.cloud_cover
    data['ee_collection'] = collection.collection
    data['indices'] = collection.indices

    if (collection.spacecraft == 'LANDSAT'):
        data['sensor'] = collection.sensor
        data['process'] = collection.process
        data['tier'] = collection.tier
        data['number'] = collection.number

    if (collection.spacecraft == 'SENTINEL'):
        data['number'] = collection.number
        if collection.number == 2:
            data['process'] = collection.process
            if collection.process == 'SR':
                data['SCL'] = collection.SCL_data
    return data


def infoEE(collection):
    """ Return an information ee.Dictionary """
    information = info(collection)
    information.pop('algorithms')
    information.pop('ee_collection')
    information.pop('indices')
    information.pop('visualization')

    return ee.Dictionary(information)


def getCommonBands(*datasets, match='name'):
    """ Get the common bands of the parsed collections

    :param match: the field to match, can be: name or alias
    :type match: str
    """
    renamed = True if match == 'alias' else False
    allbands = [set(ds.bandNames(renamed)) for ds in datasets]
    set0 = allbands[0]
    for bset in allbands[1:]:
        set0 = set0.intersection(bset)

    return list(set0)