# coding=utf-8
""" Google Earth Engine Sentinel Collections """
from .visualization import *
from .datasets import OpticalSatellite, ImageCollection
from .bands import OpticalBand, BitBand, ClassificationBand, ExpressionBand, \
    RangeBand, Precisions
from .masks import Mask
from . import helpers
import geetools
from . import register
import ee


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
    _common = {'scale':0.0001, 'precision': Precisions.uint16}
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
                   precision=Precisions.uint16)
    qa20 = BitBand(name='QA20', alias='qa20', bits=dict(), resolution=20,
                   precision=Precisions.uint32)
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

    visualizers = Visualizers(
        NSR = Visualization.NSR([nir, swir, red]),
        TrueColor = Visualization.trueColor([red, green, blue]),
        FalseColor = Visualization.falseColor([nir, red, green]),
        RGB = Visualization.trueColor([red, green, blue]),
      )

    extra_bands = (ndvi, nbr, ndre)

    def add_cloud_probability(self, collection=None, name='cloud_probability'):
        """ add band from ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
        to the specified collection. If collection is None the complete S2
        collection will be used. You can use any collection as long as it keeps
        the image ids """
        if collection is None:
            collection = self.collection()

        clouds = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
        propertyField = 'system:index'
        Filter = ee.Filter.equals(leftField=propertyField,
                                  rightField=propertyField)
        join = ee.Join.saveFirst(matchKey='match', outer=False)
        joined = ee.ImageCollection(join.apply(collection, clouds, Filter))

        def overJoined(image):
            # properties = image.propertyNames()
            # retain = properties.remove('match')
            match = ee.Image(image.get('match')).rename(name)
            return image.addBands(match)

        return joined.map(overJoined)

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
                  negatives=('cloud', 'snow', 'shadow', 'water', 'cirrus'),
                  decoder=decode_hollstein))

    def add_scl(self, collection=None, name='SCL', renamed=False):
        """ Add the SCL band from S2SR collection to every image of the
        given collection. If collection is None it will use the complete S2TOA
        collection """
        s2sr = Sentinel2SR()
        S2SR_start = s2sr.start_date
        if collection is None:
            collection = self.collection().filterDate(S2SR_start, helpers.today())
        srcol = s2sr.collection()
        propertyField = 'system:index'
        Filter = ee.Filter.equals(leftField=propertyField,
                                  rightField=propertyField)
        join = ee.Join.saveFirst(matchKey='match', outer=False)
        joined = ee.ImageCollection(join.apply(collection, srcol, Filter))

        def overJoined(image):
            # properties = image.propertyNames()
            # retain = properties.remove('match')
            match = ee.Image(image.get('match'))\
                      .select(s2sr.scl.alias if renamed else s2sr.scl.name)\
                      .rename(name)
            return image.addBands(match)

        return joined.map(overJoined)

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
    aot = OpticalBand('AOT', 'aerosol_thickness', precision=Precisions.uint16,
                      resolution=10)
    wvp = OpticalBand('WVP', 'water_vapor_pressure', precision=Precisions.uint16,
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
        negatives=['saturated', 'dark', 'shadow', 'water',
                   'clouds_mid_probability', 'clouds_high_probability',
                   'cirrus', 'snow']
    )
    tci_r = OpticalBand(name='TCI_R', alias='true_color_red',
                        precision=Precisions.uint8, resolution=10)
    tci_g = OpticalBand(name='TCI_G', alias='true_color_blue',
                        precision=Precisions.uint8, resolution=10)
    tci_b = OpticalBand(name='TCI_B', alias='true_color_green',
                        precision=Precisions.uint8, resolution=10)

    cloud_prob = RangeBand('MSK_CLDPRB', 'cloud_probability', 'percentage',
                             precision=Precisions.uint8, resolution=20)
    snow_prob = RangeBand('MSK_SNWPRB', 'snow_probability', 'percentage',
                             precision=Precisions.uint8, resolution=10)


    # Add visualizers
    visualizers = Visualizers(
        NSR = Visualization.NSR([Sentinel2.nir, Sentinel2.swir, Sentinel2.red]),
        TrueColor = Visualization.trueColor([Sentinel2.red, Sentinel2.green, Sentinel2.blue]),
        FalseColor = Visualization.falseColor([Sentinel2.nir, Sentinel2.red, Sentinel2.green]),
        SCL = Visualization('SCL', [scl], 1, 11, palette=[
            '#ff0004', '#868686', '#774b0a', '#10d22c', '#ffff52', '#0000ff',
            '#818181', '#c0c0c0', '#f1f1f1', '#bac5eb', '#52fff9'
        ])
    )

    bands = (
        Sentinel2.aerosol, Sentinel2.blue, Sentinel2.green, Sentinel2.red,
        Sentinel2.red_edge_1, Sentinel2.red_edge_2, Sentinel2.red_edge_3,
        Sentinel2.nir, Sentinel2.red_edge_4, Sentinel2.water_vapor,
        Sentinel2.swir, Sentinel2.swir2, aot, wvp, scl, tci_r, tci_g, tci_b,
        cloud_prob, snow_prob, Sentinel2.qa10, Sentinel2.qa20, Sentinel2.qa60
    )

    masks = (
        Mask.fromBand('SCL', scl,
                      negatives=['cloud', 'saturated', 'dark', 'shadow', 'water',
                                 'clouds_mid_probability', 'clouds_high_probability',
                                 'cirrus', 'snow']),
        Mask.fromBand('QA60', Sentinel2.qa60)
    )

    def __init__(self, **kwargs):
        super(Sentinel2SR, self).__init__(**kwargs)


@register
class Sentinel1(ImageCollection):
    """ Sentinel 1 """
    id = 'COPERNICUS/S1_GRD'
    short_name = 'S1'
    start_date = '2014-10-03'
    end_date = None
    polarizations = ('HH', 'HH-HV', 'VV', 'VV-VH')
    HH = OpticalBand('HH', 'HH', 'db', precision=Precisions.double, resolution=10)
    HV = OpticalBand('HV', 'HV', 'db', precision=Precisions.double, resolution=10)
    VH = OpticalBand('VH', 'VH', 'db', precision=Precisions.double, resolution=10)
    VV = OpticalBand('VV', 'VV', 'db', precision=Precisions.double, resolution=10)
    angle = OpticalBand('angle', 'angle', 'degrees',
                        precision=Precisions.float, resolution=10)
    bands = (HH, HV, VH, VV, angle)
    masks = (Mask.empty(), )

    # https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel-1/radar_vegetation_index_code_dual_polarimetric/supplementary_material.pdf
    RVI = ExpressionBand(
        'RVI', 'rvi',
        """
        ( 4 * (10**(VH/10)) ) / 
        ( (10**(VV/10)) + (10**(VH/10)) )
        """,
        [VV, VH],
        precision=Precisions.float,
        resolution=10
    )

    extra_bands = (RVI,)

    def __init__(self, **kwargs):
        super(Sentinel1, self).__init__(**kwargs)

    def filter_polarisation(self, polarisation, collection=None):
        """ To get an homogeneous collection you have to filter using one of:

        - HH, HH-HV, VV or VV-VH

        HH and HH-HV collections cover the north and south pole. The rest of
        the world is covered by VV and VV-VH

        See: https://sentinel.esa.int/web/sentinel/missions/sentinel-1/observation-scenario

        :param polarisation: can be HH HH-HV VV or VV-VH
        :type polarisation: str
        :param collection: an image collection to filter. If None it will use
        the complete S1 collection
        :type collection: ee.ImageCollection
        """
        col = collection or self.collection()
        property = 'transmitterReceiverPolarisation'
        def filt(option):
            return ee.Filter.listContains(property, option)

        hh = filt('HH')
        hv = filt('HV')
        vv = filt('VV')
        vh = filt('VH')

        if isinstance(polarisation, str):
            filters = {
                'HH': ee.Filter.And(hh, hv.Not(), vv.Not(), vh.Not()),
                'HH-HV': ee.Filter.And(hh, hv, vv.Not(), vh.Not()),
                'VV': ee.Filter.And(hh.Not(), hv.Not(), vv, vh.Not()),
                'VV-VH': ee.Filter.And(hh.Not(), hv.Not(), vv, vh)
            }
            if polarisation not in filters:
                raise ValueError('Polarisation {} not available, options are: {}'.format(polarisation, filters.keys()))

            filter = filters[polarisation]
        else:
            if list(polarisation) == ['HH']:
                filter = ee.Filter.And(hh, hv.Not(), vv.Not(), vh.Not())
            elif list(polarisation) == ['HH', 'HV'] or list(polarisation) == ['HV', 'HH']:
                filter = ee.Filter.And(hh, hv, vv.Not(), vh.Not())
            elif list(polarisation) == ['VV']:
                filter = ee.Filter.And(hh.Not(), hv.Not(), vv, vh.Not())
            else:
                filter = ee.Filter.And(hh.Not(), hv.Not(), vv, vh)

        return col.filter(filter)


    def filter_orbit(self, ascending=True, collection=None):
        """ Filter by orbit ascending or descending """
        col = collection or self.collection()
        property = 'orbitProperties_pass'
        if ascending:
            return col.filterMetadata(property, 'equals', 'ASCENDING')
        else:
            return col.filterMetadata(property, 'equals', 'DESCENDING')

    def filter_instrument(self, instrument, collection=None):
        """ Options are 'IW' (Interferometric Wide Swath),
        'EW' (Extra Wide Swath) or 'SM' (Strip Map).
        See https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/acquisition-modes
        for details."""
        options = ('IW', 'EW', 'SM')
        if instrument not in options:
            raise ValueError('instrument must be one of {}'.format(options))
        property = 'instrumentMode'
        col = collection or self.collection()
        return col.filter(ee.Filter.eq(property, instrument))

    def convert_to_intensity(self, image):
        """ Convert values as they come in GEE (db) to intensity (DN) """
        rest = image.bandNames().remove('angle')
        intensity = ee.Image().expression('10**(v/10)', dict(v=image.select(rest)))
        return image.addBands(intensity, overwrite=True)

    def convert_to_db(self, image):
        """ Convert values from intensity to db """
        rest = image.bandNames().remove('angle')
        intensity = ee.Image().expression('log10(v)*10', dict(v=image.select(rest)))
        return image.addBands(intensity, overwrite=True)

    @staticmethod
    def _lee_filter(image):
        """ Lee speckle filter from

        https://groups.google.com/g/google-earth-engine-developers/c/puhSXaJ7NmI/m/TNOd0knyBQAJ
        """
        bandNames = image.bandNames()

        # image must be in natural units, i.e. not in dB!
        # Set up 3x3 kernels
        weights3 = ee.List.repeat(ee.List.repeat(1, 3), 3)
        kernel3 = ee.Kernel.fixed(3, 3, weights3, 1, 1, False)

        mean3 = image.reduceNeighborhood(ee.Reducer.mean(), kernel3)
        variance3 = image.reduceNeighborhood(ee.Reducer.variance(), kernel3)

        # Use a sample of the 3x3 windows inside a 7x7 windows to determine
        # gradients and directions
        sample_weights = ee.List(
            [[0,0,0,0,0,0,0],
             [0,1,0,1,0,1,0],
             [0,0,0,0,0,0,0],
             [0,1,0,1,0,1,0],
             [0,0,0,0,0,0,0],
             [0,1,0,1,0,1,0],
             [0,0,0,0,0,0,0]])

        sample_kernel = ee.Kernel.fixed(7, 7, sample_weights, 3, 3, False)

        # Calculate mean and variance for the sampled windows and store as 9 bands
        sample_mean = mean3.neighborhoodToBands(sample_kernel)
        sample_var= variance3.neighborhoodToBands(sample_kernel)

        # Determine the 4 gradients for the sampled windows
        gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
        gradients = gradients.addBands(
            sample_mean.select(6).subtract(sample_mean.select(2)).abs())
        gradients = gradients.addBands(
            sample_mean.select(3).subtract(sample_mean.select(5)).abs())
        gradients = gradients.addBands(
            sample_mean.select(0).subtract(sample_mean.select(8)).abs())

        # And find the maximum gradient amongst gradient bands
        max_gradient = gradients.reduce(ee.Reducer.max())

        # Create a mask for band pixels that are the maximum gradient
        gradmask = gradients.eq(max_gradient)

        # duplicate gradmask bands: each gradient represents 2 directions
        gradmask = gradmask.addBands(gradmask)

        # Determine the 8 directions
        directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
        directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
        directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
        directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
        # The next 4 are the not() of the previous 4
        directions = directions.addBands(directions.select(0).Not().multiply(5))
        directions = directions.addBands(directions.select(1).Not().multiply(6))
        directions = directions.addBands(directions.select(2).Not().multiply(7))
        directions = directions.addBands(directions.select(3).Not().multiply(8))

        # Mask all values that are not 1-8
        directions = directions.updateMask(gradmask)

        # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
        directions = directions.reduce(ee.Reducer.sum())

        sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

        # Calculate localNoiseVariance
        sigmaV = sample_stats.toArray().arraySort().arraySlice(0, 0, 5).arrayReduce(ee.Reducer.mean(), [0])

        # Set up the 7*7 kernels for directional statistics
        rect_weights = ee.List.repeat(ee.List.repeat(0, 7), 3).cat(ee.List.repeat(ee.List.repeat(1, 7), 4))

        diag_weights = ee.List(
            [[1,0,0,0,0,0,0],
             [1,1,0,0,0,0,0],
             [1,1,1,0,0,0,0],
             [1,1,1,1,0,0,0],
             [1,1,1,1,1,0,0],
             [1,1,1,1,1,1,0],
             [1,1,1,1,1,1,1]])

        rect_kernel = ee.Kernel.fixed(7, 7, rect_weights, 3, 3, False)
        diag_kernel = ee.Kernel.fixed(7, 7, diag_weights, 3, 3, False)

        # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
        dir_mean = image.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
        dir_var = image.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

        dir_mean = dir_mean.addBands(image.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
        dir_var = dir_var.addBands(image.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

        # and add the bands for rotated kernels
        for i in range(1, 4):
            dir_mean = dir_mean.addBands(image.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2 * i + 1)))
            dir_var= dir_var.addBands(image.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2 * i + 1)))
            dir_mean = dir_mean.addBands(image.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2 * i + 2)))
            dir_var= dir_var.addBands(image.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2 * i + 2)))


        # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
        dir_mean = dir_mean.reduce(ee.Reducer.sum())
        dir_var= dir_var.reduce(ee.Reducer.sum())

        # A finally generate the filtered value
        varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))

        b = varX.divide(dir_var)

        result = dir_mean.add(b.multiply(image.subtract(dir_mean)))

        result = result.arrayFlatten([bandNames])  #convert the resultant Array back to a multi-band image

        return(result)

    def lee_filter(self, image, renamed=False):
        """ Lee speckle filter from

        https://groups.google.com/g/google-earth-engine-developers/c/puhSXaJ7NmI/m/TNOd0knyBQAJ
        """
        # convert to natural units
        image = self.convert_to_intensity(image)

        angleB = self.getBandByName('angle')
        angle = angleB.alias if renamed else angleB.name
        bands = image.bandNames()
        bands = ee.List(
            ee.Algorithms.If(bands.contains(angle), bands.remove(angle), bands)
        )
        def wrap(band, im):
            im = ee.Image(im)
            band = ee.String(band)
            f = self._lee_filter(im.select(band))
            return im.addBands(f, overwrite=True)

        filtered = ee.Image(bands.iterate(wrap, image))\
                     .copyProperties(source=image)\
                     .set('system:time_start', image.date().millis())
        db = self.convert_to_db(ee.Image(filtered))
        return db


@register
class S2Cloudless(ImageCollection):
    id = 'COPERNICUS/S2_CLOUD_PROBABILITY'
    short_name = 'S2CLOUD'
    start_date = '2015-06-23'
    end_date = None
    probability = RangeBand('probability', 'probability', 0, 100,
                            Precisions.uint8, 10, 'percentage', 1)
    bands = (probability,)

    def addProbability(self, collection, rename=False):
        """ Combine this collection with the given S2 collection """
        # Join S2 SR with cloud probability dataset to add cloud mask.
        name = self.probability.alias if rename else self.probability.name
        matched = ee.Join.saveFirst(name).apply(
            primary=collection,
            secondary=self.collection(),
            condition=ee.Filter.equals(leftField='system:index',
                                       rightField='system:index')
        )

        def wrap(image):
            match = ee.Image(image.get(name))
            match = match.rename(name)
            return ee.Image(image).addBands(match)

        return ee.ImageCollection(matched.map(wrap))

    def apply(self, image, threshold, renamed=False):
        """ Apply the mask to an image using the
            given threshold. If the probability is higher than the threshold, it's
            considered as cloud and therefore masked """
        image = ee.Image(image)
        name = self.probability.alias if renamed else self.probability.name
        mask = image.select(name).lte(threshold)
        return ee.Image(image.updateMask(mask))

    def applyCollection(self, collection, threshold, renamed=False):
        """ Apply the mask to the images in a collection using the
        given threshold. If the probability is higher than the threshold, it's
        considered as cloud and therefore masked """
        # name = self.probability.alias if renamed else self.probability.name
        matched = self.addProbability(collection, renamed)
        return ee.ImageCollection(matched.map(lambda i: self.apply(i, threshold, renamed)))
