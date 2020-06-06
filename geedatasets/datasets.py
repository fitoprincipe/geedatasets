# coding=utf-8
import ee
from .bands import *
import geetools


class Dataset:
    """ Parent class for common operations """
    id = None
    short_name = None
    type = None

    def __init__(self, **kwargs):
        pass

    def eeObject(self):
        """ Create the EE object using the class Type and ID"""
        return self.type(self.id)

    def info(self):
        return None

    def __str__(self):
        return self.info()


class Image(Dataset):
    INFO = """ID: {id}
Bands: {bands}
Visualizers: {visualizers}
Date: {date}
"""
    type = ee.Image
    bands = tuple()
    date = None
    region = None
    visualizers = tuple()

    def test(self, image=None, renamed=False, verbose=False):
        if not image:
            renamed = False
            image = self.image
        # Check bands
        bands_gee = image.bandNames().getInfo()
        bands_here = [b.alias for b in self.bands] if renamed else [b.name for b in self.bands]
        for band in bands_here:
            if band not in bands_gee:
                raise AssertionError('Band {} not present in GEE'.format(band))

        for band in self.bands:
            band.test(image, renamed, verbose)

        # Check date
        date_gee = image.date().format('yMMdd').getInfo()
        try:
            assert date_gee == self.date
        except AssertionError:
            raise AssertionError('date in GEE is {} and in this instance {}'.format(date_gee, self.date))

    def info(self):
        if self.bands:
            bands = ["{} ({})".format(b.name, b.alias) for b in self.bands]
        else:
            bands = None

        if self.visualizers:
            visualizers = [vis.name for vis in self.visualizers]
        else:
            visualizers = None

        params = {attr: getattr(self, attr) for attr in dir(self)}
        params['bands'] = bands
        params['visualizers'] = visualizers
        return self.INFO.format(**params)

    def image(self):
        return self.eeObject()

    @property
    def opticalBands(self):
        return self.getBandsByType(OpticalBand)

    @property
    def bitBands(self):
        return self.getBandsByType(BitBand)

    @property
    def classificationBands(self):
        return self.getBandsByType(ClassificationBand)

    def addBand(self, image, band, renamed=False):
        """ Add an expression band """
        if not isinstance(band, ExpressionBand):
            raise ValueError('Only a ExpressionBand can be added')
        return band.apply(image, renamed)

    def getBandsByType(self, reference):
        return [band for band in self.bands if isinstance(band, reference)]

    def bandNames(self, renamed=False):
        """ List of band names """
        if not renamed:
            bands = [band.name for band in self.bands]
        else:
            bands = [band.alias for band in self.bands]

        return bands

    def precisions(self, renamed=False):
        """ dict of precisions using band ids as keys if not renamed or
        band names if renamed """
        precisions_dict = {}
        for band in self.bands:
            name = band.alias if renamed else band.name
            precisions_dict[name] = band.precision
        return precisions_dict

    def scales(self, renamed=False):
        """ dict of ranges (min and max values as dict) using band ids as keys
        if not renamed or band names if renamed """
        scales_dict = {}
        for band in self.bands:
            name = band.alias if renamed else band.name
            scales_dict[name] = band.scale

        return scales_dict

    def resolutions(self, renamed=False):
        """ dict of scales using band ids as keys if not renamed or
        band names if renamed """
        res_dict = {}
        for band in self.bands:
            name = band.alias if renamed else band.name
            res_dict[name] = band.resolution
        return res_dict

    def getBandByAlias(self, alias):
        data = None
        for b in self.bands:
            bid = b.alias
            if bid == alias:
                data = b
        return data

    def getBandByName(self, name):
        data = None
        for b in self.bands:
            bid = b.name
            if bid == name:
                data = b
        return data


class ImageCollection(Dataset):
    INFO = """ID: {id}
Start Date: {start_date}
End Date: {end_date}
Bands: {bands}
Visualizers: {visualizers}
"""
    type = ee.ImageCollection
    bands = tuple()
    start_date = None
    end_date = None
    region = None
    visualizers = tuple()

    def test(self, image=None, renamed=False, verbose=False):
        if not image:
            renamed = False
            image = ee.Image(self.collection.first())
        # Check bands
        bands_gee = image.bandNames().getInfo()
        bands_here = [b.alias for b in self.bands] if renamed else [b.name for b in self.bands]
        for band in bands_here:
            if band not in bands_gee:
                raise AssertionError('Band {} not present in GEE'.format(band))

        if verbose:
            print('Image bands: {} \nDataset bands: {}'.format(bands_gee, bands_here))

        for band in self.bands:
            band.test(image, renamed, verbose)

    def info(self):
        if self.bands:
            bands = ["{} ({})".format(b.name, b.alias) for b in self.bands]
        else:
            bands = None

        if self.visualizers:
            visualizers = [vis.name for vis in self.visualizers]
        else:
            visualizers = None

        params = {attr: getattr(self, attr) for attr in dir(self)}
        params['bands'] = bands
        params['visualizers'] = visualizers
        return self.INFO.format(**params)

    def collection(self):
        """ Google Earth Engine Original Image Collection """
        return self.eeObject()

    @property
    def opticalBands(self):
        return self.getBandsByType(OpticalBand)

    @property
    def bitBands(self):
        return self.getBandsByType(BitBand)

    @property
    def classificationBands(self):
        return self.getBandsByType(ClassificationBand)

    def addBand(self, image, band, renamed=False):
        """ Add an expression band """
        if not isinstance(band, ExpressionBand):
            raise ValueError('Only a ExpressionBand can be added')
        return band.apply(image, renamed)

    def getBandsByType(self, reference):
        return [band for band in self.bands if isinstance(band, reference)]

    def bandNames(self, renamed=False):
        """ List of band names """
        if not renamed:
            bands = [band.name for band in self.bands]
        else:
            bands = [band.alias for band in self.bands]

        return bands

    def precisions(self, renamed=False):
        """ dict of precisions using band ids as keys if not renamed or
        band names if renamed """
        precisions_dict = {}
        for band in self.bands:
            name = band.alias if renamed else band.name
            precisions_dict[name] = band.precision
        return precisions_dict

    def scales(self, renamed=False):
        """ dict of ranges (min and max values as dict) using band ids as keys
        if not renamed or band names if renamed """
        scales_dict = {}
        for band in self.bands:
            name = band.alias if renamed else band.name
            scales_dict[name] = band.scale

        return scales_dict

    def resolutions(self, renamed=False):
        """ dict of scales using band ids as keys if not renamed or
        band names if renamed """
        res_dict = {}
        for band in self.bands:
            name = band.alias if renamed else band.name
            res_dict[name] = band.resolution
        return res_dict

    def bitImage(self, image, band, renamed=False):
        """ Get an image from the bit information from the qa band

        :param band: the quality band name
        :type band: str
        :param image: the image to decode with the qa band
        :type image: ee.Image
        :return: the image with the decode bands added
        """
        b = self.getBandByAlias(band) if renamed else self.getBandByName(band)
        if not isinstance(b, (BitBand, ClassificationBand)):
            raise ValueError('The band must be a bit or classification')

        return b.decodeImage(image, renamed)

    def applyMask(self, image, mask_band=None, negatives=None, positives=None,
                  renamed=False):
        """ Get a mask image using the passed mask band and apply it to the
        passed image

        :param image: the image to get the mask from
        :type image: ee.Image
        :param mask_band: the mask band
        :type mask_band: str
        :param classes: the classes for the mask. For example: ['cloud', 'shadow']
        :type classes: list
        :return: The image passed with pixels masked
        :rtype: ee.Image
        """
        if mask_band is None:
            f = self.common_masks[0]
            return f(image, negatives, positives, renamed)
        else:
            band = self.getBandByAlias(mask_band) if renamed else self.getBandByName(mask_band)
            if not isinstance(band, (BitBand, ClassificationBand)):
                raise ValueError('The band must be a bit or classification band')

            return band.applyMask(image, negatives, positives, renamed)

    def getBandByAlias(self, alias):
        data = None
        for b in self.bands:
            bid = b.alias
            if bid == alias:
                data = b
        return data

    def getBandByName(self, name):
        data = None
        for b in self.bands:
            bid = b.name
            if bid == name:
                data = b
        return data

    def visualization(self, option, renamed=False):
        """ Return visualization parameters for ui.Map.addLayer.

        :param option: see Visualization class
        """
        options = [vis.name for vis in self.visualizers]
        if option not in options:
            raise ValueError('{} not in visualizers. Options are {}'.format(option, options))
        visualizers = [vis for vis in self.visualizers if vis.name == option]
        return visualizers[0].params(renamed)

    def rename(self, image, reference='all'):
        """ Rename bands according to the parsed reference. It can be:
        optical, thermal, bits, all """
        if reference == 'all':
            original_names = {band.name: band.alias for band in self.bands}
        else:
            original_names = {band.name: band.alias for band in self.bands if band.reference == reference}

        return geetools.tools.image.renameDict(image, original_names)

    def proxyImage(self, renamed=False):
        """ Create an Image with the band names, type and scale but empty """
        precisions = self.precisions(renamed=renamed)
        first_band = self.bands[0]
        if not renamed:
            name = first_band.name
        else:
            name = first_band.alias

        init = ee.Image.constant(0).rename(name)
        init = convertPrecision(init, precisions[name])
        for i, band in enumerate(self.bands):
            if i == 0: continue
            if not renamed:
                name = band.name
            else:
                name = band.alias

            img = ee.Image.constant(0).rename(name)
            img = convertPrecision(img, precisions[name])
            init = init.addBands(img)

        return init

class Table(Dataset):
    INFO = """ID: {id}
"""
    type = ee.FeatureCollection

    def collection(self):
        return self.eeObject()


class OpticalSatellite:
    spacecraft = None
    short_name = None
    cloud_cover = None
    algorithms = None
    # common masks elements must be functions with 4 args: image, negatives, positives, renamed
    common_masks = [lambda i, classes, renamed, negatives: i]