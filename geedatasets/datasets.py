# coding=utf-8
import ee
from .bands import *
import geetools
import datetime


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
    extra_bands = tuple()

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
    masks = tuple()
    extra_bands = tuple()

    def start(self):
        """ Get the start date

        :return: start date as a date object
        :rtype: datetime.date
        """
        if self.start_date:
            return datetime.date.fromisoformat(self.start_date)
        else:
            return None

    def end(self):
        """ Get the end date. If None, it will return the date for today

        :return: end date as a date object
        :rtype: datetime.date
        """
        if self.end_date is None:
            return datetime.date.today()
        else:
            return datetime.date.fromisoformat(self.end_date)

    @property
    def all_bands(self):
        return list(self.bands)+list(self.extra_bands)

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

        for band in self.all_bands:
            band.test(image, renamed, verbose)

    def info(self):
        # BANDS
        if self.all_bands:
            bands = ["{} ({})".format(b.name, b.alias) for b in self.bands]
        else:
            bands = None

        # VISUALIZERS
        if self.visualizers:
            visualizers = [vis.name for vis in self.visualizers]
        else:
            visualizers = None

        # MASKS
        if self.masks:
            masks = [mask.name for mask in self.masks]
        else:
            masks = None

        params = {attr: getattr(self, attr) for attr in dir(self)}
        params['bands'] = bands
        params['visualizers'] = visualizers
        params['masks'] = masks
        return self.INFO.format(**params)

    def collection(self, bounds=None, date=None):
        """ Google Earth Engine Original Image Collection.
        You can filter by bounds and date using parameters
        """
        col = self.eeObject()
        if isinstance(bounds, (ee.Feature, ee.FeatureCollection, ee.Image)):
            bounds = bounds.geometry()
        if bounds:
            unbounded = bounds.isUnbounded()
            col = ee.ImageCollection(
                ee.Algorithms.If(unbounded, col, col.filterBounds(bounds)))

        if isinstance(date, (list, tuple)):
            col = col.filterDate(ee.Date(date[0]), ee.Date(date[1]))
        elif isinstance(date, (str, ee.Date)):
            col = col.filterDate(ee.Date(date))
        return col

    def getMask(self, name):
        """ Get a mask by its name. The mask object must have a method called
        `apply(image, negatives, positives, renamed)`
        """
        if not self.masks:
            return None
        for mask in self.masks:
            if mask.name == name:
                return mask

    def _getBandsByType(self, reference):
        return [band for band in self.all_bands if isinstance(band, reference)]

    @property
    def opticalBands(self):
        return self._getBandsByType(OpticalBand)

    @property
    def bitBands(self):
        return self._getBandsByType(BitBand)

    @property
    def classificationBands(self):
        return self._getBandsByType(ClassificationBand)

    def bandNames(self, renamed=False):
        """ List of band names """
        bands = [band.alias for band in self.bands] if renamed else [band.name for band in self.bands]
        return bands

    def bandAlias(self):
        """ List of band alias """
        return [band.alias for band in self.bands]

    def extraBandNames(self, renamed=False):
        """ List of extra band names """
        bands = [band.alias for band in self.extra_bands] if renamed else [band.name for band in self.extra_bands]
        return bands

    def precisions(self, renamed=False):
        """ dict of precisions using band ids as keys if not renamed or
        band names if renamed """
        precisions_dict = {}
        for band in self.all_bands:
            name = band.alias if renamed else band.name
            precisions_dict[name] = band.precision
        return precisions_dict

    def scales(self, renamed=False):
        """ dict of scales using band ids as keys
        if not renamed or band names if renamed """
        scales_dict = {}
        for band in self.all_bands:
            name = band.alias if renamed else band.name
            scales_dict[name] = band.scale

        return scales_dict

    def resolutions(self, renamed=False):
        """ dict of scales using band ids as keys if not renamed or
        band names if renamed """
        res_dict = {}
        for band in self.all_bands:
            name = band.alias if renamed else band.name
            res_dict[name] = band.resolution
        return res_dict

    def minResolution(self):
        """ Get minimal resolution between all bands """
        res = [band.resolution for band in self.bands]
        return min(res)

    def getBandByAlias(self, alias):
        data = None
        for b in self.all_bands:
            bid = b.alias
            if bid == alias:
                data = b
        return data

    def getBandByName(self, name):
        data = None
        for b in self.all_bands:
            bid = b.name
            if bid == name:
                data = b
        return data

    def getBand(self, band, by='name'):
        """ Get a band by name or alias """
        return self.getBandByAlias(band) if by.lower() == 'alias' else self.getBandByName(band)

    def visualization(self, option, renamed=False):
        """ Return visualization parameters for ui.Map.addLayer.

        :param option: see Visualization class
        """
        options = [vis.name for vis in self.visualizers]
        if option not in options:
            raise ValueError('{} not in visualizers. Options are {}'.format(option, options))
        visualizers = [vis for vis in self.visualizers if vis.name == option]
        return visualizers[0].params(renamed)

    def addBand(self, image, band, renamed=False):
        """ Add an expression band to an image """
        if not isinstance(band, ExpressionBand):
            raise ValueError('Only a ExpressionBand can be added')
        return band.apply(image, renamed)

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

    def applyMask(self, image, mask=None, negatives=None, positives=None,
                  renamed=False):
        """ Get a mask image using the passed mask band and apply it to the
        passed image

        :param image: the image to get the mask from
        :type image: ee.Image
        :param mask: the mask band
        :type mask: str
        :param classes: the classes for the mask. For example: ['cloud', 'shadow']
        :type classes: list
        :return: The image passed with pixels masked
        :rtype: ee.Image
        """
        if mask is None:
            f = self.masks[0]
            return f(image, negatives, positives, renamed)
        else:
            f = self.getMask(mask)
            return f.apply(image, negatives, positives, renamed)

    def rename(self, image):
        """ Rename bands """
        original_names = {band.name: band.alias for band in self.bands}
        return geetools.tools.image.renameDict(image, original_names)

    def proxyImage(self, renamed=False, masked=True):
        """ Create an Image with the band names, type and scale but empty """
        precisions = self.precisions(renamed=renamed)
        first_band = self.bands[0]
        if not renamed:
            name = first_band.name
        else:
            name = first_band.alias

        init = ee.Image.constant(0).rename(name)
        if masked:
            init = init.selfMask()
        init = convertPrecision(init, precisions[name])
        for i, band in enumerate(self.bands):
            if i == 0: continue
            if not renamed:
                name = band.name
            else:
                name = band.alias
            img = ee.Image.constant(0).rename(name)
            if masked:
                img = img.selfMask()
            img = convertPrecision(img, precisions[name])
            init = init.addBands(img)

        return init

    def rescale(self, image, target, renamed=False, match='alias', drop=False):
        """ Re-scale the values of image which must belong to collection so the
            values match the ones from collection_from

        :param dataset: The Dataset to which belongs the image
        :type dataset: ImageCollection
        :param target: the Dataset to get the range from
        :type target: ImageCollection
        :param renamed: is the parsed image renamed?
        :type renamed: bool
        :param match: whether to match by 'name' or 'alias'
        :type match: str
        :param drop: drop the non matching bands?
        :type drop: bool
        """
        # Create comparative collection
        common_bands = getCommonBands(self, target, match=match)
        if renamed:
            bands = [self.getBand(b, match).alias for b in common_bands]
        else:
            bands = [self.getBand(b, match).name for b in common_bands]

        scales = dict()
        scales_match = dict()
        precisions = dict()
        precisions_match = dict()

        # fill common bands scales and precisions
        for band in common_bands:
            b1 = self.getBand(band, match)
            bset = b1.alias if renamed else b1.name

            scales[bset] = b1.scale
            precisions[bset] = b1.precision

            # bset comes from dataset because image comes from dataset
            b2 = target.getBand(band, match)
            scales_match[bset] = b2.scale
            precisions_match[bset] = b2.precision

        # bands in the image
        bands = ee.List(bands)
        scales = ee.Dictionary(scales)
        scales_match = ee.Dictionary(scales_match)

        def iteration(band, ini):
            ini = ee.Image(ini)
            band = ee.String(band)

            iband = ini.select(band)
            scale = ee.Number(scales.get(band))
            scale_match = ee.Number(scales_match.get(band))
            newband = iband.multiply(scale).divide(scale_match)

            return ini.addBands(newband, overwrite=True)

        final = ee.Image(bands.iterate(iteration, image))
        final = final.cast(precisions_match)

        if drop:
            final = final.select(bands)
        return final


class Table(Dataset):
    INFO = """ID: {id}
"""
    type = ee.FeatureCollection

    def collection(self):
        return self.eeObject()


class OpticalSatellite:
    spacecraft = None
    cloud_cover = None