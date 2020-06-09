# coding=utf-8
from .helpers import *
import ee
import geetools


class Band:
    """ Bands """
    BAND_STR = """Band {name}({alias}):
resolution: {resolution}
precision: {precision}
units: {units}"""

    accepted_precisions = ('float', 'double', 'int8', 'uint8', 'uint16',
                           'int16', 'uint32', 'int32', 'int64', 'byte', None)

    def __init__(self, name, alias, precision=None, resolution=None,
                 units=None):
        """ Band object

        :param name: band original name (from GEE)
        :param alias: alias of the band, for example: Landsat 8 B4 is red
        :param precision: band type
        :param resolution: spatial resolution in meters
        :param units: the units of expression for the band. For example, 'db'
            for Sentinel 1
        """
        self.name = name
        self.alias = alias
        self.precision = precision
        self.resolution = resolution
        self.units = units

    def __str__(self):
        return self.BAND_STR.format(**self.__dict__, precision=self.precision)

    def test(self, image, renamed=False, verbose=False):
        """ Make a call to GEE to check if information matches """
        name = self.alias if renamed else self.name
        band = image.select(name)
        info = band.getInfo()
        precision = info['bands'][0]['data_type']['precision']
        localinfo = dataTypeInfo(self.precision)

        if precision == 'int':
            min = info['bands'][0]['data_type']['min']
            max = info['bands'][0]['data_type']['max']
            if verbose:
                print(localinfo['min'], min)
                print(localinfo['max'], max)
            try:
                assert localinfo['min'] == min
            except AssertionError:
                raise AssertionError('Band "{}" in GEE has a min of "{}" but "{}" in this instance'.format(self.name, min, localinfo['min']))

            try:
                assert localinfo['max'] == max
            except AssertionError:
                raise AssertionError('Band "{}" in GEE has a max of "{}" but "{}" in this instance'.format(self.name,max, localinfo['max']))

        if verbose:
            print(localinfo['precision'], precision)
        try:
            assert localinfo['precision'] == precision
        except AssertionError:
            raise AssertionError('Band "{}" in GEE has a precision of "{}" but "{}" in this instance'.format(self.name, precision, localinfo['precision']))


    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, value):
        if value not in self.accepted_precisions:
            msg = 'Presicion must be one of {}'
            raise ValueError(msg.format(self.accepted_precisions))
        self._precision = value

    def changePrecision(self, image, precision, renamed=False):
        if precision not in self.accepted_precisions:
            msg = 'Presicion must be one of {}'
            raise ValueError(msg.format(self.accepted_precisions))
        band = self.alias if renamed else self.name
        return image.cast(dict({band:ee.PixelType(precision)}))


class OpticalBand(Band):
    """ Special optical band """
    BAND_STR = """Band {name}({alias}):
resolution: {resolution}
precision: {precision}
units: {units}
scale: {scale}
wavelength: {wavelength}
"""
    def __init__(self, name, alias, units='reflectance', scale=1,
                 wavelength=None, **kwargs):
        super(OpticalBand, self).__init__(name, alias, units=units, **kwargs)
        self.wavelength = wavelength
        self.scale = scale


class BitBand(Band):
    """ Special bit band """
    BAND_STR = """Band {name}({alias}):
resolution: {resolution}
precision: {precision}
units: {units}
bits: {bits}
positives: {positives}
negatives: {negatives}
"""
    def __init__(self, name, alias, bits, positives=None, negatives=None,
                 **kwargs):
        """ By default all classes are considered as positive classes,
            but usually classes like 'clouds' or 'shadow' are negative. So you
            must specify it using the negatives parameter
        """
        super(BitBand, self).__init__(name, alias, **kwargs)
        self.bits = bits
        self.positives = positives
        self.negatives = negatives

        # check consistency
        if self.positives:
            for pos in self.positives:
                if pos not in self.options:
                    raise ValueError('Value {} not in {}'.format(pos, self.options))

        if self.negatives:
            for neg in self.negatives:
                if neg not in self.options:
                    raise ValueError('Value {} not in {}'.format(neg, self.options))

        if not self.positives and not self.negatives:
            self.positives = self.options

        if not self.positives:
            self.positives = [val for val in self.options if val not in self.negatives]

        if not self.negatives:
            self.negatives = [val for val in self.options if val not in self.positives]

    @property
    def options(self):
        """ A list of all classes """
        opt = list()
        for allclss in self.bits.values():
            for _, clss in allclss.items():
                opt.append(clss)
        return opt

    @property
    def reader(self):
        """ bit reader """
        return geetools.bitreader.BitReader(self.bits)

    def decodeImage(self, image, renamed=False):
        """ Decode the passed image using the bits information of this band """
        band = self.alias if renamed else self.name
        return self.reader.decodeImage(image, band)


class ClassificationBand(Band):
    """ band for classification """
    BAND_STR = """Band {name}({alias}):
resolution: {resolution}
precision: {precision}
units: {units}
classes: {classes}
positives: {positives}
negatives: {negatives}
"""
    def __init__(self, name, alias, classes=None, positives=None,
                 negatives=None, **kwargs):
        """" By default all classes are considered as positive classes,
        but usually classes like 'clouds' or 'shadow' are negative

        :param classes: classes are the keys (str). Values must be a string
            using JavaScript logical operators and the special word 'value'.
            For example (ndvi):
            {'soil': 'value>=0&&value<=0.4', 'vegetation': 'value>0.4'}

        :type classes: dict
        """
        super(ClassificationBand, self).__init__(name, alias, **kwargs)
        self.classes = classes
        self.positives = positives
        self.negatives = negatives

        # check consistency
        if self.positives:
            for pos in self.positives:
                if pos not in self.options:
                    raise ValueError('Value {} not in {}'.format(pos, self.options))

        if self.negatives:
            for neg in self.negatives:
                if neg not in self.options:
                    raise ValueError('Value {} not in {}'.format(neg, self.options))

        if not self.positives and not self.negatives:
            self.positives = self.options

        if not self.positives:
            self.positives = [val for val in self.options if val not in self.negatives]

        if not self.negatives:
            self.negatives = [val for val in self.options if val not in self.positives]

    @property
    def options(self):
        return list(self.classes.keys())

    def addClass(self, name, value, replace=False):
        """ add a class temporarily to the object """
        if not replace and name in self.options:
            raise ValueError('value {} already present'.format(name))
        self.classes[name] = value

    def decodeImage(self, image, renamed=False):
        """ Decode the classification band

        :return: an image in which each class is a band
        :rtype: ee.Image
        """
        band = self.alias if renamed else self.name
        image = image.select(band)
        classes = ee.Dictionary(self.classes)

        def wrap(key, value):
            exp = ee.String(value).cat('?1:0')
            mask = ee.Image(0).expression(exp, {'value': image})
            return mask.rename([key])

        newdict = classes.map(wrap)
        return geetools.tools.image.mixBands(newdict.values())


class ExpressionBand(Band):
    """ A band to add a new band """
    def __init__(self, name, alias, expression=None, bands=None, extra=None,
                 **kwargs):
        """ An expression band. Bands in list of bands must be instances of Band.
        Example (NDVI S2):

        ```
        s2sr = geetools.datasets.Sentinel2SR()
        ndviBand = ExpressionBand('NDVI', 'ndvi', '(nir-red)/(nir+red)',
                                  [s2sr.getBandByAlias('red'), s2sr.getBandByAlias('nir')],
                                  precision='float')

        # suppose we have a S2SR renamed image in variable i
        ndvi = ndviBand.apply(i, renamed=True)
        ```

        :param expression: the expression to apply when adding the band
        :type expression: str
        :param bands: the bands present in the expression. Elements in the list must be instances of Band
        :type bands: list
        :param extra: extra parameters for the expression
        :type extra: dict
        """
        super(ExpressionBand, self).__init__(name, alias, **kwargs)
        self.expression = expression
        self.bands = bands
        self.extra = extra
        self.resolution = kwargs.get(
            'resolution',
            min([b.resolution for b in self.bands if isinstance(b.resolution, (int, float))])
        )

    def apply(self, image, renamed=False):
        """ Apply the band expression to the parsed image """
        params = dict()
        for band in self.bands:
            name = band.alias if renamed else band.name
            params[band.alias] = image.select(name)
            params[band.name] = image.select(name)

        if self.extra:
            for key, value in self.extra.items():
                params.setdefault(key, value)
        name = self.alias if renamed else self.name
        final = ee.Image(0).expression(self.expression, params).rename(name)
        if self.precision:
            final = self.changePrecision(final, self.precision, renamed)
        return final
