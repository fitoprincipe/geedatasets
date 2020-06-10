# coding=utf-8
from .bands import BitBand, ClassificationBand
import ee


class Mask:
    def __init__(self, name, options=None, negatives=None, positives=None,
                 decoder=None):
        """ Mask object. It can be created using a bit band or a
        classification band, or can be a custom made via a `decoder`.
        When providing a `decoder` it must be a function take takes the
        following arguments:

        - image: the image to apply the mask
        - renamed: a bool indicating if the image has been renamed

        and returns an image with one band per category.

        The `band` argument has priority over the `decoder`.

        The arguments `negatives` and `positives` must be provided for setting
        the default behaviour, but it can be changed in the `apply` method
        """
        self.name = name
        self.negatives = negatives
        self.positives = positives
        self.decodeImage = decoder
        self._options = options

    def __str__(self):
        return """'{name}' Mask
options: {options}
default negatives: {negatives}
default positives: {positives}
""".format(name=self.name, options=self.options,
           negatives=self.negatives, positives=self.positives)

    @property
    def options(self):
        negatives = self.negatives or []
        positives = self.positives or []
        if self._options:
            return self._options
        else:
            return list(positives)+list(negatives)

    def getNegative(self, image, classes=None, renamed=False):
        """ Pixels with value = 1 will be masked """
        if not classes:
            classes = self.negatives

        decoded = self.decodeImage(image, renamed).select(classes)
        mask = decoded.reduce(ee.Reducer.sum()).rename('mask').Not()
        return mask

    def getPositive(self, image, classes=None, renamed=False):
        """ Pixels with value = 0 will be masked """
        if not classes:
            classes = self.positives

        decoded = self.decodeImage(image, renamed).select(classes)
        mask = decoded.reduce(ee.Reducer.sum()).rename('mask')
        return mask

    def get(self, image, negatives=None, positives=None, renamed=False):
        """ Get a mask """
        if positives and negatives:
            positives = [pos for pos in positives if pos in self.options]
            negatives = [neg for neg in negatives if neg in self.options]
            pmask = self.getPositive(image, positives, renamed)
            nmask = self.getNegative(image, negatives, renamed)
            mask = pmask.And(nmask)
        elif positives:
            positives = [pos for pos in positives if pos in self.options]
            mask = self.getPositive(image, positives, renamed)
        elif negatives:
            negatives = [neg for neg in negatives if neg in self.options]
            mask = self.getNegative(image, negatives, renamed)
        else:
            if self.positives and self.negatives:
                pmask = self.getPositive(image, self.positives, renamed)
                nmask = self.getNegative(image, self.negatives, renamed)
                mask = pmask.And(nmask)
            elif self.positives:
                mask = self.getPositive(image, self.positives, renamed)
            elif self.negatives:
                mask = self.getNegative(image, self.negatives, renamed)
            else:
                mask = ee.Image(1).rename('mask')
        return mask

    def apply(self, image, negatives=None, positives=None, renamed=False):
        """ Apply the mask for positives and negatives """
        mask = self.get(image, negatives, positives, renamed)
        return image.updateMask(mask)

    @classmethod
    def fromBand(cls, name, band, options=None, negatives=None, positives=None):
        """ Make a Mask from a QA Band """
        return cls(name, options or band.options,
                   negatives or band.negatives,
                   positives or band.positives,
                   band.decodeImage)
