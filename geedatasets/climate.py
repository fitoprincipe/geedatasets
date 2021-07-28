# coding=utf-8
""" Google Earth Engine Climate Collections """
from .visualization import *
from .datasets import OpticalSatellite, ImageCollection, Image
from .bands import OpticalBand, BitBand, ClassificationBand, ExpressionBand, \
    RangeBand, Precisions
from .masks import Mask
from . import helpers
import geetools
from . import register
import ee


@register
class WorldClim(ImageCollection):
    INFO = """ID: {id}
short name: {short_name}"""
    id = 'WORLDCLIM/V1/MONTHLY'
    short_name = 'WCLIM'

    tavg = OpticalBand('tavg', 'tavg', units='°C', scale=0.1,
                       resolution=1000, precision=Precisions.float,
                       description='Average Temperature')
    tmin = OpticalBand('tmin', 'tmin', units='°C', scale=0.1,
                       resolution=1000, precision=Precisions.float,
                       description='Min Temperature')
    tmax = OpticalBand('tmax', 'tmax', units='°C', scale=0.1,
                       resolution=1000, precision=Precisions.float,
                       description='Max Temperature')
    prec = OpticalBand('prec', 'prec', units='mm', scale=1,
                       resolution=1000, precision=Precisions.float,
                       description='Precipitation')

    bands = [tavg, tmin, tmax, prec]

    _months = {
        'january': '01', 'february': '02', 'march': '03', 'april': '04',
        'may': '05', 'june': '06', 'july': '07', 'august': '08',
        'september': '09', 'october': '10', 'november': '11', 'december': '12'
    }

    def getMonth(self, month):
        """ Get a month image passing the month name (ex. 'january') or number """
        if isinstance(month, (str,)):
            m = self._months[month]
        else:
            m = '0{}'.format(month) if month < 10 else '{}'.format(month)

        return ee.Image('{}/{}'.format(self.id, m))

    def annualPrecipitation(self):
        """ Get annual precipitation """
        reduced = self.collection().select(self.prec.name).reduce(ee.Reducer.sum())
        return ee.Image(reduced).rename('mm')

    