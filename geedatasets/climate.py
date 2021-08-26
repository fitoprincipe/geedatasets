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
import datetime


@register
class WorldClim(ImageCollection):
    INFO = """ID: {id}
short name: {short_name}"""
    id = 'WORLDCLIM/V1/MONTHLY'
    short_name = 'WorldClim'

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


@register
class ChirpsDaily(ImageCollection):
    INFO = """ID: {id}
short name: {short_name}"""
    id = "UCSB-CHG/CHIRPS/DAILY"
    short_name = "ChirpsDaily"

    prec = OpticalBand('precipitation', 'precipitation', units='mm', scale=1,
                       resolution=5000, precision=Precisions.float,
                       description='Precipitation')
    start_date = '1981-01-01'

    bands = [prec]

    def _filter_year(self, year):
        start = ee.Date.fromYMD(year, 1, 1)
        end = ee.Date.fromYMD(year, 12, 31)
        collection = self.collection(date=(start, end))
        return collection

    def annualPrecipitation(self, year):
        """ Annual precipitation for the given year """
        collection = self._filter_year(year)
        reduced = collection.select(self.prec.name).reduce(ee.Reducer.sum())
        return ee.Image(reduced).rename('mm').set('YEAR', year)

    def annualPrecipitationSeries(self, start_year=None, end_year=None):
        today = datetime.date.today()
        this_year = today.year
        if end_year is None:
            end_year = this_year-1
        if start_year is None:
            start_year = 1981

        years = ee.List.sequence(start_year, end_year)
        return ee.ImageCollection.fromImages(
            years.map(lambda y: self.annualPrecipitation(y)))

    def monthlyPrecipitation(self, year, month):
        """ Monthly precipitation """
        collection = self._filter_year(year)
        month = ee.Number(month)
        filter = ee.Filter.calendarRange(month, month, 'month')
        filtered = collection.filter(filter)
        image = filtered.reduce(ee.Reducer.sum()).rename('mm')
        date = ee.Date.fromYMD(year, month, 15)
        return image.set('system:time_start', date.millis())

    def monthlyPrecipitationSeries(self, year, month_property='MONTH'):
        """ Monthly precipitation Series. Return an ImageCollection with one image
        per month. Each image has the month name in its properties (MONTH) """
        months = ee.List.sequence(1, 12)
        name_property = '{}_NAME'.format(month_property)

        rel = ee.Dictionary({
            '1': 'January', '2': 'February', '3': 'March', '4': 'April',
            '5': 'May', '6': 'June', '7': 'July', '8': 'August',
            '9': 'September', '10': 'October', '11': 'November',
            '12': 'December'
        })

        def wrap(month):
            image = self.monthlyPrecipitation(year, month)
            month_str = month.toInt().format()
            return image.set(month_property, month) \
                        .set(name_property, ee.String(rel.get(month_str)))

        return ee.ImageCollection.fromImages(months.map(wrap))


