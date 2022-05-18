# coding=utf-8
from .helpers import getCommonBands
from ._version import __version__

ALL = list()


def register(cls):
    """ Register the class in ALL """
    classes = [c.__name__ for c in ALL]
    if cls.__name__ in classes:
        raise ValueError('Class {} already registered'.format(cls.__name__))
    shorts = [c.short_name for c in ALL]
    if cls.short_name in shorts:
        raise ValueError('Short name {} already registered'.format(cls.short_name))
    ids = [c.id for c in ALL]
    if cls.id in ids:
        raise ValueError('ID {} already registered'.format(cls.id))

    ALL.append(cls)
    return cls


from . import landsat_c1, landsat, sentinel, modis, climate


def fromId(id):
    """ Create a collection from a parsed ID """
    instances = [instance for instance in ALL if instance.id==id]
    if len(instances)>0:
        dataset = instances[0]
        return dataset()
    else:
        raise ValueError('ID {} not in {}'.format(id, [cls.id for cls in ALL]))


def fromShortName(short_name):
    """ Create a collection from a short name """
    instances = [instance for instance in ALL if instance.short_name==short_name]
    if len(instances)>0:
        dataset = instances[0]
        return dataset()
    else:
        raise ValueError('Name {} not in {}'.format(short_name, [cls.short_name for cls in ALL]))


def options():
    """ Print all available options """
    for ds in ALL:
        print('{} ({})'.format(ds.id, ds.short_name))
