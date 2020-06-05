# coding=utf-8

ALL = list()


def register(cls):
    """ Register the class in ALL """
    classes = [c.__name__ for c in ALL]
    if cls.__name__ not in classes:
        ALL.append(cls)
    else:
        raise ValueError('Class {} already registered'.format(cls.__name__))
    return cls
