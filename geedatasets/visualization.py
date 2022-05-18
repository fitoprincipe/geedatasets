# coding=utf-8


class Visualizers:
    """ Proxy class to hold many visualizers """
    def __init__(self, **kwargs):
        self.__all = kwargs
        for name, vis in kwargs.items():
            setattr(self, name, vis)

    def get(self, name):
        return self.__all[name]

class Visualization:
    def __init__(self, name, bands, min, max, palette=None, **kwargs):
        # bands must be 1 or 3
        if len(bands)>3:
            self.bands = bands[0:3]
        elif len(bands)<3:
            self.bands = [bands[0]]
        else:
            self.bands = bands

        self.name = name
        self.min = min
        self.max = max
        self.palette = palette
        self.kwargs = kwargs

    def params(self, renamed=False):
        bands = [b.alias for b in self.bands] if renamed else [b.name for b in self.bands]
        p = dict(bands=bands, min=self.min, max=self.max, **self.kwargs)
        if self.palette:
            p['palette'] = self.palette
        return p

    @classmethod
    def by_factor(cls, name, bands, min_factor, max_factor, palette=None, **kwargs):
        mins = list()
        for band in bands:
            if min_factor:
                mins.append(((1-band.offset)/band.scale) / min_factor)
            else:
                mins.append(0)
        maxs = list()
        for band in bands:
            if max_factor:
                maxs.append(((1-band.offset)/band.scale) / max_factor)
            else:
                maxs.append(0)
        return cls(name, bands, mins, maxs, palette, **kwargs)

    @classmethod
    def trueColor(cls, bands):
        return cls.by_factor('trueColor', bands, 0, 3)

    @classmethod
    def falseColor(cls, bands):
        return cls.by_factor('falseColor', bands, 0, 3)

    @classmethod
    def NSR(cls, bands):
        return cls.by_factor('NSR', bands, 0, 2)