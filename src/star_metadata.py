from typing import List
from astropy.coordinates import SkyCoord
from star_description import StarMetaData


class UpsilonData(StarMetaData):
    def __init__(self, key='UPSILON', var_type=None, probability=None, flag=None, period=None):
        super().__init__(key)
        self.var_type = var_type
        self.probability = probability
        self.flag = flag
        self.period = period


    def get_upsilon_string(self):
        # extract upsilon strings from star_descr
        # might have to move outside of UpsilonMatch
        return "\nVar: prob={0:.2f}({1}),type={2}".format(self.probability, self.flag, self.var_type)


    def __repr__(self):
        return f'Key:{self.key}, Var Type:{self.var_type}, Probability:{self.probability},' \
               f' flag:{self.flag}, Period:{self.period}'


    def __str__(self):
        return f'Key:{self.key}, Var Type:{self.var_type}, Probability:{self.probability},' \
               f' flag:{self.flag}, Period:{self.period}'


class CompStarData(StarMetaData):
    def __init__(self, key='COMPSTARS', compstar_ids=List[int]):
        super().__init__(key)
        self.compstar_ids = compstar_ids


    def __repr__(self):
        return f'Key:{self.key}, Compstars: {self.compstar_ids}'


    def __str__(self):
        return f'Key:{self.key}, Compstars: {self.compstar_ids}'


class SelectedStarData(StarMetaData):
    def __init__(self, key='SELECTED'):
        super().__init__(key)


class StarFileData(StarMetaData):
    def __init__(self, local_id: int, minmax: str = None, var_min=None, var_max=None, var_type: str = None,
                 our_name: str = None, period: float = None, period_err: float = None, epoch: float = None,
                 key='STARFILE'):
        super().__init__(key)
        self.local_id = local_id
        self.minmax = self._strip_if_not_none(minmax)
        self.var_min = var_min
        self.var_max = var_max
        self.var_type = self._strip_if_not_none(var_type)
        self.our_name = self._strip_if_not_none(our_name)
        self.period = period
        self.period_err = period_err
        self.epoch = epoch


    @staticmethod
    def _strip_if_not_none(arg):
        return arg.strip() if arg is not None else None


    def __repr__(self):
        return f'Key:{self.key}, {self.local_id} {self.minmax} {self.var_min}-{self.var_max} ' \
               f'{self.var_type} {self.our_name} {self.period} {self.period_err} {self.epoch}'


    def __str__(self):
        return self.__repr__()


class CatalogData(StarMetaData):
    def __init__(self, key=None, catalog_id=None, name=None, coords: SkyCoord = None, separation=-1, extradata=None):
        # the name of the catalog
        super().__init__(key)
        # the id in the catalog
        self.catalog_id = catalog_id
        # the name of the object in this catalog
        self.name = name
        # the coords in the catalog
        self.coords = coords
        # the separation between the vast coords and the catalog coords
        self.separation = separation
        # optional extra info on the catalog object
        self.extradata = extradata


    def get_name_and_separation(self):
        return self.name, self.separation


    def __repr__(self):
        return f'Catalog:{self.key}, CatalogId:{self.catalog_id}, Name:{self.name}, ' \
               f'Coords:{self.coords}, Separation:{self.separation}'


    def __str__(self):
        return f'Catalog:{self.key}, CatalogId:{self.catalog_id}, Name:{self.name}, ' \
               f'Coords:{self.coords}, Separation:{self.separation}'
