from typing import List
from astropy.coordinates import SkyCoord
from star_description import StarMetaData


class UpsilonData(StarMetaData):
    def __init__(
        self, key="UPSILON", var_type=None, probability=None, flag=None, period=None
    ):
        super().__init__(key)
        self.var_type = var_type
        self.probability = probability
        self.flag = flag
        self.period = period

    def get_upsilon_string(self):
        # extract upsilon strings from star_descr
        # might have to move outside of UpsilonMatch
        return "\nVar: prob={0:.2f}({1}),type={2}".format(
            self.probability, self.flag, self.var_type
        )

    def __repr__(self):
        return (
            f"Key:{self.key}, Var Type:{self.var_type}, Probability:{self.probability},"
            f" flag:{self.flag}, Period:{self.period}"
        )

    def __str__(self):
        return (
            f"Key:{self.key}, Var Type:{self.var_type}, Probability:{self.probability},"
            f" flag:{self.flag}, Period:{self.period}"
        )


class CompStarData(StarMetaData):
    def __init__(self, compstar_ids: List[int], key="COMPSTARS", extra_id=-1):
        super().__init__(key)
        self.compstar_ids = compstar_ids
        # extra star which is also constant, and can be used as K star for AAVSO ensemble calculations
        self.extra_id = extra_id

    def __repr__(self):
        return f"Key:{self.key}, Compstars: {self.compstar_ids}"

    def __str__(self):
        return f"Key:{self.key}, Compstars: {self.compstar_ids}"


class SelectedFileData(StarMetaData):
    def __init__(self, key="SELECTEDTAG"):
        super().__init__(key)


class SiteData(StarMetaData):
    def __init__(
        self,
        minmax: str = None,
        vsx_var_flag=None,
        separation: float = None,
        var_min=None,
        var_max=None,
        var_type: str = None,
        our_name: str = None,
        period: float = None,
        period_err: str = None,
        source: str = None,
        epoch: str = None,
        comments: str = None,
        key="SITE",
    ):
        super().__init__(key)
        self.minmax = self._strip_if_string(minmax)
        self.vsx_var_flag = self._strip_if_string(vsx_var_flag)
        self.separation = separation
        self.var_min = var_min
        self.var_max = var_max
        self.var_type = self._strip_if_string(var_type)
        self.our_name = self._strip_if_string(our_name)
        self.period = period
        self.period_err = period_err
        self.source = source
        self.epoch = epoch
        self.comments = comments

    @staticmethod
    def _strip_if_string(arg, is_nan=True):
        return arg.strip() if (arg is not None and isinstance(arg, str)) else arg

    def __repr__(self):
        return (
            f"Key:{self.key}, {self.minmax} {self.var_min}-{self.var_max} "
            f"{self.var_type} {self.our_name} {self.period} {self.period_err} {self.epoch}"
        )

    def __str__(self):
        return self.__repr__()


class CatalogData(StarMetaData):
    def __init__(
        self,
        key=None,
        catalog_id=None,
        name=None,
        coords: SkyCoord = None,
        separation=-1,
        vmag=None,
        vmag_err=None,
        extradata=None,
    ):
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
        # the visual magnitude
        self.vmag = vmag
        # the visual magnitude error
        self.vmag_err = vmag_err
        # optional extra info on the catalog object
        self.extradata = extradata

    def get_name_and_separation(self):
        return self.name, self.separation

    def __repr__(self):
        return (
            f"Catalog:{self.key}, CatalogId:{self.catalog_id}, Name:{self.name}, "
            f"Coords:{self.coords}, Separation:{self.separation}"
        )

    def __str__(self):
        return (
            f"Catalog:{self.key}, CatalogId:{self.catalog_id}, Name:{self.name}, "
            f"Coords:{self.coords}, Separation:{self.separation}"
        )


#
# class DataTypes:
#     upsilondata = UpsilonData.key
#     sitedata = SiteData.key
#     compstardata = CompStarData.key
#     selectedfiledata = SelectedFileData.key
