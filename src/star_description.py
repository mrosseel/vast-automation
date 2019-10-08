from astropy.coordinates import SkyCoord

class StarDescription:
    def __init__(self, local_id=None, aavso_id=None, coords: SkyCoord=None, vmag=None, e_vmag=None, match=None, label=None,
                 xpos=None, ypos=None, path=None, obs=None):
        # star id given by munipack
        self.local_id = local_id
        # id used to identify to aavso (could be real variable name or UCAC4?)
        self.aavso_id = aavso_id
        # coords as measured by munipack
        self.coords = coords
        # vmag as provided by one of the matching catalogs
        self.vmag = vmag
        # error on vmag as provided by one of the matching catalogs
        self.e_vmag = e_vmag
        # arrray of catalogs matching this star
        self._match = match if match is not None else []
        # x position on the reference frame
        self.xpos = xpos
        # y position on the reference frame
        self.ypos = ypos
        # label to be used when charting this star
        self.label = self.local_id if label is None else label
        # the path where the star is defined
        self.path = path
        # the number of observations for this star
        self.obs = obs


    @property
    def match(self):
        return self._match


    # val is a CatalogMatch object
    @match.setter
    def match(self, val):
        self._match.append(val)
        print(f"appending, is now {self._match}")


    # does this star has a catalog with a certain name?
    def has_catalog(self, catalog: str) -> bool:
        if self.match is not None:
            catalog_match_list = [x for x in self.match if x.name_of_catalog == catalog]
            if len(catalog_match_list) >= 1:
                return True
        else:
            return False


    def get_catalog(self, catalog: str, strict=False):
        if self.match is not None:
            catalog_match_list = [x for x in self.match if x.name_of_catalog == catalog]
            if len(catalog_match_list) != 1:
                if strict:
                    raise AssertionError("star_description.py: Searching for {} in {}, received {} matches, expected 1"
                                         .format(catalog, self, len(catalog_match_list)))
            else:
                return catalog_match_list[0]
        return None

    def get_catalog_list(self):
        return [x.name_of_catalog for x in self.match]

    # extract matching strings from star_descr
    def get_match_string(self, catalog, strict=False):
        # will give an assertion error if the catalog match is not unique
        result = self.get_catalog(catalog, strict)
        if result is None:
            return None, None
        return result.catalog_id, result.separation


    def __repr__(self):
        return "StarDescription({0},{1},{2},{3},{4},{5})".format(
            self.local_id, self.aavso_id, self.coords, self.vmag, self._match, self.path)


    def __str__(self):
        return f"local_id: {self.local_id}, aavso_id: {self.aavso_id}, coords: {self.coords}, xy: {self.xpos}, " \
               f"{self.ypos}, vmag: {self.vmag}, nr matches: {len(self._match)}, " \
               f"matches: {self._match}, path: {self.path}"


class CatalogMatch():
    def __init__(self, name_of_catalog=None, catalog_id=None, name=None, coords: SkyCoord=None, separation=-1):
        # the name of the catalog
        self.name_of_catalog = name_of_catalog
        # the id in the catalog
        self.catalog_id = catalog_id
        # the name of the object in this catalog
        self.name = name
        # the coords in the catalog
        self.coords = coords
        # the separation between the munipack-found coords and the catalog coords
        self.separation = separation


    def __repr__(self):
        return f'Catalog:{self.name_of_catalog}, CatalogId:{self.catalog_id}, Name:{self.name}, Coords:{self.coords}, Separation:{self.separation}'


    def __str__(self):
        return f'Catalog:{self.name_of_catalog}, CatalogId:{self.catalog_id}, Name:{self.name}, Coords:{self.coords}, Separation:{self.separation}'


class UpsilonMatch():
    def __init__(self, name_of_catalog='Upsilon', var_type=None, probability=None, flag=None, period=None):
        self.name_of_catalog = name_of_catalog
        self.var_type = var_type
        self.probability = probability
        self.flag = flag
        self.period = period


    # extract upsilon strings from star_descr
    # might have to move outside of UpsilonMatch
    def get_upsilon_string(self):
        return "\nVar: prob={0:.2f}({1}),type={2}".format(self.probability, self.flag, self.var_type)


    def __repr__(self):
        return f'Catalog:{self.name_of_catalog}, Var Type:{self.var_type}, Probability:{self.probability}, flag:{self.flag}, Period:{self.period}'


    def __str__(self):
        return f'Catalog:{self.name_of_catalog}, Var Type:{self.var_type}, Probability:{self.probability}, flag:{self.flag}, Period:{self.period}'
