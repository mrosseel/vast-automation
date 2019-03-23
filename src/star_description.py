class StarDescription:
    def __init__(self, local_id=None, aavso_id=None, coords=None, vmag=None, e_vmag=None, match=None, upsilon=None, label=None,
                 xpos=None, ypos=None):
        self.local_id = local_id
        self.aavso_id = aavso_id
        self.coords = coords
        self.vmag = vmag
        self.e_vmag = e_vmag
        self._match = match if match != None else []
        self._upsilon = upsilon
        self.xpos = xpos
        self.ypos = ypos
        self.label = self.local_id if label is None else label

    @property
    def upsilon(self):
        return self._upsilon

    @upsilon.setter
    def upsilon(self, val):
        vartype, probability, flag, period = val
        self._upsilon = {'vartype': vartype, 'probability': probability, 'flag': flag, 'period': period}

    @property
    def match(self):
        return self._match

    # val is a CatalogMatch object
    @match.setter
    def match(self, val):
        self._match.append(val)

    # extract matching strings from star_descr
    def get_match_string(self, catalog, strict=True):
        name = None
        separation = None

        if self.match is not None:
            catalog_match_list = [x for x in self.match if x.name_of_catalog == catalog]
            if len(catalog_match_list) != 1:
                if strict:
                    raise AssertionError("star_description.py: Searching for {} in {}, received {} matches, expected 1"
                                         .format(catalog, self, len(catalog_match_list)))
            else:
                name = catalog_match_list[0].catalog_id
                separation = catalog_match_list[0].separation

        return name, separation

    def __repr__(self):
        return "StarDescription({0},{1},{2},{3},{4},{5})".format(
            self.local_id, self.aavso_id, self.coords, self.vmag, self._match, self._upsilon)

    def __str__(self):
        return "local_id: {0}, aavso_id: {1}, coords: {2}, vmag: {3}".format(
            self.local_id, self.aavso_id, self.coords, self.vmag, self._match, self._upsilon)

class CatalogMatch():
    def __init__(self, name_of_catalog=None, catalog_id=None, name=None, coords=None, separation=-1):
        self.name_of_catalog = name_of_catalog
        self.catalog_id = catalog_id
        self.name = name
        self.coords = coords
        self.separation = separation

    def __repr__(self):
        return f'Catalog:{self.name_of_catalog}, CatalogId:{self.catalog_id}, Name:{self.name}, Coords:{self.coords}, Separation:{self.separation}'

    def __str__(self):
        return f'Catalog:{self.name_of_catalog}, CatalogId:{self.catalog_id}, Name:{self.name}, Coords:{self.coords}, Separation:{self.separation}'


# extract upsilon strings from star_descr
# should be migrated to an upsilon catalog_match. Problem is format is quite different, could use extra_dict
def get_upsilon_string(star_description):
    upsilon = star_description.upsilon
    if not upsilon == None:
        upsilon_text = "\nVar: prob={0:.2f}({1}),type={2}".format(upsilon['probability'], upsilon['flag'],
                                                                  upsilon['vartype'])
    else:
        upsilon_text = ''
    return upsilon_text
