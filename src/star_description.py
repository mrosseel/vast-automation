from typing import List
from astropy.coordinates import SkyCoord


class StarMetaData:
    def __init__(self, metadata_id=None):
        # the name of the catalog
        self.metadata_id = metadata_id


    def __repr__(self):
        return f'metadata_id:{self.metadata_id}'


    def __str__(self):
        return self.__repr__()


class StarDescription:
    def __init__(self, local_id: int = None, aavso_id: str = None, coords: SkyCoord = None, vmag: float = None,
                 e_vmag: float = None, metadata: List[StarMetaData] = None,
                 label: str = None, xpos: float = None, ypos: float = None, path: str = None, obs: int = None):
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
        self._metadata = metadata if metadata is not None else []
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
    def metadata(self):
        return self._metadata


    @metadata.setter
    def metadata(self, val: StarMetaData):
        self._metadata.append(val)


    def has_metadata(self, metadata_id: str) -> bool:
        """ does this star have a catalog with a certain name? """
        if self.metadata is not None:
            catalog_match_list = [x for x in self.metadata if x.metadata_id == metadata_id]
            if len(catalog_match_list) >= 1:
                return True
        else:
            return False


    def get_metadata(self, catalog: str, strict=False):
        if self.metadata is None:
            return None
        catalog_match_list = [x for x in self.metadata if x.metadata_id == catalog]
        if len(catalog_match_list) != 1:
            if strict:
                raise AssertionError("star_description.py: Searching for {} in {}, received {} matches, expected 1"
                                     .format(catalog, self, len(catalog_match_list)))
        else:
            return catalog_match_list[0]


    def get_metadata_list(self):
        return [x.metadata_id for x in self.metadata]


    def get_metadata_string(self, catalog, strict=False):
        # extract matching strings from star_descr
        # will give an assertion error if the catalog match is not unique
        result = self.get_metadata(catalog, strict)
        if result is None:
            return None, None
        return result.catalog_id, result.separation


    def __repr__(self):
        return "StarDescription({0},{1},{2},{3},{4},{5})".format(
            self.local_id, self.aavso_id, self.coords, self.vmag, self._metadata, self.path)


    def __str__(self):
        return f"local_id: {self.local_id}, aavso_id: {self.aavso_id}, coords: {self.coords}, xy: {self.xpos}, " \
               f"{self.ypos}, vmag: {self.vmag}, nr matches: {len(self._metadata)}, " \
               f"matches: {self._metadata}, path: {self.path}"
