from typing import List, Dict
from astropy.coordinates import SkyCoord
import logging
import traceback


class StarMetaData:
    def __init__(self, key: str = None):
        # the name of the catalog
        self.key = key


    def __repr__(self):
        return f'key:{self.key}'


    def __str__(self):
        return self.__repr__()


class StarDescription:
    def __init__(self, local_id: int = None, aavso_id: str = None, coords: SkyCoord = None, vmag: float = None,
                 e_vmag: float = None, metadata: Dict[str, StarMetaData] = None,
                 label: str = None, xpos: float = None, ypos: float = None, path: str = None, obs: int = None):
        # star id given by munipack
        self.local_id = local_id
        # id used to identify to aavso (could be real variable name or UCAC4?)
        self.aavso_id = aavso_id
        # coords as measured by vast
        self.coords = coords
        # vmag as provided by one of the matching catalogs
        self.vmag = vmag
        # error on vmag as provided by one of the matching catalogs
        self.e_vmag = e_vmag
        # arrray of catalogs matching this star
        self._metadata = metadata if metadata is not None else {}
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
        # where can the result files be found. {"phase": "/data/result/bla.png", "light": "...", "lightpa": "..."}
        self.result = {}


    @property
    def metadata(self):
        return self._metadata


    @metadata.setter
    def metadata(self, val: StarMetaData):
        self.set_metadata(val, strict=False)


    def set_metadata(self, val: StarMetaData, strict=True):
        if strict and val.key in self._metadata:
            error = f"Tried to add metadata with key {val.key} but this is already present."
            print("traceback:", traceback.print_exc())
            logging.error(error)
            raise ValueError(error)
        elif val.key in self._metadata:
            logging.warning(f"Overwriting metadata of star {self.local_id} with key {val.key}, strict is False. "
                            f"{self._metadata}")

        self._metadata[val.key] = val


    def has_metadata(self, key: str) -> bool:
        """ does this star have a catalog with a certain name? """
        return key in self._metadata


    def get_metadata(self, key: str):
        if self._metadata is None or not self.has_metadata(key):
            return None
        return self._metadata[key]


    def get_metadata_list(self):
        return list(self.metadata.keys())


    def __repr__(self):
        return "StarDescription({0},{1},{2},{3},{4},{5})".format(
            self.local_id, self.aavso_id, self.coords, self.vmag, self._metadata, self.path)


    def __str__(self):
        return f"local_id: {self.local_id}, aavso_id: {self.aavso_id}, coords: {self.coords}, xy: {self.xpos}, " \
               f"{self.ypos}, vmag: {self.vmag}, nr matches: {len(self._metadata)}, " \
               f"matches: {self._metadata}, path: {self.path}"
