from typing import List, Dict, Tuple
import numpy as np
from star_description import StarDescription


# observations dict: # {JD, (mag, magerr)}
class ComparisonStars:
    def __init__(self, ids, star_descriptions: List[StarDescription], observations: Dict[str, Tuple[float, float]],
                 comp_catalogmags, comp_catalogerr):
        self.ids = ids
        # one StarDescription per comparison star
        self.star_descriptions = star_descriptions
        # one array of observations per comparison star
        self.observations = observations
        # one catalog magnitude per comparison star
        self.comp_catalogmags = comp_catalogmags
        # one catalog error per comparison star
        self.comp_catalogerr = comp_catalogerr


    # return a subset of this ComparisonStars object


    def get_filtered_comparison_stars(self, ids: List[int]):
        mask = np.isin(self.ids, ids)
        return ComparisonStars(np.array(self.ids)[mask], np.array(self.star_descriptions)[mask],
                               np.array(self.observations)[mask],
                               np.array(self.comp_catalogmags)[mask], np.array(self.comp_catalogerr)[mask])
