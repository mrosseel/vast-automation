from typing import List
import numpy as np
from star_description import StarDescription


# observations dict: # {JD, (mag, magerr)}
class ComparisonStars:
    def __init__(self, ids, star_descriptions: List[StarDescription], observations, comp_catalogmags, comp_catalogerr):
        self.ids = ids
        self.star_descriptions = star_descriptions
        self.observations = observations
        self.comp_catalogmags = comp_catalogmags
        self.comp_catalogerr = comp_catalogerr

    # return a subset of this ComparisonStars object
    def get_filtered_comparison_stars(self, ids: List[int]):
        mask = np.isin(self.ids, ids)
        return ComparisonStars(np.array(self.ids)[mask], np.array(self.star_descriptions)[mask],
                               np.array(self.observations)[mask],
                               np.array(self.comp_catalogmags)[mask], np.array(self.comp_catalogerr)[mask])


