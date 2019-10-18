from typing import List
import numpy as np


# observations dict: # {JD, (mag, magerr)}
class ComparisonStars:
    def __init__(self, ids, star_descriptions, observations, comp_catalogmags, comp_catalogerr):
        self.ids = ids
        self.star_descriptions = star_descriptions
        self.observations = observations
        self.comp_catalogmags = comp_catalogmags
        self.comp_catalogerr = comp_catalogerr

    def get_filtered_comparison_stars(self, ids: List[int]):
        mask = np.isin(self.ids, ids)
        return ComparisonStars(self.ids[mask], self.star_descriptions[mask], self.observations[mask],
                               self.comp_catalogmags[mask], self.comp_catalogerr[mask])


