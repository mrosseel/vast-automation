from typing import List, Dict, Tuple
import numpy as np
from star_description import StarDescription


# observations dict: # {JD, (mag, magerr)}
class ComparisonStars:
    def __init__(self, ids, star_descriptions: List[StarDescription], observations: Dict[str, Tuple[float, float]],
                 comp_catalogmags, comp_catalogerr):
        self.ids = ids
        # one StarDescription per comparison star: [index_of_comp_star] = StarDescription
        self.star_descriptions = star_descriptions
        # observations per comparison star: [index_of_comp_star][JD] = Tuple(mag, err)
        self.observations = observations
        # one catalog magnitude per comparison star: [index_of_comp_star] = catalog_mag
        self.comp_catalogmags = comp_catalogmags
        # one catalog error per comparison star: [index_of_comp_star] = catalog_err
        self.comp_catalogerr = comp_catalogerr


    # return a subset of this ComparisonStars object


    def get_filtered_comparison_stars(self, ids: List[int]):
        mask = np.isin(self.ids, ids)
        return ComparisonStars(np.array(self.ids)[mask], np.array(self.star_descriptions)[mask],
                               np.array(self.observations)[mask],
                               np.array(self.comp_catalogmags)[mask], np.array(self.comp_catalogerr)[mask])

    def get_brightest_comparison_star_index(self):
        return np.argmin(self.comp_catalogmags)

    def get_star_id_index(self, star_id):
        return np.where(self.ids == star_id)[0][0]

    def __str__(self):
        return f'ComparisonStars class: ids={self.ids}, #sds={len(self.star_descriptions)}, ' \
               f'#observations={len(self.observations)}, #catalogmags={len(self.comp_catalogmags)}, ' \
               f'#catalogerr={len(self.comp_catalogerr)}.'
