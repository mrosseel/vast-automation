# observations dict: # {JD, (mag, magerr)}
class ComparisonStars:
    def __init__(self, ids, star_descriptions, observations, comp_catalogmags, comp_catalogerr):
        self.ids = ids
        self.star_descriptions = star_descriptions
        self.observations = observations
        self.comp_catalogmags = comp_catalogmags
        self.comp_catalogerr = comp_catalogerr
