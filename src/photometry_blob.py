# collect = np.full([len(apertures), detectablestars, nrfiles, 2],np.inf, dtype=float)
# fwhm = np.full([nrfiles, 3], np.inf, dtype=float)
# jd = np.full([nrfiles], np.inf, dtype=float)
class PhotometryBlob:
    def __init__(
        self,
        jd=None,
        fwhm=None,
        collect=None,
        apertures=None,
        stddevs=None,
        counts=None,
    ):
        self.jd = jd
        self.fwhm = fwhm
        self.collect = collect
        self.apertures = apertures
        self.stddevs = stddevs
        self.counts = counts

    # def __repr__(self):
    #     return "StarDescription({0},{1},{2},{3},{4})".format(
    #         self.local_id, self.aavso_id, self.coords, self.vmag, self._match)
    #
    # def __str__(self):
    #     return "local_id: {0}, aavso_id: {1}, coords: {2}, vmag: {3}, nr matches: {4}".format(
    #         self.local_id, self.aavso_id, self.coords, self.vmag, len(self._match))
