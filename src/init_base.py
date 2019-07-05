#### DO NOT TOUCH ####
# standard directories

#### DO NOT TOUCH ####

# various dir settings
class Settings:
    vsx_catalog_name = "vsx_catalog.bin"
    def __init__(self, datadir, reference_file='new-image.fits'):
        self.datadir = datadir
        self.codedir = './'
        self.basedir = self.datadir
        self.reference_file = reference_file
        self.reference_header = self.basedir + self.reference_file
        self.vsxcatalogdir = self.codedir + Settings.vsx_catalog_name
        self.testdir = self.basedir + 'search/'
        self.conf_phot=self.basedir+'muniphot.conf'
        self.conf_match=self.basedir+'match.conf'
        self.fitsdir = self.basedir + "fits/"
        self.convfitsdir = self.basedir + "converted_fits/"
        self.photometrydir = self.basedir + "photometry/"
        self.matchedphotometrydir = self.basedir + "matched_photometry/"

        # results
        self.resultdir = self.basedir + "results/"
        self.lightcurvedir = self.resultdir + "lightcurves/"
        self.posdir = self.resultdir + "positions/"
        self.worldposdir = self.resultdir + "world_positions/"
        self.chartsdir = self.resultdir + "charts/"
        self.phasedir = self.resultdir + "phasediagrams/"
        self.aavsoreportsdir = self.resultdir + "aavso_reports/"
        self.fieldchartsdirs = self.resultdir + "field_charts/"
