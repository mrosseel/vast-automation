import csv


class FormatException(Exception):
    """
    Raised when the data does not conform to AAVSO format specification.
    """


class ExtendedFormatWriter(object):
    """
    A class responsible for writing observation data in AAVSO
    `Visual File Format`_.

    The API here mimics the ``csv`` module in Python standard library.

    To write your observations into the data file, you first need to create
    the writer, passing to it the destination file and your observer code.
    Then call :py:meth:`~VisualFormatWriter.writerow` for every single
    observation, for example:

        >>> with open('data.txt', 'wb') as fp:
        ...     writer = VisualFormatWriter(fp, 'XYZ')
        ...     writer.addrow({
        ...         'name': 'SS CYG',
        ...         'date': '2450702.1234',
        ...         'magnitude': '<11.1',
        ...         'comment_code': '',
        ...         'comp1': '110',
        ...         'comp2': '113',
        ...         'chart': '070613',
        ...         'notes': 'This is a test',
        ...     })

    .. _`Visual File Format`: http://www.aavso.org/aavso-visual-file-format
    """


    def __init__(self, fp, observer_code, *, delimiter=',', date_format='JD', type='EXTENDED',
                 obstype='CCD', software='pyaavso', location=("unknown", "unknown", "unknown")):
        """
        Creates the writer which will write observations into the file-like
        object given in first parameter. The only other required parameter
        is the official AAVSO-assigned observer code.

        #TYPE=Extended
        #OBSCODE=
        #SOFTWARE=
        #DELIM=
        #DATE=
        #OBSTYPE=

        :param fp: file-like object to write observations into
        :param observer_code: AAVSO observer code
        :param delimiter: field delimiter (set as DELIM header)
        :param date_format: observation date format (one of *JD* or *Excel*)
        :param obstype: observation type (*Visual* or *PTG*)
        """
        self.observer_code = observer_code
        self.date_format = date_format
        self.obstype = obstype
        self.fp = fp
        self.data = []
        self.data.append(f"#TYPE={type}")
        self.data.append(f'#OBSCODE={observer_code}')
        self.data.append(f"#SOFTWARE={software}")
        self.data.append(f"#DELIM={delimiter}")
        self.data.append(f"#DATE={date_format.upper()}")
        self.data.append(f"#OBSTYPE={obstype}")
        self.data.append(f"#Observation site Latitude = {location[0]}")
        self.data.append(f"#Observation site Longitude = {location[1]}")
        # self.data.append(f"#Observation site Altitude = {location[2]} m")
        self.data.append("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES")


    def flush(self):
        string_to_write = "\n".join(self.data)
        self.fp.write(string_to_write)
        self.data = []


    def addrow(self, observation_data):
        """
        Writes a single observation to the output file.

        If the ``observation_data`` parameter is a dictionary, it is
        converted to a list to keep a consisted field order (as described
        in format specification). Otherwise it is assumed that the data
        is a raw record ready to be written to file.

        :param observation_data: a single observation as a dictionary or list
        """
        if isinstance(observation_data, str):
            row = observation_data
        else:
            row = self.dict_to_row(observation_data)
        self.data.append(row)


    @classmethod
    def dict_to_row(cls, observation_data):
        """
        Takes a dictionary of observation data and converts it to a list
        of fields according to AAVSO visual format specification.

        #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES

        :param cls: current class
        :param observation_data: a single observation as a dictionary
        """
        return f"{observation_data['name']:.30},{observation_data['date']:.16}," \
               f"{observation_data['magnitude']:.3f},{observation_data['magnitude_error']:.3f}," \
               f"{observation_data['filter']:.5},{observation_data['transformed']:.3}," \
               f"{observation_data['magnitude_type']:.3}," \
               f"{observation_data.get('comparison_name', 'na'):.20}," \
               f"{observation_data.get('comparison_magnitude', 'na'):.2}," \
               f"{observation_data['check_name']:.20}," \
               f"{observation_data['check_magnitude']}," \
               f"{observation_data['airmass']:.2f}," \
               f"{observation_data.get('group', 'na'):.5}," \
               f"{observation_data.get('chart', 'na'):.20}," \
               f"{observation_data.get('notes', 'na'):.2000}"
