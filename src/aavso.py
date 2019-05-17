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
        ...     writer.writerow({
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
                 obstype='CCD', software='pyaavso'):
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
        n__format = "#TYPE={}\n".format(type)
        fp.write(n__format)
        fp.write('#OBSCODE={}\n'.format(observer_code))
        fp.write("#SOFTWARE={}\n".format(software))
        fp.write("#DELIM={}\n".format(delimiter))
        fp.write("#DATE={}\n".format(date_format.upper()))
        fp.write("#OBSTYPE={}\n".format(obstype))
        fp.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")
        self.writer = csv.writer(fp, delimiter=delimiter)

    def writerow(self, observation_data):
        """
        Writes a single observation to the output file.

        If the ``observation_data`` parameter is a dictionary, it is
        converted to a list to keep a consisted field order (as described
        in format specification). Otherwise it is assumed that the data
        is a raw record ready to be written to file.

        :param observation_data: a single observation as a dictionary or list
        """
        if isinstance(observation_data, (list, tuple)):
            row = observation_data
        else:
            row = self.dict_to_row(observation_data)
        self.writer.writerow(row)

    @classmethod
    def dict_to_row(cls, observation_data):
        """
        Takes a dictionary of observation data and converts it to a list
        of fields according to AAVSO visual format specification.

        #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES

        :param cls: current class
        :param observation_data: a single observation as a dictionary
        """
        row = []
        row.append("{:.30}".format(observation_data['name']))
        row.append("{:.16}".format(observation_data['date']))
        row.append("{:.8}".format(observation_data['magnitude']))
        row.append("{:.6}".format(observation_data['magnitude_error']))
        row.append("{:.5}".format(observation_data['filter']))
        row.append("{:.3}".format(observation_data['transformed']))
        row.append("{:.3}".format(observation_data['magnitude_type']))
        row.append("{:.20}".format(observation_data.get('comparison_name', 'na')))
        row.append("{:.8}".format(observation_data.get('comparison_magnitude', 'na')))
        row.append("{:.20}".format(observation_data.get('check_name', 'na')))
        row.append("{:.8}".format(observation_data.get('check_magnitude', 'na')))
        row.append("{:.7}".format(observation_data.get('airmass', 'na')))
        row.append("{:.5}".format(observation_data.get('group', 'na')))
        row.append("{:.20}".format(observation_data.get('chart', 'na')))
        row.append("{:.2000}".format(observation_data.get('notes', 'na')))
        return row
