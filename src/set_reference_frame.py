import init
import argparse
import do_calibration


def run(filename, extension):
    reference_frame = filename
    reference_frame_index = do_calibration.find_index_of_file(init.fitsdir, reference_frame, extension)
    lines = [reference_frame, str(reference_frame_index)]
    with open(init.basedir + 'reference_frame.txt', 'w') as f:
        f.write('\n'.join(lines))
    print("Reference frame: {}, index: {}".format(reference_frame, reference_frame_index))
    assert reference_frame_index == do_calibration.find_index_of_file(init.fitsdir, reference_frame, extension)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Set the reference frame to be used.')
    parser.add_argument('filename', help='the filename of the reference frame')
    parser.add_argument('-e', '--extension', help='the wildcard to select all files in the fits directory. Default is *.fits')
    args = parser.parse_args()
    extension = '*.fit'
    if args.extension:
        extension = args.extension
    print(args)
    run(args.filename, extension)
