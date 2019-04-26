import os
import argparse
import utils
import init_loader
import logging

def run(datadir, fits_dir, filename, extension):
    reference_frame = filename
    reference_frame_index = utils.find_index_of_file(fits_dir, reference_frame, extension)
    lines = [reference_frame, str(reference_frame_index)]
    with open(datadir + 'reference_frame.txt', 'w') as f:
        f.write('\n'.join(lines))
    logging.info("Reference frame: {}, index: {}".format(reference_frame, reference_frame_index))
    assert reference_frame_index == utils.find_index_of_file(fits_dir, reference_frame, extension)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Set the reference frame to be used.')
    parser.add_argument('filename', help='the full path to the reference frame, should be in the dir with all the fits files.')
    parser.add_argument('datadir', help='the root dir for this dataset, this is where the file will be written.')
    args = parser.parse_args()
    fits_dir = os.path.dirname(args.filename)
    fits_dir = utils.add_trailing_slash(fits_dir)
    datadir = utils.add_trailing_slash(args.datadir)
    filename, file_extension = os.path.splitext(args.filename)
    if not os.path.isfile(args.filename):
        logging.info(f"The given file does not exist: {args.filename}")
        exit(1)
    init_loader.meta_init(datadir)
    run(datadir, fits_dir, args.filename, file_extension)

