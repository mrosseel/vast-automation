import argparse
import logging
import utils_sd


def main(cmdline_args):
    sds = utils_sd.construct_star_descriptions(cmdline_args.datadir, list(map(lambda x: int(x), cmdline_args.stars)))
    for star in sds:
        print(star, '\n')


def star_to_dat(star: int):
    return f"out{star:05}.dat"


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description='Print all data on one star')
    parser.add_argument('-d', '--datadir',
                        help="The directory where the data can be found (usually the vast dir)",
                        nargs='?', required=True)
    parser.add_argument('-s', '--stars', help="List the star id's to plot", nargs='+')
    parser.add_argument('-p', '--plot', help="Plot the stars using pgfv")
    args = parser.parse_args()
    main(args)
