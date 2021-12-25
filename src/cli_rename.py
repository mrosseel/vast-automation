import re
import logging
import argparse

def renamer(inputcsv, startcounter, regx):
    p = re.compile(regx)
    counter = startcounter
    with open(inputcsv, 'r') as infile, open(inputcsv+".out", 'w') as outfile:
        lines = infile.readlines()
        for line in lines:
            outline = re.sub(regx, f"RMH-HMB{counter}", line, count=1)
            if p.match(line) is not None:
                counter = counter + 1
            outfile.write(outline)
    print(f"Done, start={startcounter}, end={counter}")

if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description="vast automation to rename the temporary star names in the csv file to their final names")
    parser.add_argument(
        "-i",
        "--inputcsv",
        help="the input csv file",
        nargs="?",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--regex",
        help="regex to detect the temporary star names",
        nargs="?",
        required=False,
    )
    parser.add_argument(
        "-s",
        "--startcounter",
        help="the start counter",
        required=True,
    )
    parser.add_argument(
        "-x", "--verbose", help="Set logging to debug mode", action="store_true"
    )
    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if not args.regex:
        regx = "MJ-NEW[0-9]+"
    else:
        regx = args.regex
    renamer(args.inputcsv, int(args.startcounter), regx=regx)
