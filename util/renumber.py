# load csv module
import csv
import sys
import argparse
import logging

def process(inputfile, outputfile, prefix, start):
    # open file for reading
    with open(inputfile) as input, open(outputfile, 'w') as output:

        # read file as csv file
        csvReader = csv.reader(input)
        writer = csv.writer(output, delimiter=',', quotechar='', quoting=csv.QUOTE_NONE)

        counter = int(start)
        # for every row, print the row
        for row in csvReader:
            if len(row) == 0:
                continue
            elif row[0][0] == '#':
                writer.writerow(row)
                continue

            row[0] = prefix + f"{counter:02}"
            writer.writerow(row)
            counter = counter + 1


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s %(name)s: %(levelname)s %(message)s")
    process(sys.argv[1], sys.argv[2],sys.argv[3], sys.argv[4])
