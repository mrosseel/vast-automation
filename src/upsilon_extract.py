from init_loader import init, settings
import logging

def getResults(number):
    with open(settings.basedir+'upsilon_output.txt') as f:
        result = []
        for i in range(0,number):
           line = f.readline().split(',')
           result.append(line)
        del result[0]
        return result

def getCopyPng(number):
    logging.info(getFileList(number, 'cp', '.png'))

def getEogPng(number):
    logging.info(getFileList(number, 'eog', '.png'))

def getCopyCurves(number):
    logging.info(getFileList(number, 'cp', '.txt', fileprefix='curve_'))

def getFileList(number, prefix, suffix, fileprefix=''):
    input = getResults(number)
    result = prefix + " "
    for star in input:
        result = result + fileprefix + str(star[0]).zfill(5) + suffix + " "
    return result

def getArray(number):
    input = getResults(number)
    result = []
    for star in input:
        result.append(int(star[0]))
    logging.info(result)


num = 100
getCopyPng(num)
getEogPng(num)
getCopyCurves(num)
getArray(num)
