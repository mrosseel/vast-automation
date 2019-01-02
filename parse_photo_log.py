import re

file = open("output.txt", 'r')
lines = file.readlines()

# patternfirst = r"(\w+.pht)\n(\w+)"
patternfirst = r".*(phot\w+\.pht)"
regexfirst = re.compile(patternfirst)
patternsecond = r"(\w+) "
regexsecond = re.compile(patternsecond)

testline = "./inputfiles/WWCrA_allflat/converted_fits/kout017938.fts => ./inputfiles/WWCrA_allflat/photometry/phot017938.pht"
testmatch = regexfirst.match(testline)
print(testmatch)

# raise SystemExit
print("length of lines is ", len(lines))

anom = 0
for index, line in enumerate(lines):
    firstmatch = regexfirst.match(line)
    # print(f"{index} - match: {firstmatch}")
    if firstmatch and index < len(lines):
        filename = firstmatch.group(1)
        matchy = regexsecond.match(lines[index+1])
        number = matchy.group(1) if matchy else -1
        if int(number) == -1:
            pass
            # print(lines[index+1])
        if int(number) != 10000 and int(number) != -1:
            anom=anom+1
            print(f"For file {filename} there are {number} matches.")
print(f"Number of anomalies: {anom}")