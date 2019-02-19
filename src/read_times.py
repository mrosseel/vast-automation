import os
import pandas as pd

import init

# putting write times of lightcurve files in a CSV t check there are no delays/deadlocks
the_dir = os.listdir(init.lightcurvedir)
rows_list = []
for entry in the_dir:
	dict1 = {}
        the_time = os.path.getmtime(init.lightcurvedir + entry)
	dict1.update({'name': entry, 'time': the_time})
	rows_list.append(dict1)
df  = pd.DataFrame(rows_list)
df.to_csv('read_times.csv')
