import pandas as pd

################################
year = "2017-18"
dt = pd.read_csv("../sampled_parameter_1000_set_year_"+year+"_10May2019.csv")
df = {}
file_num =500
for num in xrange(file_num+1):
	filename = "results_calibrated_parameters_year_"+year+"num_"+str(num)+".csv"
	try:
		df[num] = pd.read_csv(filename)
		#print ("reading file number ="), num
	except: continue
	

#combine data-frames
df2 = df[0].copy()
for num in xrange(1,file_num+1):
	try: df2 = df2.append(df[num])
	except: continue

	

##check for missing iters
iter_list = list(df2["iter"].unique())
missing_list = []
for num in xrange(1000):
	if num not in iter_list:
		#print ("missing =="), num
		missing_list.append(num)
		


df2.drop_duplicates(subset='iter', inplace=True)
print ("missing list"), year, missing_list

###################################
df2.drop(['incidence(millions)'], axis=1, inplace=True)
dt.drop(['vac_eff_hospitalization', 'vac_eff_mortality'], axis=1, inplace=True)
merged_df = df2.merge(dt, how = 'inner', on = ['iter'])
merged_df.drop(['year_y'], axis=1, inplace=True)
merged_df.rename(columns={'year_x':'year'}, inplace=True)
################################
#to drop rows
print merged_df.shape
merged_df.dropna(inplace=True)
#to drop rows
print merged_df.shape
merged_df = merged_df[:1000]
print merged_df.shape
#####################################
merged_df.to_csv("2.results_calibrated_parameters_year_"+year+"_COMBINED_June10_2019.csv")