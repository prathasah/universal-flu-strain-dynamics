#!/usr/bin/python
import sys
sys.path.insert(0, r'../Influenza')
import numpy as np
import Simulation
import Optimization
import random
import subprocess


def optimize_universal_vaccine_distribution(season, proportion_universal_vacDoses, objective, index):
	o = Optimization.optimization(season = season, proportion_universalVaccine_doses = proportion_universal_vacDoses, objective = objective, index = index)
	o.optimize()
	
	PVbest, seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise, outcome = o.optimization_output()
	
	return PVbest, seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise, list(outcome)
	
	
##########################################################################
def start_pool(outcome, mpi_index,mpi_chunksize):
	
	year = mpi_index[0]
	uv_cov = mpi_index[1]
	write_index = mpi_index[2]

	for num in range(mpi_index[3], mpi_index[4]):
		print ("starting optimization"), year, outcome, uv_cov, num
		PVbest, seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise, burden_outcome = optimize_universal_vaccine_distribution(year, uv_cov, objective, num)
		elem1 = [year, num, outcome, total_doses, seasonal_vacDoses, universal_vacDoses] +  PVbest + burden_outcome
		print ("optimization done =="), year, outcome, uv_cov, num
		writer[year][write_index].writerow(elem1)
		myfile[year][write_index].flush()
	
	
#######################################################################
def create_files(outcome, yearlist, start_index, num_files):

	writer = {}
	myfile = {}
	    
	header = ["year", "iter", "outcome" ,"total_doses", "seasonal_vacdoses", "universal_vacdoses", "PV_lowrisk_0.5", "PV_lowrisk_5", "PV_lowrisk_10", "PV_lowrisk_15", "PV_lowrisk_20", "PV_lowrisk_25", "PV_lowrisk_30", "PV_lowrisk_35", "PV_lowrisk_40", "PV_lowrisk_45", "PV_lowrisk_50", "PV_lowrisk_55", "PV_lowrisk_60", "PV_lowrisk_65", "PV_lowrisk_70", "PV_lowrisk_75", "PV_highrisk_0.5", "PV_highrisk_5", "PV_highrisk_10", "PV_highrisk_15", "PV_highrisk_20", "PV_highrisk_25", "PV_highrisk_30", "PV_highrisk_35", "PV_highrisk_40", "PV_highrisk_45", "PV_highrisk_50", "PV_highrisk_55", "PV_highrisk_60", "PV_highrisk_65", "PV_highrisk_70", "PV_highrisk_75", "seasonal_vacDoses_0", "seasonal_vacDoses_0.5", "seasonal_vacDoses_5", "seasonal_vacDoses_10", "seasonal_vacDoses_15", "seasonal_vacDoses_20", "seasonal_vacDoses_25", "seasonal_vacDoses_30", "seasonal_vacDoses_35", "seasonal_vacDoses_40", "seasonal_vacDoses_45", "seasonal_vacDoses_50", "seasonal_vacDoses_55", "seasonal_vacDoses_60", "seasonal_vacDoses_65", "seasonal_vacDoses_70", "seasonal_vacDoses_75", "universal_vacDoses_0", "universal_vacDoses_0.5", "universal_vacDoses_5", "universal_vacDoses_10", "universal_vacDoses_15", "universal_vacDoses_20", "universal_vacDoses_25", "universal_vacDoses_30", "universal_vacDoses_35", "universal_vacDoses_40", "universal_vacDoses_45", "universal_vacDoses_50", "universal_vacDoses_55", "universal_vacDoses_60", "universal_vacDoses_65", "universal_vacDoses_70", "universal_vacDoses_75", "total_vacDoses_0", "total_vacDoses_0.5", "total_vacDoses_5", "total_vacDoses_10", "total_vacDoses_15", "total_vacDoses_20", "total_vacDoses_25", "total_vacDoses_30", "total_vacDoses_35", "total_vacDoses_40", "total_vacDoses_45", "total_vacDoses_50", "total_vacDoses_55", "total_vacDoses_60", "total_vacDoses_65", "total_vacDoses_70", "total_vacDoses_75",
  "outcome_0", "outcome_0.5", "outcome_5", "outcome_10", "outcome_15", "outcome_20", "outcome_25", "outcome_30", "outcome_35", "outcome_40", "outcome_45", "outcome_50", "outcome_55", "outcome_60", "outcome_65", "outcome_70", "outcome_75"]
		  
		  
		
	for year in yearlist:
	    writer[year] = {}
	    myfile[year] = {}
	
	    for num in xrange(start_index, start_index+num_files):
		    myfile[year][num] = open('optimized_'+outcome+'_year_'+year+'num_'+str(num)+'.csv','wb')
		    writer[year][num] = csv.writer(myfile[year][num])
		    writer[year][num].writerow(header)
	
	return myfile, writer
		
		
######################################################################################		
if __name__ == "__main__":
	
	
	yearlist = ['2011-12', '2012-13', '2013-14', '2014-15']
	UV_coverage_list = [0, 0.25, 0.5, 0.75, 1]
	#outcome_list = ['totalInfections', 'totalHospitalizations', 'totalDeaths','totalCost']
	outcome = 'totalInfections'
	rank = mpi4py.MPI.COMM_WORLD.Get_rank()
	size = mpi4py.MPI.COMM_WORLD.Get_size()
	mpi_chunks = 5
	mpi_chunksize = 1000/mpi_chunks
	start_index = 0
	mpi_indexes = [((year, UV_coverage, num, num*mpi_chunksize, num*mpi_chunksize+ mpi_chunksize)) for num in xrange(start_index,start_index+mpi_chunks) for UV_coverage in UV_coverage_list for year in yearlist]
	myfile, writer = create_files(outcome, yearlist, start_index, mpi_chunks)
	
	
	for (index, mpi_index) in enumerate(mpi_indexes):
		if index%size!=rank: continue
		print "Task number %d being done by processor %d of %d" % (index, rank, size)
		start_pool(outcome, mpi_index, mpi_chunksize)



    
    
	
		
	
		
