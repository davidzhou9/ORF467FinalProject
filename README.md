# ORF467FinalProject
Final Project of ORF467 Fall'18: Analyzing the ride-sharing potential and empty re-positioning requirements of a Nationwide aTaxi System

## Task 1: Assessing the Rideshare Potential of Each County In Each County, State and Region

Steps:

1.	Download all relevant CSV files for your assigned state. To do this, first find your state’s FIPS code (for example, New Mexico’s is 35). Then go to NationWideTrips’18Kyle and download all the csv files that start with that state’s FIPS code.
2.	Next, store these files in a subdirectory of the MATLAB file. I called this subdirectory “NewMexicoPersonTrips”.
3.	Make any necessary changes to variables in the MATLAB code and run it in sections. Note that for some of the larger states (especially Texas), it will take hours to run.
a.	A good way to estimate runtime is to use the estimated runtimes I left in the comments of each code section and scale it by number of trips you have. I had like ~7 million trips for New Mexico.
4.	For the “pTripsP_{State}.xlsx”, “pTripsC_{State}.xlsx” and “pTripsS_{State}.xlsx” outputted, take the data from the 3 files and combine into one xlsx file called “{State}_pTrips.xlsx” with three sheets for each granularity level. Store in folder called “{State}_pTrip Files”.
5.	For the “WalkP_{State}.xlsx”, “WalkC_{State}.xlsx” and “WalkS_{State}.xlsx” outputted files do the same combining thing into a file called “RegionWalk.xlsx”. Store in folder called “{State}_Walk Files”.
6.	For Part 3, the algorithm will produce a file called “{State}_{StateFIPS}_vehicleTrip”. Store in a folder called “{State}_vehicleTrip Files”.
7.	For Part 4, take the three outputted files and store in folder called “{State}_aTaxiStatistic Files”.
8.	For Part 5, take the one outputted file and store in folder called “{State}_airportTrip”.
