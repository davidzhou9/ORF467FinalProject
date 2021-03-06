# ORF467FinalProject
Final Project of ORF467 Fall'18: Analyzing the ride-sharing potential and empty re-positioning requirements of a Nationwide aTaxi System.

Contributors to this repository's code (in alphabetical order):
- David Zhou (dz4@princeton.edu)
- Jasmine Young (jgyoung@princeton.edu)
- Millian Gehrer (mgehrer@princeton.edu)
- Stewart Stroebel (sps3@princeton.edu)

## Task 0: Download this repo

Steps:

1. Download this repo as a zip file by clicking "Clone or download" > "Download ZIP".
2. Alternatively use Git and clone with the following HTTPS: https://github.com/davidzhou9/ORF467FinalProject.git

## Task 1: Assessing the Rideshare Potential of Each County In Each County, State and Region
#### Contact for this Section: David Zhou

Steps:

NOTE: Make sure you download the correct Part1 MATLAB file. If your state has a FIPS code less than 10, use the one titled finalProject_Part1_FIPSLessThan10_vF.m. Otherwise, you use finalProject_Part1_vF.m.

1.	Download all relevant CSV files for your assigned state. To do this, first find your state’s FIPS code (for example, New Mexico’s is 35). Then go to NationWideTrips’18Kyle and download all the csv files that start with that state’s FIPS code.
2.	Next, store these files in a subdirectory of the MATLAB file. I called this subdirectory “NewMexicoPersonTrips”.
3.	Make any necessary changes to variables in the MATLAB code and run it in sections. Note that for some of the larger states (especially Texas), it will take hours to run.
a.	A good way to estimate runtime is to use the estimated runtimes I left in the comments of each code section and scale it by number of trips you have. I had like ~7 million trips for New Mexico.
4.	For the “pTripsP_{State}.xlsx”, “pTripsC_{State}.xlsx” and “pTripsS_{State}.xlsx” outputted, take the data from the 3 files and combine into one xlsx file called “{State}_pTrips.xlsx” with three sheets for each granularity level. Store in folder called “{State}_pTrip Files”.
5.	For the “WalkP_{State}.xlsx”, “WalkC_{State}.xlsx” and “WalkS_{State}.xlsx” outputted files do the same combining thing into a file called “RegionWalk.xlsx”. Store in folder called “{State}_Walk Files”.
6.	For Part 3, the algorithm will produce a file called “{State}_{StateFIPS}_vehicleTrip”. Store in a folder called “{State}_vehicleTrip Files”.
7.	For Part 4, take the three outputted files and store in folder called “{State}_aTaxiStatistic Files”.
8.	For Part 5, take the one outputted file and store in folder called “{State}_airportTrip”.

## Task 2: Assess Fleet Size Requirements for Your Region
#### Contact for this Section: Stewart Stroebel. Jasmine, David and Millian edited/reviewed the code as well.

Steps:

NOTE: For large states, one may need to increase MATLAB's memory's requirements (See Notes below for how to do so). May also need to run on a computer with at least 16 GB of RAM.

1. Create a new folder for Task 2.
2. Copy over "{State}_{StateFIPS}_vehicleTrip.csv" file from the previous part into this new folder.
3. Change the readtable command on the 4th line such that the file name reflects the name of your state's vehicleTrip file.
5. Make any necessary changes to file ouput names, etc. in the code.
4. Run the code and collect the output.

## Task 3: Assess the Empty Vehicle Reposition Implications
#### Contact for this Section: Millian Gehrer, David Zhou

Steps:

1. Create a new folder for Task 3. 
2. Copy over "SuperPixel2SuperPixelVehicleTripOorD_{State}.csv" and "{State}MinCounts.csv" 
3. Change the filenames in lines 4, 44, 72, 88, 94, and 246 to reflect the name of your state.
4. Change the value of regionFleetSize in line 97 to reflect the (1.1 * minFleetSize) as calculated in Part 2.
5. Run the code in sections, saving the graphs and output text along the way. 
6. Upload all necessary generated files to your group's Dropbox folder.

## Notes
- Memory constraints may become prohibitive on some computers due to RAM constraints. As a work around go to Preferences > Workspace > Increase Maximum Array Size to 10,000 and uncheck "Limit the maximum array size to a percentage of RAM". Also go to Preferences > General > Java Heap Memory > Increase to Maximum. Always run on the laptops with the most RAM.
 - Task 2 has been modified to have 6 minute time blocks instead of the 10 minute time blocks specified in the assignment. This is because repositioning in Task 3 is done every 6 minutes (10x a hour).
- Inputs from each step are dependent on previous steps (except for Task 1) so one may need to run all of this repository's code to generated correct inputs/outputs. In other words, if you have generated files from your code, it may not be compatible with this repo's code.
- We performed the repositioning analysis with superpixels of 25x25 due to space and time constraints.
- Feel free to contribute to this repo as well if you spot errors/bugs!
