%% finalProject_Part 1 (Assessing the Rideshare potential...)
% note that all estimated runtimes are for New Mexico (~6.7 million total
% trips) adjust accordingly for other states.
%% Reading in all trip data for a certain state
% takes <2 mins, resizing table is not an issue
% taken from HW8 code

numOfFiles = 76; % enter # of csv files your state has
allTrips = table;
subFile = 1; 
FIPSCode = '08001'; % enter the lowest FIPS code of the CSV files
i = 0;
%%
% this part is for reading in the csv files - make necessary changes for
% state names, etc.
while i < numOfFiles
   fileName = strcat('ColoradoPersonTrips/FinalOriginPixel', FIPSCode, '_', num2str(subFile), '.csv');
   %fileName
   %break
   if isfile(fileName)
       data = readtable(fileName);
       subFile = subFile + 1;
       i = i + 1;
   else
       subFile = 1;
       tempFIPS = str2num(FIPSCode);
       tempFIPS = tempFIPS + 1;
       FIPSCode = strcat('0', num2str(tempFIPS));
       continue;
   end
   allTrips = vertcat(allTrips, data);
end

%% Split the trips into those greater 350 miles and those less
% takes a few seconds

normalTrips = allTrips(allTrips.GCDistance < 350, :);
airplaneTrips = allTrips(allTrips.GCDistance >= 350, :);

%% Subsection 1 - "Running Count of personTrips for each Pixel..."

%% Get lengths of certain vars (< 10 seconds)
normalTripsOXCoord = normalTrips.OXCoord;
normalTripsOYCoord = normalTrips.OYCoord;
allXY = horzcat(normalTripsOXCoord, normalTripsOYCoord);

uniquePixels = unique(allXY, 'rows');
numUniquePixels = length(uniquePixels);
uniqueCounties = unique(normalTrips.OFIPS);
numUniqueCounties = length(uniqueCounties);

pixelTripsTable = table;
countyTripsTable = table;

%% Get running counts for personTrips for each Pixel, County and the overall state

%% create pixelTripsTable
xPixelCoord = zeros(numUniquePixels, 1);
yPixelCoord = zeros(numUniquePixels, 1);
countTrips = zeros(numUniquePixels, 1);
%% This loop takes about 10 mins
for i = 1:numUniquePixels
   currX = uniquePixels(i, 1);
   currY = uniquePixels(i, 2);

   numOfCurrPixel = height(normalTrips(normalTrips.OXCoord == currX & normalTrips.OYCoord == currY, 1));
   xPixelCoord(i, 1) = currX;
   yPixelCoord(i, 1) = currY;
   countTrips(i, 1) = numOfCurrPixel;
   
end
%%
pixelTripsTable.xPixel = xPixelCoord;
pixelTripsTable.yPixel = yPixelCoord;
pixelTripsTable.count = countTrips;

%% create countyTripsTable
% takes < 10 seconds

countyFIPS = zeros(numUniqueCounties, 1);
countCountyTrips = zeros(numUniqueCounties, 1);

for i = 1:numUniqueCounties
   currCounty = uniqueCounties(i, 1);

   numOfCurrCounty = height(normalTrips(normalTrips.OFIPS == currCounty, 1));
   countyFIPS(i, 1) = currCounty;
   countCountyTrips(i, 1) = numOfCurrCounty;
   
end

countyTripsTable.CountyFIPS = countyFIPS;
countyTripsTable.count = countCountyTrips;

%% Create stateTripsTable
stateTripsTable = table;
stateTripsTable.StateFIPS = 35;
stateTripsTable.count = height(normalTrips);


%% Write to Workbook files
% takes < 10 seconds

writetable(pixelTripsTable, 'pTripsP_NewMexico.xlsx');
writetable(countyTripsTable, 'pTripsC_NewMexico.xlsx');
writetable(stateTripsTable, 'pTripsS_NewMexico.xlsx');

%% Subsection 2 - "If GCD < 0.707 it is assigned to walking..."


%% Get all trips with GCD < 0.707 and taxi trips
% Takes < 10 seconds

walkingTrips = normalTrips(normalTrips.GCDistance < 0.707, :);
taxiTrips = normalTrips(normalTrips.GCDistance >= 0.707, :);

%% Get lengths of certain var for determining walking trip counts
% takes < 10 seconds
walkingTripsOXCoord = walkingTrips.OXCoord;
walkingTripsOYCoord = walkingTrips.OYCoord;
allXYWalking = horzcat(walkingTripsOXCoord, walkingTripsOYCoord);

uniquePixelsWalking = unique(allXYWalking, 'rows');
numUniquePixelsWalking = length(uniquePixelsWalking);
uniqueCountiesWalking = unique(walkingTrips.OFIPS);
numUniqueCountiesWalking = length(uniqueCountiesWalking);

pixelTripsTableWalking = table;
countyTripsTableWalking = table;

%% Get running counts for walking personTrips for each Pixel, County and the overall state

%% create pixelTripsTable for Walking
xPixelCoordWalking = zeros(numUniquePixelsWalking, 1);
yPixelCoordWalking = zeros(numUniquePixelsWalking, 1);
countTripsWalking = zeros(numUniquePixelsWalking, 1);
% This loop takes about 30 seconds
for i = 1:numUniquePixelsWalking
   currX = uniquePixelsWalking(i, 1);
   currY = uniquePixelsWalking(i, 2);

   numOfCurrPixel = height(walkingTrips(walkingTrips.OXCoord == currX & walkingTrips.OYCoord == currY, 1));
   xPixelCoordWalking(i, 1) = currX;
   yPixelCoordWalking(i, 1) = currY;
   countTripsWalking(i, 1) = numOfCurrPixel;
   
end
%%
pixelTripsTableWalking.xPixel = xPixelCoordWalking;
pixelTripsTableWalking.yPixel = yPixelCoordWalking;
pixelTripsTableWalking.count = countTripsWalking;

%% create countyTripsTable for Walking
% takes < 10 seconds

countyFIPSWalking = zeros(numUniqueCountiesWalking, 1);
countCountyTripsWalking = zeros(numUniqueCountiesWalking, 1);

for i = 1:numUniqueCountiesWalking
   currCounty = uniqueCountiesWalking(i, 1);

   numOfCurrCounty = height(walkingTrips(walkingTrips.OFIPS == currCounty, 1));
   countyFIPSWalking(i, 1) = currCounty;
   countCountyTripsWalking(i, 1) = numOfCurrCounty;
   
end

countyTripsTableWalking.CountyFIPS = countyFIPSWalking;
countyTripsTableWalking.count = countCountyTripsWalking;

%% Create stateTripsTable for Walking
stateTripsTableWalking = table;
stateTripsTableWalking.StateFIPS = 35;
stateTripsTableWalking.count = height(walkingTrips);

%% Write to Workbook files
% takes < 10 seconds

writetable(pixelTripsTableWalking, 'WalkP_NewMexico.xlsx');
writetable(countyTripsTableWalking, 'WalkC_NewMexico.xlsx');
writetable(stateTripsTableWalking, 'WalkS_NewMexico.xlsx');

%% Subsection 3 - "If 0.707 <= GCD < 350 it is shared aTaxi trip..."

%% Get Destination Zone
% takes ~ 3 mins

destZones = zeros(height(taxiTrips), 1);

for i = 1:height(taxiTrips)
    oXTemp = taxiTrips.OXCoord(i);
    oYTemp = taxiTrips.OYCoord(i);
    
    dXTemp = taxiTrips.DXCoord(i);
    dYTemp = taxiTrips.DYCoord(i);
    
    
    
    deltaX = dXTemp - oXTemp;
    deltaY = dYTemp - oYTemp;
    
    arcTan2Degrees = atan2(deltaY, deltaX) * 180 / pi;
    if arcTan2Degrees < 0
        arcTan2Degrees = arcTan2Degrees + 360;
    end
    
    gcDistCurr = taxiTrips.GCDistance(i);
    
    if gcDistCurr < 3
        destZones(i) = 1 + floor(arcTan2Degrees / 90);
    else
        destZones(i) = 5 + floor(arcTan2Degrees / 45);
    end
end

%% Set destination zones as a column on the taxi table
taxiTrips.DestZone = destZones;

%% Sort the trips by OXCoord, OYCoord, DestZone then ODeparture Time
% takes about <10 seconds
sortedTaxiTrips = sortrows(taxiTrips, [9 10 20 11]);
%% Set up for assigning the customers to the aTaxis
taxiDepartNum = zeros(height(sortedTaxiTrips), 1);
taxiDepartTime = zeros(height(sortedTaxiTrips), 1);
taxiNumCustomers = zeros(height(sortedTaxiTrips), 1);

% set up just the first customer's aTaxi as base case
taxiDepartNum(1,1) = 1;
firstGCDist = sortedTaxiTrips.GCDistance(1);
if firstGCDist < 3.0
    taxiDepartTime(1,1) = sortedTaxiTrips.ODepartureTime(1) + 300;
else
    taxiDepartTime(1,1) = sortedTaxiTrips.ODepartureTime(1) + 450;
end
taxiNumCustomers(1,1) = 1;

%% ASSIGNING CUSTOMERS TO aTAXIS - THIS IS THE HOT LOOP
% takes ~ 4 mins

numCustInCurrTaxi = 1;

for i = 2:height(sortedTaxiTrips)
   if sortedTaxiTrips.OXCoord(i) == sortedTaxiTrips.OXCoord(i - 1) && ...
       sortedTaxiTrips.OYCoord(i) == sortedTaxiTrips.OYCoord(i - 1) && ...
       sortedTaxiTrips.DestZone(i) == sortedTaxiTrips.DestZone(i - 1) && ...
       sortedTaxiTrips.ODepartureTime(i) < taxiDepartTime(i - 1) && ...
       numCustInCurrTaxi < 6
        
        taxiDepartNum(i, 1) = taxiDepartNum(i - 1, 1);
        numCustInCurrTaxi = numCustInCurrTaxi + 1;
        taxiNumCustomers(i, 1) = numCustInCurrTaxi;
        taxiDepartTime(i, 1) = taxiDepartTime(i - 1, 1);
     
   else
        taxiDepartNum(i, 1) = taxiDepartNum(i - 1, 1) + 1;
        numCustInCurrTaxi = 1;
        taxiNumCustomers(i, 1) = numCustInCurrTaxi;
        if sortedTaxiTrips.GCDistance(i) < 3
            taxiDepartTime(i, 1) = sortedTaxiTrips.ODepartureTime(i) + 300;
        else
            taxiDepartTime(i, 1) = sortedTaxiTrips.ODepartureTime(i) + 450;
        end
   end
   
end

%% Set the new columns in the table
sortedTaxiTrips.aTaxiDepartureNum = taxiDepartNum;
sortedTaxiTrips.aTaxiDepartureTime = taxiDepartTime;
sortedTaxiTrips.aTaxiCustomerNum = taxiNumCustomers;

%% Calculate the pixel distance for all these trips
pixelDist = 0.5 * ((sortedTaxiTrips.DXCoord - sortedTaxiTrips.OXCoord).^2 + ...
    (sortedTaxiTrips.DYCoord - sortedTaxiTrips.OYCoord).^2).^0.5;

sortedTaxiTrips.pixelDistance = pixelDist;
%% Get only unique taxi departures with necessary fields
%1st col = OFIPS, 2nd = OXCoord, 3rd = OYCoord, 4th = DXCoord, 5th =
%DYCoord, 6th = aTaxiDepartureNum, 7th = aTaxiDepartureTime, 8th = aTaxiCustomerNum, 
%9th = pixelDistance
matrixTrips = sortedTaxiTrips{:, [6 9 10 17 18 21 22 23 24]};

finalTaxiTrips = cell(max(sortedTaxiTrips.aTaxiDepartureNum), 9);
finalUniqueDepartNums = unique(sortedTaxiTrips.aTaxiDepartureNum);

%% Loop for calculating cartesianPersonTripMiles, cartesianVehicleMiles, etc.
% columns for finalTaxiTrips: 1 = oFIPS, 2 = OXPixel, 3 = OYPixel, 4 = oTime, 5 = DXPixel, 6 =
% DYPixel, 7 = DepartingOccupancy, 8 = SumCartesianPersonTripMiles, 9 =
% CartesianVehicleMiles

% Fun fact: This loop would've initially taken ~ 7 days to run if you used tables.
% The trick is to use cell matrices rather than tables for manipulation (they are
% about 100x faster). This loop takes about ~1 min

indexCurr = 1;

for i = 1:length(finalUniqueDepartNums)
    currDepartNum = finalUniqueDepartNums(i);
    tempSumPersonTripMiles = 0;
    tempSumVehicleMiles = 0;
    prevXCoord = matrixTrips(indexCurr, 2);
    prevYCoord = matrixTrips(indexCurr, 3);
    
    while indexCurr <= height(sortedTaxiTrips) && matrixTrips(indexCurr, 6) == currDepartNum
        tempSumPersonTripMiles = tempSumPersonTripMiles + matrixTrips(indexCurr, 9);
       
        currXCoord = matrixTrips(indexCurr, 4);
        currYCoord = matrixTrips(indexCurr, 5);
        tempSumVehicleMiles = tempSumVehicleMiles + 0.5 * ...
            ((currXCoord - prevXCoord)^2 ...
            + (currYCoord - prevYCoord)^2)^0.5;
        prevXCoord = currXCoord;
        prevYCoord = currYCoord;
        indexCurr = indexCurr + 1;
    end
    
    finalTaxiTrips(i, 8) = num2cell(tempSumPersonTripMiles);
    finalTaxiTrips(i, 9) = num2cell(tempSumVehicleMiles);
   
    
    % we just get the last one of these values since it's the same for all
    % trips in the current taxi
    finalTaxiTrips(i, 1) = num2cell(matrixTrips(indexCurr - 1, 1)); %OFIPS
    finalTaxiTrips(i, 2) = num2cell(matrixTrips(indexCurr - 1, 2)); %OXCoord
    finalTaxiTrips(i, 3) = num2cell(matrixTrips(indexCurr - 1, 3)); %OYCoord
    finalTaxiTrips(i, 4) = num2cell(matrixTrips(indexCurr - 1, 7)); %departure time
    finalTaxiTrips(i, 5) = num2cell(matrixTrips(indexCurr - 1, 4)); %DXCoord
    finalTaxiTrips(i, 6) = num2cell(matrixTrips(indexCurr - 1, 5)); %DYCoord
    finalTaxiTrips(i, 7) = num2cell(matrixTrips(indexCurr - 1, 8)); %aTaxiCustomerNum

end
%% Convert finalTaxiTrips to a table for output
% takes about 15 seconds
finalTaxiTripsTable = cell2table(finalTaxiTrips);
finalTaxiTripsTable.Properties.VariableNames = {'OFIPS' 'OXPixel' 'OYPixel' ...
'OTime' 'DXPixel' 'DYPixel' 'DepartingOccupancy' 'SumCartesianPersonTripMiles' ...
'CartesianVehicleMiles'};

writetable(finalTaxiTripsTable, 'NewMexico_35_vehicleTrip.csv');
%% Subsection 4 - "Create Summary Stats for each Pixel..."


%% Write our Region_oFIPS_vehicleTrip and aTaxiC files
% about ~2 mins
uniqueCountiesTaxiTrips = unique(finalTaxiTripsTable.OFIPS);
numUniqueCountiesTaxiTrips = length(uniqueCountiesTaxiTrips);

% col 1 = CountyFIPS, 2 = personTripCount, 3 = vehicleTripCount, 4 =
% personTripMiles, 5 = vehicleTripMiles
aTaxiCTable = cell2table(cell(numUniqueCountiesTaxiTrips, 5));
aTaxiCTable.Properties.VariableNames = {'CountyFIPS' 'personTripCount' ...
    'vehicleTripCount' 'personTripMiles' 'vehicleTripMiles'};

for i = 1:numUniqueCountiesTaxiTrips
   currCounty = uniqueCountiesTaxiTrips(i, 1);

   currCountyTable = finalTaxiTripsTable(finalTaxiTripsTable.OFIPS == currCounty, :);
   
   % creating the aTaxiC file
   aTaxiCTable.CountyFIPS(i) = num2cell(currCounty);
   aTaxiCTable.personTripCount(i) = num2cell(sum(currCountyTable.DepartingOccupancy));
   aTaxiCTable.vehicleTripCount(i) = num2cell(height(currCountyTable));
   aTaxiCTable.personTripMiles(i) = num2cell(sum(currCountyTable.SumCartesianPersonTripMiles));
   aTaxiCTable.vehicleTripMiles(i) = num2cell(sum(currCountyTable.CartesianVehicleMiles));
   
   % outputting each Region_oFIPS_vehicleTrip file
   % fileName = strcat('NewMexico_', num2str(currCounty), '_vehicleTrip.xlsx');
   % writetable(currCountyTable, fileName);
   
end
%% Output the aTaxiC file
writetable(aTaxiCTable, 'NewMexico_aTaxiC.xlsx');

%% Create the aTaxiP file
% get unique pixels from the finalTaxiTripsTable
taxiTripsOXCoord = finalTaxiTripsTable.OXPixel;
taxiTripsOYCoord = finalTaxiTripsTable.OYPixel;
allXYTaxi = horzcat(taxiTripsOXCoord, taxiTripsOYCoord);

uniqueTaxiPixels = unique(allXYTaxi, 'rows');
numUniqueTaxiPixels = length(uniqueTaxiPixels);
%% Now create the aTaxiP file
% takes ~ 20 seconds

% col 1 = xPixel, 2 = yPixel, 3 = personTripCount, 4 = vehicleTripCount, 5 =
% personTripMiles, 6 = vehicleTripMiles
aTaxiPMatrix = cell(numUniqueTaxiPixels, 6);
matrixIndex = 1;

for i = 1:numUniqueTaxiPixels
   currX = uniqueTaxiPixels(i, 1);
   currY = uniqueTaxiPixels(i, 2);
   aTaxiPMatrix(i, 1) = num2cell(currX);
   aTaxiPMatrix(i, 2) = num2cell(currY);
   
   tempPersonTripCount = 0;
   tempVehicleTripCount = 0;
   tempPersonTripMiles = 0;
   tempVehicleTripMiles = 0;
   
   while matrixIndex <= length(finalTaxiTrips) && finalTaxiTrips{matrixIndex, 2} == currX ...
       && finalTaxiTrips{matrixIndex, 3} == currY
        tempVehicleTripCount = tempVehicleTripCount + 1;
        tempPersonTripCount = tempPersonTripCount + finalTaxiTrips{matrixIndex, 7};
        tempPersonTripMiles = tempPersonTripMiles + finalTaxiTrips{matrixIndex, 8};
        tempVehicleTripMiles = tempVehicleTripMiles + finalTaxiTrips{matrixIndex, 9};
        matrixIndex = matrixIndex + 1;
   end
   
   % creating the aTaxiC file
   aTaxiPMatrix(i, 3) = num2cell(tempPersonTripCount); %personTripCount
   aTaxiPMatrix(i, 4) = num2cell(tempVehicleTripCount); %vehicleTripCount
   aTaxiPMatrix(i, 5) = num2cell(tempPersonTripMiles); %personTripMiles
   aTaxiPMatrix(i, 6) = num2cell(tempVehicleTripMiles); %vehicleTripMiles
end

%% Output the aTaxiP file
aTaxiPTable = cell2table(aTaxiPMatrix);
aTaxiPTable.Properties.VariableNames = {'xPixel' 'yPixel' ...
    'personTripCount' 'vehicleTripCount' 'personTripMiles' 'vehicleTripMiles'};
writetable(aTaxiPTable, 'NewMexico_aTaxiP.xlsx');

%% Create the aTaxiS file
aTaxiSTable = table;
aTaxiSTable.StateFIPS = 35;
aTaxiSTable.personTripCount = sum(finalTaxiTripsTable.DepartingOccupancy);
aTaxiSTable.vehicleTripCount = height(finalTaxiTripsTable);
aTaxiSTable.personTripMiles = sum(finalTaxiTripsTable.SumCartesianPersonTripMiles);
aTaxiSTable.vehicleTripMiles = sum(finalTaxiTripsTable.CartesianVehicleMiles);

writetable(aTaxiSTable, 'NewMexico_aTaxiS.xlsx');

%% Subsection 5 - "New if GCD >= 350 miles"
% this is for the airplane trips
airplaneTripsTable = table;
airplaneTripsTable.oXPixel = airplaneTrips.OXCoord;
airplaneTripsTable.oYPixel = airplaneTrips.OYCoord;
airplaneTripsTable.oTime = airplaneTrips.ODepartureTime;
%% Write the airplaneTripsTable
writetable(airplaneTripsTable, 'airportTrip35.xlsx');