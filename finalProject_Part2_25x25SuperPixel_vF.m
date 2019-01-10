%% Part 2 of ORF 467 Final Project 

% Read in Region_oFIPS_vehicleTrip data from part 1
data = readtable('California_06_vehicleTrip.csv');

% Sort so corresponding oFIPS are next to eachother
data = sortrows(data, 1);
%% Create Region_oFIPS_vehiclesInUse table

% Total number of vehicle trips
totalTrips = height(data);
%%
vehicleTrips = zeros(totalTrips, 242);

%% Populate oFips column

% Copy oFIPS column into vehicleTrips talbe
vehicleTrips(1:totalTrips, 1) = data.OFIPS;

%% Populate vehicleTripCount column

% Initial oFIPS value
oFIPS_current = data.OFIPS(1);

% Initial Count value
count = 1;
vehicleTrips(1, 2) = count;

% Iterate through each vehicleTrip
for i = 2:totalTrips
    % Get next oFIPS value
    oFIPS_new = data.OFIPS(i);
    % while oFIPS values are the same, increment vehicleCount
    if oFIPS_current == oFIPS_new
        count = count + 1;
        vehicleTrips(i, 2) = count;
    else % oFIPS values are different
        % Reset count value
        count = 1;
        vehicleTrips(i, 2) = count;
        % Reset current oFIPS value
        oFIPS_current = oFIPS_new;
    end
end

%% Populate time columns

% For first trip, since there are no prior counts
% Get oTime integer
oTime = ceil(data.OTime(1)/360);
% Get vehicle Miles
vehicleMiles = data.CartesianVehicleMiles(1);
% Get trip finish time integer
endTime = ceil(data.OTime(1)/360 + 10*vehicleMiles/(30));

% If the oTime is past 24 hours, adjust both times to 
% take place in the current 24 hours
if oTime >= 240
   oTime = oTime - floor(oTime/240)*240;
   endTime = endTime - floor(endTime/240)*240;
elseif endTime >= 240 % End time is past 24 hours
    endTime = endTime - floor(endTime/240)*240; % adjust end time
   % Times before end of day
   for j = (oTime + 2):242
       vehicleTrips(1, j) = 1;
   end
   % Times after end of day (adjusted)
   for j = 3:(endTime + 2)
       vehicleTrips(1, j) = 1;
   end
else % oTime and endTime are within 24 hour timeframe (normal)
    for j = (oTime + 2):(endTime + 2)
        vehicleTrips(1, j) = 1;
    end
end

% For the rest of trips, iterate through all remaining trips
for i = 2:totalTrips
    
    % Copy data from previous row
    vehicleTrips(i, 3:242) = vehicleTrips(i - 1, 3:242);
    
    % Get oTime integer
    oTime = ceil(data.OTime(i)/360);
    
    % Get vehicle Miles
    vehicleMiles = data.CartesianVehicleMiles(i);
    
    % Get trip finish time integer
    endTime = ceil(data.OTime(i)/360 + 10*vehicleMiles/(30));
    
    % Increment all integers in between oTime and totalTime
    % If the oTime is past 24 hours, adjust both times to 
    % take place in the current 24 hours
    if oTime >= 240
        oTime = oTime - floor(oTime/240)*240;
        endTime = endTime - floor(endTime/240)*240;
    elseif endTime >= 240 % End time is past 24 hours
        endTime = endTime - floor(endTime/240)*240; % adjust end time
        % Times before end of day
        for j = (oTime + 2):242
        vehicleTrips(i, j) = vehicleTrips(i, j) + 1;
        end
        % Times after end of day (adjusted)
        for j = 3:(endTime + 2)
            vehicleTrips(i, j) = vehicleTrips(i, j) + 1;
        end
    else % oTime and endTime are within 24 hour timeframe (normal)
        for j = (oTime + 2):(endTime + 2)
            vehicleTrips(i, j) = vehicleTrips(i, j) + 1;
        end
    end
end

%% Find the FIPS Min Fleet Size

% Get unique pixels
uniqueFIPS = unique(vehicleTrips(:, 1));

% Vector that keeps track of min fleet size per FIPS
FIPS_minVec = zeros(height(table(uniqueFIPS)), 1);

% Min Fleet size for lowest FIPS Code
FIPS_minVec(1) = max(max(vehicleTrips(vehicleTrips(:, 1) == uniqueFIPS(1), 3:242)));

% Iterate through rest of FIPS Codes
for i = 2:height(table(uniqueFIPS))
    % Get cumulative trips in ten minute periods for FIPS Code
    cumTrips = max(vehicleTrips(vehicleTrips(:, 1) == uniqueFIPS(i), 3:242));
    % Check that there is more than one row
    if size(cumTrips, 2) == 1
        % If only one row, just take the vector
        cumTrips = vehicleTrips(vehicleTrips(:, 1) == uniqueFIPS(i), 3:242);
    end
    % Get cumulative trips in ten minute periods for previous FIPS Code
    priorTrips = max(vehicleTrips(vehicleTrips(:, 1) == uniqueFIPS(i - 1), 3:242));
    % Check that there is more than one row
    if size(priorTrips, 2) == 1
        % If only one row, just take the vector
        priorTrips = vehicleTrips(vehicleTrips(:, 1) == uniqueFIPS(i - 1), 3:242);
    end
    % Subtract to get trips contained in a FIPS Code
    FIPStrips = cumTrips - priorTrips;
    % Get min fleet size
    FIPS_minVec(i) = max(FIPStrips);
end

% Sum{FIPSMinFleetSize}
sumFipsMinFleet = sum(FIPS_minVec);

%% Get Region Minimum Fleet Size

% Region Min Fleet Size is maximum of last row
RegionMinFleet = max(vehicleTrips(totalTrips, 3:242));

%% Output Region and Sum of Fips Min Fleet Sizes

fprintf('Region Min Fleet Size is: %d \n', RegionMinFleet);
fprintf('Sum of FIPS Min Fleet Sizes is: %d \n', sumFipsMinFleet);

% Percent difference
percentDiff = 100*(sumFipsMinFleet - RegionMinFleet)/sumFipsMinFleet;
fprintf('RegionMinFleetSize is %.2f percent smaller than the sum of FIPS Min Fleet Size \n', percentDiff);

RegionFleetSize = round(1.1*RegionMinFleet);
fprintf('Region Fleet Size is: %d \n', RegionFleetSize);

%% Write vehicles in Use table

vehiclesInUse = array2table(vehicleTrips);

% Change name to reflect State
writetable(vehiclesInUse, 'California_06_vehiclesInUse.csv');

%% Second Part

% Create Super Pixel Vehicle Trip array
superPixTrip = zeros(totalTrips*2, 5);

% Iterate through trips, populating superpixel table
for i = 1:totalTrips
    % -1 for leaving superpixel
    superPixTrip(2*i - 1, 1) = -1;
    % oXPixel Superpixel
    superPixTrip(2*i - 1, 2) = round(data.OXPixel(i)/25); 
    % oYPixel Superpixel
    superPixTrip(2*i - 1, 3) = round(data.OYPixel(i)/25); 
    if ceil(data.OTime(i)/360) >= 240
        superPixTrip(2*i - 1, 4) = floor(data.OTime(i)/360) - floor(floor(data.OTime(i)/360)/240)*240;
    else
        superPixTrip(2*i - 1, 4) = floor(data.OTime(i)/360);
    end
    
    % 1 for arriving at superpixel
    superPixTrip(2*i, 1) = 1;
    % DXPixel Superpixel
    superPixTrip(2*i, 2) = round(data.DXPixel(i)/25); 
    % oYPixel Superpixel
    superPixTrip(2*i, 3) = round(data.DYPixel(i)/25); 
    % DTime, check if past 86,400
    DTime = data.OTime(i) + 3600*data.CartesianVehicleMiles(i)/30;
    if DTime >= 86400
        superPixTrip(2*i, 4) = floor(DTime/360)- floor(floor(DTime/360)/240)*240;
    else
        superPixTrip(2*i, 4) = floor(DTime/360);
    end
end
%% Sort the data

% Sort based on X SuperPixel, Y SuperPixel, and Time
superPixTripSort = sortrows(superPixTrip, [2 3 4]);

%% Add a column to get a running sum of column 1

% Get all Unique Pixels
uniqueSuperPixels = unique(superPixTripSort(:, 2:3), 'rows');
%%
% Get Current Pixel
currentSuperPixel = uniqueSuperPixels(1, :);
% Set Running count
superPixTripSort(1, 5) = superPixTripSort(1, 1);

% Iterate through remaining pixels
for i = 2:totalTrips*2
    % New Super Pixel
    newSuperPixel = superPixTripSort(i, 2:3);
    % Next Super Pixel is the same
    if isequal(currentSuperPixel, newSuperPixel)
        % Change Count
        count = superPixTripSort(i - 1, 5) + superPixTripSort(i, 1);
        superPixTripSort(i, 5) = count;
    else % Next Super Pixel is different
        % Reset count to 1
        count = superPixTripSort(i, 1);
        superPixTripSort(i, 5) = count;
    end
    % Set new Current Super Pixel
    currentSuperPixel = newSuperPixel;
end

%% Get min values of the Count for each superPixel

% Array to keep track of min values per superpixel
minCountSuperPixel = zeros(height(table(uniqueSuperPixels)), 1);

% Iterate through each unique superpixel
for i = 1:height(table(uniqueSuperPixels))
    % Get SuperPixel Coordinates
    superXPixel = uniqueSuperPixels(i, 1);
    superYPixel = uniqueSuperPixels(i, 2);
    
    % Find Min value of sorted SuperPixel array corresponding
    % to these superpixel coordinates
    
    % Extract SuperPixel Count values corresponding to these superpixels
    superPixelCounts = ...
        superPixTripSort(superPixTripSort(:, 2) == superXPixel & superPixTripSort(:, 3) == superYPixel, 5);
    
    % Find min count
    minCountSuperPixel(i) = min(superPixelCounts);
end

%% save for part 3
minCounts = array2table(minCountSuperPixel);
writetable(minCounts, 'CaliforniaMinCounts.csv');
%% Find how many vehicles we need at midnight to not "stock-out"

% Sum the minimum values across superPixels
% Only use the values that are negative, if the min value is 1,
% that means the superPixel never has a deficit and therefore needs
% zero vehicles at the start of the day.

sumMin = -sum(minCountSuperPixel(minCountSuperPixel < 0));

fprintf(['Number of vehicles we would need at midnight to \n' ...
    'not stock out an any time during the day, assuming \n' ...
    'no repositioning during the day at SuperPixel level: %d \n'], sumMin);

%% Output table to a csv file

SuperPixelVehicleTripOorD = array2table(superPixTripSort);

% Change name to reflect State
writetable(SuperPixelVehicleTripOorD, 'SuperPixel2SuperPixelVehicleTripOorD_California.csv');