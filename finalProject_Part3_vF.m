%% PART 3A of 467 Final Project
%% i) create SummarySuperPixelVehicleActivity_Region{} 
% read in SuperPixel2SuperPixelFile
SuperPixelVehicleTripOorD = readtable('SuperPixel2SuperPixelVehicleTripOorD_Nevada.csv');
% sort by xPixel, yPixel
SuperPixelVehicleTripOorD = sortrows(SuperPixelVehicleTripOorD, [2 3]);
% find unique xPixels, yPixels
uniqueSuperPixels = unique(SuperPixelVehicleTripOorD(:,2:3),'rows');
numUniqueSuperPixels = height(uniqueSuperPixels); 
%% sum -1s (aTaxiVehicleDepartures)
departures = SuperPixelVehicleTripOorD...
    (SuperPixelVehicleTripOorD.superPixTripSort1 == -1, 2:3);
aTaxiVehicleDepartures = zeros(numUniqueSuperPixels, 1);
for i = 1:numUniqueSuperPixels
    currentX = uniqueSuperPixels{i,1};
    currentY = uniqueSuperPixels{i,2};
    xRows = departures.superPixTripSort2 == currentX;
    yRows = departures.superPixTripSort3 == currentY;
    aTaxiVehicleDepartures(i) = sum(and(xRows, yRows));
end 
%% sum +1s (aTaxiVehicleMadeEmpty)
arrivals = SuperPixelVehicleTripOorD...
    (SuperPixelVehicleTripOorD.superPixTripSort1 == 1, 2:3);
aTaxiVehicleMadeEmpty = zeros(numUniqueSuperPixels, 1);
for i = 1:numUniqueSuperPixels
    currentX = uniqueSuperPixels{i,1};
    currentY = uniqueSuperPixels{i,2};
    xRows = arrivals.superPixTripSort2 == currentX;
    yRows = arrivals.superPixTripSort3 == currentY;
    aTaxiVehicleMadeEmpty(i) = sum(and(xRows, yRows));
end 
%% create and sort summary table
varNames = {'xSuperPixel', 'ySuperPixel', 'SumaTaxiVehicleDepartures',...
    'SumaTaxiVehicleMadeEmpty'};
summaryTable = table(uniqueSuperPixels.superPixTripSort2, ...
    uniqueSuperPixels.superPixTripSort3, aTaxiVehicleDepartures,...
    aTaxiVehicleMadeEmpty, 'VariableNames', varNames);
sortSummaryTable = sortrows(summaryTable, [-3 -4]);
%% assign sequential integers starting with 1 to the sorted array
for i = 1:numUniqueSuperPixels
    sortSummaryTable.index(i) = i; 
end
%% save sorted array
writetable(sortSummaryTable, 'SortedSummarySuperPixelVehicleActivity_Nevada.csv');
%% a-c) calculate and print requested statistics
% a
originPixels = sortSummaryTable.SumaTaxiVehicleDepartures > 0;
madeEmptyPixels = sortSummaryTable.SumaTaxiVehicleMadeEmpty > 0;
numoriginPixels = sum(originPixels);
nummadeEmptyPixels = sum(madeEmptyPixels);
numOriginandmadeEmpty = sum(and(originPixels, madeEmptyPixels)); 
fprintf('# 25x25 superPixel Originating at least one aTaxi: %d \n', numoriginPixels);
fprintf('# 25x25 superPixel where at least one aTaxi is made empty: %d \n', nummadeEmptyPixels);
fprintf('# 25x25 superPixels where at least one aTaxi is originated and one made empty: %d \n', numOriginandmadeEmpty);
% b 
fprintf('Total aTaxiVehicleDepartures: %d \n', sum(aTaxiVehicleDepartures));
fprintf('Total aTaxiVehicleMadeEmpty: %d \n', sum(aTaxiVehicleMadeEmpty));
% c
superPixel1 = strcat('(', num2str(sortSummaryTable{1,1}), ',', num2str(sortSummaryTable{1,2}), ')');
disp('SuperPixel #1:')
disp(superPixel1);
fprintf('SuperPixel #1 aTaxiVehicleDepartures : %d \n', sortSummaryTable{1,3});
fprintf('Percentage of Region Total : %.2f \n', 100*sortSummaryTable{1,3}/sum(aTaxiVehicleDepartures));
fprintf('SuperPixel #1 aTaxiVehicleMadeEmpty : %d \n', sortSummaryTable{1,4});
fprintf('Percentage of Region Total : %.2f \n', 100*sortSummaryTable{1,4}/sum(aTaxiVehicleMadeEmpty));
%% d) Draw aTaxiVehicleDepartures v SuperPixel # // Superimpose aTaxiVehicleMadeEmpty v SuperPixel #
plot(sortSummaryTable.index, sortSummaryTable.SumaTaxiVehicleDepartures, 'b.');
hold on
plot(sortSummaryTable.index, sortSummaryTable.SumaTaxiVehicleMadeEmpty, 'r.');
xlabel('Super Pixel # (in Descending Order of aTaxi Activity)');
ylabel('# aTaxis');
title('aTaxi Activity at Each SuperPixel in Nevada');
legend('SumVehicleDepartures', 'SumVehicleMadeEmpty');
%% e) create hashtable - entering x,y as string gives you index in uniqueSuperPixels
keys = zeros(numUniqueSuperPixels,1);
for i = 1:numUniqueSuperPixels
    keys(i) = uniqueSuperPixels{i,1} * numUniqueSuperPixels + uniqueSuperPixels{i,2};
end 
values = zeros(numUniqueSuperPixels,1);
for i = 1:numUniqueSuperPixels
    values(i) = i;
end 
hashIndex = containers.Map(keys, values);
%% PART 3B of 467 Final Project
%% i) Sort SuperPixelVehileTripOorD{}by the last column, time to create 
% SortedByTimeSuperPixelVehicleTripOorD_Region{} 
% read in SuperPixel2SuperPixelFile
SuperPixelVehicleTripOorD = readtable('SuperPixel2SuperPixelVehicleTripOorD_Nevada.csv');
% sort by Time
sortbyTime = sortrows(SuperPixelVehicleTripOorD, 4);
%% compute StateOfVehicles@SuperPixel
% read in vector of min counts at each pixel (from 2.2.iii) - set 0 if
% positive
minCounts = readtable('NevadaMinCounts.csv');
minCounts = table2array(minCounts);
% set this from part 2, regionFleetSize = ceil(1.1 * MinFleetSize)
regionFleetSize = 134662;
maxFleetSize = -sum(minCounts(:,1));
% assign regionFleet based on need according to maxFleetSize positions
midnightVehicles = -min(minCounts,0)./maxFleetSize .* regionFleetSize;
% take only integer parts, then add remaining vehicles based on largest
% remainder
integerPortion = floor(midnightVehicles);
remainder = midnightVehicles - integerPortion;
while sum(integerPortion) < regionFleetSize
    [maxValue, idx] = max(remainder); 
    integerPortion(idx) = integerPortion(idx) + 1; 
    remainder(idx) = 0; 
end 
stateofVehiclesMidnight = integerPortion;
%% ii) compute StateOfVehicles@SuperPixel{} for all t in 6min increments
stateofVehicles = zeros(numUniqueSuperPixels, 242);
% make first two columns xSuperPixel and ySuperPixel
stateofVehicles(:,1:2) = table2array(uniqueSuperPixels);
% make third column midnight
stateofVehicles(:,3) = stateofVehiclesMidnight(:,1);
%% convert sortbyTime to array for performance
sbt = table2array(sortbyTime);
sbt(:,4) = ceil(sbt(:,4)/360);
%% read through SortedByTimeSuperPixelVehicleTripOorD_Region{} in increments
% of 6 minutes, compute number of vehicles for each superPixel at that time
for i = 1:239
    vector = find(sbt(:,4) == i);
    for j = min(vector):max(vector)
        change = sbt(j,1);
        value = hashIndex(sbt(j,2)*numUniqueSuperPixels + sbt(j,3));
        stateofVehicles(value, 3 + i) = stateofVehicles(value, 3 + i) + change;
    end 
end
%% compute costs for transportation problem
vec1 = zeros(numUniqueSuperPixels,1);
for i = 1:numUniqueSuperPixels
    vec1(i) = i;
end 
[p,q] = meshgrid(vec1, vec1);
pairs = [p(:) q(:)];
cost = hypot(uniqueSuperPixels{pairs(:,1), 1} - uniqueSuperPixels{pairs(:,2), 1}, ...
              uniqueSuperPixels{pairs(:,1), 2} - uniqueSuperPixels{pairs(:,2),2});
% 25x25 super pixels = 12.5x12.5 miles, so 1 superPixel apart = 12.5 miles apart
cost = 12.5*cost; 
%% compute inequality matrix Aeq
sparseData = zeros(2*numUniqueSuperPixels^2, 3);
%
for i = 1:numUniqueSuperPixels
    in = pairs(:,2) == i;
    inIndex = find(in);
    out = pairs(:,1) == i;
    outIndex = find(out);
    rows = i*ones(2*numUniqueSuperPixels,1);
    
    start = (i-1)*numUniqueSuperPixels * 2 + 1;
    middle = start + numUniqueSuperPixels - 1;
    finish = i * numUniqueSuperPixels * 2;
    
    % 1st col = row, 2nd col = xCoord, 3rd col == yCoord
    sparseData(start:finish, 1) = rows;
    sparseData(start:middle, 2) = inIndex;
    sparseData(middle + 1:finish, 2) = outIndex;
    sparseData(start:middle, 3) = 1;
    sparseData(middle + 1:finish, 3) = -1;
end 
%
Aeq = sparse(sparseData(:,1), sparseData(:,2), sparseData(:,3));
% set flow to self equal to 0
for i = 1:numUniqueSuperPixels
    coord = i + numUniqueSuperPixels * (i-1);
    Aeq(i, coord) = 0;
end
%% create data structures to track repositioning
emptyMovement = table;
repositionedMiles = zeros(240,1);
repositionedVehicles = zeros(240,1);
%% repositioning 
disp('starting repositioning');
neededIncreaseinFleet = 0;
movementIndex = 1;
for i = 1:239
    stateofVehicles(:,3+i) = stateofVehicles(:,3+i) + stateofVehicles(:,2+i);
    supply = max(stateofVehicles(:,3+i), 0);
    demand = abs(min(stateofVehicles(:,3+i), 0));
    beq = supply - demand;
    if min(beq) < 0
        finalsol = zeros(numUniqueSuperPixels, numUniqueSuperPixels);
        if sum(beq) < 0
            addedIn = - sum(beq);
            neededIncreaseinFleet = neededIncreaseinFleet + addedIn;
            for j = 1:addedIn
                [minValue, idx] = min(beq); 
                beq(idx) = beq(idx) + 1;
                stateofVehicles(idx,3:3+i) = stateofVehicles(idx,3:3+i) + 1;
            end 
        end 
        [sol, fval] = intlinprog(cost, numUniqueSuperPixels^2, -Aeq, beq, [], [], zeros(numUniqueSuperPixels^2,1));
        repositionedMiles(i,1) = fval;
        repositionedVehicles(i,1) = sum(sol);
        for j = 1:numUniqueSuperPixels
            for k = 1:numUniqueSuperPixels
                finalsol(j, k) = sol(numUniqueSuperPixels*(j-1) + k); 
                if finalsol(j,k) > 0
                    emptyMovement.numaTaxis(movementIndex) = finalsol(j,k);
                    emptyMovement.xEmptyFromSuperPixel(movementIndex) = uniqueSuperPixels{j,1};
                    emptyMovement.yEmptyFromSuperPixel(movementIndex) = uniqueSuperPixels{j,2};
                    emptyMovement.xEmptyToSuperPixel(movementIndex) = uniqueSuperPixels{k,1};
                    emptyMovement.yEmptyToSuperPixel(movementIndex) = uniqueSuperPixels{k,2}; 
                    emptyMovement.Time(movementIndex) = i*360;
                    movementIndex = movementIndex + 1;
                end 
            end 
        end
        for j = 1:numUniqueSuperPixels
           newtaxis = sum(finalsol(:, j));
           leavingtaxis = sum(finalsol(j, :)); 
           stateofVehicles(j,3+i) = stateofVehicles(j, 3+i) + newtaxis - leavingtaxis;
        end 
    end 
    disp(i);
end 
%% EoD Repositioning
vector = find(sbt(:,4) == 0);
midnightChange = zeros(numUniqueSuperPixels, 1);
for j = min(vector):max(vector)
    change = sbt(j,1);
    value = hashIndex(sbt(j,2)*numUniqueSuperPixels + sbt(j,3));
    midnightChange(value,1) = midnightChange(value,1) + change;
end 
beq = (stateofVehicles(:,242) - stateofVehicles(:,3) + midnightChange(:,1));
[sol, fval] = intlinprog(cost, numUniqueSuperPixels^2, -Aeq, beq, [], [], zeros(numUniqueSuperPixels^2,1));
repositionedMiles(240,1) = fval;
repositionedVehicles(240,1) = sum(sol);
finalsol = zeros(numUniqueSuperPixels, numUniqueSuperPixels);
for j = 1:numUniqueSuperPixels
    for k = 1:numUniqueSuperPixels
        finalsol(j, k) = sol(numUniqueSuperPixels*(j-1) + k);
        if finalsol(j, k) > 0
            emptyMovement.numaTaxis(movementIndex) = finalsol(j,k);
            emptyMovement.xEmptyFromSuperPixel(movementIndex) = uniqueSuperPixels{j,1};
            emptyMovement.yEmptyFromSuperPixel(movementIndex) = uniqueSuperPixels{j,2};
            emptyMovement.xEmptyToSuperPixel(movementIndex) = uniqueSuperPixels{k,1};
            emptyMovement.yEmptyToSuperPixel(movementIndex) = uniqueSuperPixels{k,2};
            emptyMovement.Time(movementIndex) = 0;
            movementIndex = movementIndex + 1;
        end
    end
end 
%% save empty movementsfile
writetable(emptyMovement, 'EmptyMovementNevada.csv');
%% print out stats
fprintf('# aTaxis that need to be repositioned at end of day between all 25x25 superPixels: %d \n', repositionedVehicles(240));
fprintf('EoD Empty  aTaxiCartesianMiles (1.1*MinFleetSize): %d \n', repositionedMiles(240));
fprintf('# repositionings during the day: %d \n', sum(repositionedVehicles));
fprintf('During the day (every 6 minutes) Empty  aTaxiCartesianMiles (1.1*MinFleetSize): %d \n', sum(repositionedMiles));
fprintf('During the day Average Empty  aTaxiCartesianMiles (1.1*MinFleetSize): %d \n', sum(repositionedMiles)/sum(repositionedVehicles));
%% vii) 1.i
figure(1);
plot(linspace(1,240,240)/10, repositionedVehicles);
xlim([1 24]);
xlabel('Time of Day (Hours, 24 = midnight)');
ylabel('# Empty aTaxis Repositioned');
title('# of Empty Taxis Repositioned By Time of Day');
%% vii) 1.ii
clf;
plot(linspace(1,240,240)/10, repositionedMiles/repositionedVehicles);
xlim([1 24]);
xlabel('Time of Day (Hours, 24 = midnight)');
ylabel('# Average Empty aTaxi VehicleMiles');
title('# of Average Empty aTaxi VehicleMiles By Time of Day');
%% vii) 2 - 5
fprintf('Average Number of Times Vehicle is Repositioned During Day: %d \n', sum(repositionedVehicles)/(regionFleetSize+neededIncreaseinFleet));
fprintf('Average Repositioning Distance: %d \n', sum(repositionedMiles)/sum(repositionedVehicles));
fprintf('Increase in Fleet Size Needed: %d \n', neededIncreaseinFleet);
%% cdf of repositionedMiles by ToD
repositionedcdf = zeros(240,1);
total = sum(repositionedMiles(:,1));
repositionedcdf(1,1) = repositionedMiles(1,1) * 100/total;
for i = 2:240
    repositionedcdf(i,1) = repositionedcdf(i-1,1) + repositionedMiles(i,1) * 100/total;
end 
clf;
plot(linspace(1,240,240)/10, repositionedcdf, 'r-');
xlim([1 24]);
xlabel('Time of Day (Hours, 24 = midnight)');
ylabel('% of Total VehicleMiles Repositioned By ToD');
title('CDF of Repositioning Distance By Time of Day');
%% Find EoD Repositioning When Using maxFleetSize
stateofVehicles = zeros(numUniqueSuperPixels, 242);
stateofVehicles(:,1:2) = table2array(uniqueSuperPixels);
stateofVehicles(:,3) = -min(minCounts,0);
sbt = table2array(sortbyTime);
sbt(:,4) = ceil(sbt(:,4)/360);
for i = 1:239
    stateofVehicles(:, 3+i) = stateofVehicles(:,2+i);
    vector = find(sbt(:,4) == i);
    for j = min(vector):max(vector)
        change = sbt(j,1);
        value = hashIndex(sbt(j,2)*numUniqueSuperPixels + sbt(j,3));
        stateofVehicles(value, 3 + i) = stateofVehicles(value, 3 + i) + change;
    end 
end
vector = find(sbt(:,4) == 0);
midnightChange = zeros(numUniqueSuperPixels, 1);
for j = min(vector):max(vector)
    change = sbt(j,1);
    value = hashIndex(sbt(j,2)*numUniqueSuperPixels + sbt(j,3));
    midnightChange(value,1) = midnightChange(value,1) + change;
end 
beq = (stateofVehicles(:,242) - stateofVehicles(:,3) + midnightChange(:,1));
[sol, fval] = intlinprog(cost, numUniqueSuperPixels^2, -Aeq, beq, [], [], zeros(numUniqueSuperPixels^2,1), []);
fprintf('EoD Repositioned Miles (Max Fleet Size): %d \n', fval);