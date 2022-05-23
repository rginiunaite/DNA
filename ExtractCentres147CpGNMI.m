%% take centres 147 of the intervals and extract sequences


%fileID = fopen('.csv','r');

% CpG intersectNMI
%Data = load('SkirmantasData/Simulations/coordinatesCpGintersectNMI.mat');   
%matrix = Data.coordinatesCpGintersectNMI;

%neither
%Data = load('SkirmantasData/Simulations/coordinatesneither.mat');     
%matrix = Data.coordinatesneither;

%CpG without NMI
%Data = load('SkirmantasData/Simulations/coordinatesCpGNOTintersectNMI.mat');   
%matrix = Data.coordinatesCpGNOTintersectNMI;

  
%Data = load('SkirmantasData/Simulations/coordinatesNMINOTintersectCpG.mat');  


for ich =1:22 %% run through all chromosomes


filename = sprintf('SkirmantasData/AllChromosomesSimulations/coordinatesCh%iNMI.mat',ich);
%filename = sprintf('../SkirmantasData/AllChromosomesSimulations/coordinatesCh24CpG.mat');

Data = load(filename);

name = sprintf('coordinatesCh%iNMI',ich);
%name = sprintf('coordinatesChYNMI');

matrix = Data.(name);


result=sum(~(all(matrix==0,2))); %% number of non-zero rows


matrixnonzero = matrix (1:result,:);


coordinates147 = zeros(length(matrixnonzero(:,1)),2);

%% select middle 147 positions
for i=1:length(matrixnonzero(:,1))-1

	%interval1 = fixed.Interval(matrixnonzero(i,1),matrixnonzero(i,2));
	if ((matrixnonzero(i,2)-matrixnonzero(i,1))>147)
	
		%% real interval, e.g. CpG
		%leftside = matrixnonzero(i,1) + round(((matrixnonzero(i,2)-matrixnonzero(i,1)))/2) - 73;% centre of interval 
		%rightside = matrixnonzero(i,1) + round(((matrixnonzero(i,2)-matrixnonzero(i,1)))/2) + 74;% centre of interval
		
		% outside intervals, e.g. nonCpG
		leftside = matrixnonzero(i,2) + round(((matrixnonzero(i+1,1)-matrixnonzero(i,2)))/2) - 73;% centre of interval 
		rightside = matrixnonzero(i,2) + round(((matrixnonzero(i+1,1)-matrixnonzero(i,2)))/2) + 74;% centre of interval 
	end
	
	coordinates147(i,1) = leftside;
	coordinates147(i,2) = rightside;

end


name2 = sprintf('SkirmantasData/AllCentresCh37/coordinatesCh%inotNMI147.mat',ich);
%name2 = sprintf('../SkirmantasData/AllCentres/coordinatesCh24notCpG147.mat');
save(convertCharsToStrings(name2),'coordinates147') 


Data = load(name2);                     
DataField = fieldnames(Data);
%name3 = sprintf('../SkirmantasData/AllCentres/coordinatesCh24notCpG147.txt');
name3 = sprintf('SkirmantasData/AllCentresCh37/coordinatesCh%inotNMI147.txt',ich);
dlmwrite(name3, Data.(DataField{1}), 'precision', 16);

end

