%% upload sequences and extract positions of different sections


%Data = load('SkirmantasData/Simulations/coordinatesCpGNOTintersectNMI147.mat');   
%Data = load('SkirmantasData/Simulations/coordinatesNMINOTintersectCpG147.mat');   
%Data = load('SkirmantasData/Simulations/coordinatesCpGintersectNMI147.mat');   


for ich =1:24 %% run through all chromosomes


filename = sprintf('SkirmantasData/AllCentresCh37/coordinatesCh%inotNMI147.mat',ich);
%filename = sprintf('SkirmantasData/AllChromosomesSimulations/coordinatesChYCpG.mat');

Data = load(filename);


%name = sprintf('coordinatesChYCpG');

matrix = Data.coordinates147;

%% load human genome sequence

fileName = sprintf('SkirmantasData/GRCh37/Chr%i.txt',ich);
%fileName = sprintf('SkirmantasData/sequenceCh24.txt');

% my sample text contains ---> apple, baseball, car, donut, & elephant in single column.
FID = fopen(fileName);
data = textscan(FID,'%s');
fclose(FID);
stringData = string(data{:});


seq = convertStringsToChars(stringData);


allsubsequences = strings([length(matrix(:,1)),1]);

str=strjoin(stringData,''); %% join into one string, before was in shorter segments


for i = 1:length(matrix(:,1))-1
	if (matrix(i,1))<strlength(str)
	subsequences = extractBetween(str,matrix(i,1),matrix(i,2)-1);
	subsequences;	
	allsubsequences(i) = subsequences;
	end
end


name2 = sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%inotNMI147.mat',ich);
%name2 = sprintf('SkirmantasData/AllCentres/allsubsequencesCh24CpG147.mat');

save(convertCharsToStrings(name2),'allsubsequences') 
%Data = load('allsubsequencesnotNMI147.mat');                     
%DataField = fieldnames(Data);


end %% finish running through all chromosomes

