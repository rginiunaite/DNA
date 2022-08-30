%% find correlations between energies and nucleosome position scores from H1_dyad_withNOS.bed (Yazdi et al. (2015))
% read coordinatesCh1CpG147.txt file, look for positions from H1_dyad_withNOS.bed that fall into those intervals and add 0 if I can't find such an interval. Add score if it does fall into an interval. Then do correlation of those values with energies. Only do correlations with score values that are non-zeros.


% load score data
filename = sprintf('../SkirmantasData/Yazdi2015/H1_dyad_withNOS.bed');

FID = fopen(filename);
data = textscan(FID,'%s');
fclose(FID);
stringData = string(data{:});
%'energiesCh1CpG147.txt'

str_ = reshape(stringData, 4, [])';

ind1 = str_(:,1) == 'chrY';

chr1data = str_(ind1,:);


% load coordinates 
filename = sprintf('../SkirmantasData/AllCentres/coordinatesCh24CpG147.txt');

coordinates = load(filename);

% load energies 


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh24CpG147.txt');

energies = load(filename);
k=1;
for i =1:length(chr1data)


	for j=1:length(coordinates)-1
		if  str2num(chr1data(i,2))> coordinates(j,1) && str2num(chr1data(i,2))<coordinates(j,2) %% CpG islands

	
			store(k,1) = chr1data(i,4);
			store(k,2) = energies(j);
			
			k=k+1;
			break
		end
	end




end



% store string file, save and convert to numeric file

%writematrix(store,'chYdata.csv')

%a = load('chYdata.csv;)'

%corrcoef(a)

%'coordinatesCh1CpG147.txt'

%fileName = sprintf('SkirmantasData/GRCh37/Chr%i.txt',ich);
%fileName = sprintf('SkirmantasData/sequenceCh24.txt');


