%% check nucleosome scores in CpG islands and NMIs 

clear all


% load coordinates CpG
filename = sprintf('../SkirmantasData/AllChromosomesSimulations/coordinatesCh1CpG.txt');
%
coordinatesCpG = load(filename);

% load coordinates notCpG
%filename = sprintf('../SkirmantasData/AllChromosomesSimulations/coordinatesCh24notCpG.txt');

%coordinatesnotCpG = load(filename);


% load score data
%filename = sprintf('../SkirmantasData/Yazdi2015/H1_dyad_withNOS.bed');

%FID = fopen(filename);
%data = textscan(FID,'%s');
%fclose(FID);
%stringData = string(data{:});
%'energiesCh1CpG147.txt'

%str_ = reshape(stringData, 4, [])';

%ind1 = str_(:,1) == 'chrY';

%chr1data = str_(ind1,:);



%T = readtable('../SkirmantasData/Yazdi2015/H1_dyad_withNOS.csv');
%T = readtable('../SkirmantasData/Harwood2019/GSM3312888_iPS_chr1-22_138-161bp.csv','Range', '1:100000');
T = readtable('../SkirmantasData/Yazdi2015/chr1_yatsi_H1_dyad_withNOS.csv');
index = find(strcmp(T.Var1 , 'chr1'));

indexChstart= index(1);

indexChend =index(end);




k=1;
for i =indexChstart:indexChend
	for j=1:length(coordinatesCpG)-1
		%% Yatzi:2015
		if  T.Var2(i)> coordinatesCpG(j,1) && T.Var3(i)<coordinatesCpG(j,2) %CpG islands
		%if  T.Var2(i)> coordinatesCpG(j,2) && T.Var3(i)<coordinatesCpG(j+1,1) %% not CpG islands	
		storeCpG(k,1) = T.Var4(i);
		
		%% Harwood:2019
		%if  T.Var2(i)> coordinatesCpG(j,1) && T.Var2(i)<coordinatesCpG(j,2) %CpG islands
		%if  T.Var2(i)> coordinatesCpG(j,2) && T.Var2(i)<coordinatesCpG(j+1,1) %% not CpG islands	
		%	storeCpG(k,1) = T.Var3(i);		

			k=k+1;
			break
		end
	end
end


writematrix(storeCpG,'ch1YazdiCpGscores.csv')


