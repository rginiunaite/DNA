%% methylate all chromosomes


%% change methylated CpG, to MpN
% find all CG and change it to MN

%data = load('Nucleosomes.mat');

for ich=1:24


name = sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpGNOTintersectNMI147.txt',ich);


fileID = fopen(name,'r');
D = textscan(fileID,'%s');
fclose(fileID);

a = D{1,1};

%filename = 'sequences/mpn_notisland_200mer_seqs.txt';  % better to use fullfile(path,name) 
name2 = sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%iMpNNOTintersectNMI147.txt',ich);
filename = name2;
%filename = 'sequences/MethylatedallsubsequencesCpGNOTintersectNMI147.txt' ;


for j = 1:length(a)
    seq = a{j};
    %seq =seq(28:end-27); % from 200 to 146
    %seq = seq(3:end-2);
    start5 = seq(1:5);
    end5 = seq(length(seq)-4:length(seq));
    middle = seq(6:length(seq)-5);
    
    newStrmiddle = strrep(middle,'CG','MN');
    newStr = append(start5,newStrmiddle,end5);
        
    %% methylated only particular indices NMI: [2594 4975]
    
    %middle1 = seq(6:2593);
    %middle2 = seq(2594:4975);
    %middle3 = seq(4976:length(seq)-5);
 
    %newStrmiddle1 = strrep(middle1,'CG','MN');
    %newStrmiddle3 = strrep(middle3,'CG','MN');
    
   % newStr = append(start5,newStrmiddle1,middle2,newStrmiddle3,end5);
    
    fid = fopen(filename,'a');
    fprintf(fid,'%s\n',newStr);          % Write the char array, interpret newline as new line
    fclose(fid);  
   
end

end

