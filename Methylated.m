%% change methylated CpG, to MpN
% find all CG and change it to MN

%data = load('Nucleosomes.mat');

%fileID = fopen('sequences/not_cpg_island_200mer_seqs.txt','r');
%fileID = fopen('optimisation/data/ChenNucleosomesALLseq.txt','r');
%fileID = fopen('sequences/allsubsequencesnotNMI147.txt','r');
%fileID = fopen('sequences/MouseSeq.txt','r');
fileID = fopen('sequences/Sequences1000CpG0Length147.txt','r');
%fileID = fopen('sequences/allsubsequencesCpGNOTintersectNMI147.txt','r');
D = textscan(fileID,'%s');
fclose(fileID);

a = D{1,1};

%filename = 'sequences/mpn_notisland_200mer_seqs.txt';  % better to use fullfile(path,name) 
filename = 'sequences/Sequences1000MpN0Length147.txt' ;
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

        


%% when from .mat file
% for j = 1:length(data.Seq)
%     
% 
%     seq = data.Seq(j).S; 
%     newStr = strrep(seq,'CG','MN');
%     data.Seq(j).S = newStr;
% end

% Seq = data.Seq;
% 
% save('NucleosomesMethylated.mat','Seq')
