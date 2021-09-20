%% change methylated CpG, to MpN
% find all CG and change it to MN

%data = load('Nucleosomes.mat');

fileID = fopen('sequences/Sequences10000CpG40.txt','r');

D = textscan(fileID,'%s');
fclose(fileID);

a = D{1,1};

filename = 'Sequences10000MpN40.txt';  % better to use fullfile(path,name) 
        
for j = 1:length(a)
    seq = a{j};
    newStr = strrep(seq,'CG','MN');
    
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