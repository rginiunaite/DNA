paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');
%seqnum = 20;

%from already generated files
% distances = load('distances100sequences.mat');
% 
% distances = cell2mat(struct2cell(distances)); %;zeros(seqnum,350);

%% select specific groups/organisms
vec = [data.Seq.group];

str = 'Virus';
idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
%idx2 = strcmp('Mouse',{data.Seq.group})
i%dxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth

%did not use, did manually
%datatosave = rmfield(datagroup, setdiff(fieldnames(datagroup), {'S'}));  
%save('datatosave')


%writetable(struct2table(datagroup), 'somefile.txt')
