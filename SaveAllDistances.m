%% Save all shapes

paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');
%seqnum = 20;

%from already generated files
% distances = load('distances100sequences.mat');
% 
% distances = cell2mat(struct2cell(distances)); %;zeros(seqnum,350);

%% select specific groups/organisms
vec = [data.Seq.group];

% str = 'Yeast';
% idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
% % idx2 = strcmp('Mouse',{data.Seq.group})
% % idxboth = or(idx,idx2);
% datagroup = data.Seq(idx); %idxboth

% 
% % %% all sequences
str = 'All';
datagroup = data.Seq; %% all sequences

seqnum = length(datagroup);

lenseq =73;

k=1; %count;

distances = zeros(167,400-lenseq);

for j = 1:seqnum

    seq = datagroup(j).S; % from nucleosomes data
    
    
    if length(seq)==400
        
        [shapes, stiff] = constructSeqParms(seq, paramset);

        abs_coord = frames(shapes); % relative to absolute coordinates

        for i=1:400-lenseq 
            distances(k,i) = norm(abs_coord(i).rc -abs_coord(i+lenseq).rc);

        end
        k=k+1;
    end


end

save('AllDistances.mat','distances')
