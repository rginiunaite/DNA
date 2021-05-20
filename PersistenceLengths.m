%% check persistence lengths


paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');
%seqnum = 20;

%from already generated files
% distances = load('distances100sequences.mat');
% 
% distances = cell2mat(struct2cell(distances)); %;zeros(seqnum,350);

%% select specific groups/organisms
vec = [data.Seq.group];
% 
str = 'Yeast';
idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
% idx2 = strcmp('Mouse',{data.Seq.group})
% idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth

% 
% % % %% all sequences
% str = 'All';
% datagroup = data.Seq; %% all sequences


%% random sequences

%  str = 'Random';
% for j =1 :50
%     datagroup(j).Seq = [random_dna(400)];
% end
% 
% %load 500 previously generated random sequences
% data500 = load('Random500Sequences400length.mat');
% 
% datagroup = data500.datagroup;

%% for all

seqnum = length(datagroup);

lenseq =73;

freqsequence = 25;

dynamic = zeros(seqnum,round((400-lenseq)/freqsequence));
apparent = zeros(seqnum,round((400-lenseq)/freqsequence));

k =1;
for j = 1:seqnum

    seq = datagroup(j).S; % from nucleosomes data
    %seq = datagroup(j).Seq; % for random sequences

    
    if length(seq)==400
        
        c=1;
        for i=1:freqsequence:400-lenseq  
            seqshort = seq(i:i+lenseq);
            seqcg = cgDNAp(seqshort);        
            montecarlo = cgDNAp_MonteCarlo_rahul(seqcg,2);
            dynamic(k,c) = montecarlo.dynamic_pl;
            apparent(k,c) = montecarlo.apparent_pl;         
            c=c+1;

        end
        
        k=k+1;

    end
    
end


save('YeastFull500Every25.mat','dynamic', 'apparent')

