clear all
% Kullbach-Leibler divergences

paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');
%seqnum = 20;

%from already generated files
% distances = load('distances100sequences.mat');
% 
% distances = cell2mat(struct2cell(distances)); %;zeros(seqnum,350);

%% select specific groups/organisms
vec = [data.Seq.group];

str = 'Yeast';
idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
% idx2 = strcmp('Mouse',{data.Seq.group})
% idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth


% % %% all sequences
% str = 'All';
% datagroup = data.Seq; %% all sequences


%% random sequences

%  str = 'Random';
% % for j =1 :50
% %     datagroup(j).Seq = [random_dna(400)];
% % end
% 
% %load 500 previously generated random sequences
% data500 = load('Random500Sequences400length.mat');
% 
% datagroup = data500.datagroup;

% % for all
% 
seqnum = length(datagroup);

lenseq =50;

%numberseq = 400/lenseq; % number of sequences

numberseq = (400-lenseq)/10; %when overlapping with 10 OVERLAPPING


k=1; %count;


AverageMatrix = zeros(numberseq,numberseq);

for j = 1:seqnum

    seq = datagroup(j).S; % from nucleosomes data
    %seq = datagroup(j).Seq; % for random sequences
    

    KullDiv = zeros(numberseq,numberseq);

        if length(seq)==400

    
        for i=1:numberseq
            for j=1:numberseq
                sequencecell = extractBetween(seq,lenseq/5*(i-1)+1,lenseq/5 *i + lenseq); % group into sequence of length n
                sequence = char(sequencecell);

                sequencecell = extractBetween(seq,lenseq/5*(j-1)+1,lenseq/5 *j + lenseq); % group into sequence of length n
                sequence2 = char(sequencecell);

                [shapes1, stiff1] = constructSeqParms(sequence, paramset);

                [shapes2, stiff2] = constructSeqParms(sequence2, paramset);

                KullDiv(i,j) = kl_div(stiff1,shapes1,stiff2,shapes2);

            end
        end
        
        save(['DataKL/KullDiv' num2str(k) '.mat'])
        AverageMatrix = AverageMatrix + KullDiv;
        end
        k=k+1;
end

AverageMatrix = AverageMatrix/10;
save(['DataKL/AverageMatrix' num2str(k) '.mat'])

%% commenting figures
% figure
% 
% b= imagesc(KullDiv)
% set(b,'AlphaData',(KullDiv~=0)) % set to white 0 values
% 
% 
% cb = colorbar
% maxValue = max(max(KullDiv));
% minValue = min(min(KullDiv(KullDiv ~= 0)));
% caxis([minValue, maxValue]);
% set(gca,'linew',3)
% 
% 
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on

% figure
% 
% b= imagesc(AverageMatrix)
% set(b,'AlphaData',(AverageMatrix~=0)) % set to white 0 values
% 
% 
% cb = colorbar
% maxValue = max(max(AverageMatrix));
% minValue = min(min(AverageMatrix(AverageMatrix ~= 0)));
% caxis([minValue, maxValue]);
% set(gca,'linew',3)
% 
% title('Average')
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on
% 
% 
% %% with error bars
% 
% % figure
% % plot(mean(distances))
% 
% figure
% plot((mean(AverageMatrix)),'LineWidth',3)
% 
% y = mean(AverageMatrix);
% x = 1:length(y);
% 
% sd_vct = std(AverageMatrix);
% hold on
% errorbar(x,y,sd_vct, 'bp')
% 
% set(gca,'linew',3)
% ylabel('Kullback-Leibler')
% ytickformat('%.1f')
% xlabel('Starting position')
% title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
% 
% %set(gca, 'YTick', round(1:4))
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on
% 
% 
% %% just means
% 
% % figure
% plot(mean(AverageMatrix))
% 
% set(gca,'linew',3)
% ylabel('Kullback-Leibler')
% ytickformat('%.1f')
% xlabel('Starting position')
% title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
% 
% %set(gca, 'YTick', round(1:4))
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on