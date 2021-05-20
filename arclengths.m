clear all
% Arc lengths of subsequences

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
% idx2 = strcmp('Mouse',{data.Seq.group})
% idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth

% 
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

%% for all

seqnum = length(datagroup);

lenseq =73;

distances = zeros(seqnum,400-lenseq-1);


k=1; %count;
for j = 1:seqnum

    seq = datagroup(j).S; % from nucleosomes data
    %seq = datagroup(j).Seq; % for random sequences

    
    if length(seq)==400
        
        [shapes, stiff] = constructSeqParms(seq, paramset);

        abs_coord = frames(shapes); % relative to absolute coordinates

        for i=1:400-lenseq -1
            sumall = 0;
            for nr = 1:lenseq
                sumall = sumall + norm(abs_coord(i+nr).rc -abs_coord(i+nr+1).rc);
            end
            distances(k,i) = sumall;
        end
        k=k+1;
    end

end


% 
% %% only means
% 
% figure
% %plot(distances(1,:))
% %hold on
% %plot(distances(2,:))
% plot((mean(distances)),'LineWidth',3)
% 
% set(gca,'linew',3)
% ylabel('Arc length')
% ytickformat('%.2f')
% xlabel('Starting position')
% title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
% 
% %set(gca, 'YTick', round(1:4))
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
%box on

%% with error bars

% figure
% plot(mean(distances))

figure
plot((mean(distances)),'LineWidth',3)

y = mean(distances);
x = 1:400-lenseq-1;

sd_vct = std(distances);
hold on
%errorbar(x,y,sd_vct, 'bp')
shadedErrorBar(x,y,sd_vct/sqrt(length(distances(:,1))));

set(gca,'linew',3)
ylabel('Arc length')
ytickformat('%.1f')
xlabel('Starting position')
title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
xlim([0,400-lenseq-1])
%set(gca, 'YTick', round(1:4))

set(gca,'FontSize',36)
ax = gca;
grid on
box on