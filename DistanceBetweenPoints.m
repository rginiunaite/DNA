clear all
% Distance from midpoints

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
% idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth


% %% all sequences
% str = 'All';
% datagroup = data.Seq; %% all sequences

% %% random sequences
% 
%  str = 'Random';
% % % for j =1 :50
% % %     datagroup(j).Seq = [random_dna(400)];
% % % end
% % 
% %load 500 previously generated random sequences
% data500 = load('Random500Sequences400length.mat');

%datagroup = data500.datagroup;

%% for all
seqnum = length(datagroup);

lenseq =73;

distances = zeros(seqnum,400-lenseq);


k=1; %count;
for j = 1:seqnum

    seq = datagroup(j).S; % from nucleosomes data
    %seq = datagroup(j).Seq; % for random sequences
    
    
    if length(seq)==400
        
        [shapes, stiff] = constructSeqParms(seq, paramset);

        abs_coord = frames(shapes); % relative to absolute coordinates

        for i=1:400-lenseq 
            distances(k,i) = norm(abs_coord(i).rc -abs_coord(i+lenseq).rc);

        end
        k=k+1;
    end

end




% figure
% plot(mean(distances))

%% Only means Figure


% figure
% plot((mean(distances)),'LineWidth',3)
% 
% 
% set(gca,'linew',3)
% ylabel('Distance')
% ytickformat('%.2f')
% xlabel('Starting position')
% title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
% 
% %set(gca, 'YTick', round(1:4))
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on



%% Figure with error bars
figure
plot((mean(distances)),'LineWidth',3)

y = mean(distances);
x = 1:400-lenseq;

sd_vct = std(distances);
hold on
shadedErrorBar(x,y,sd_vct/sqrt(length(distances(:,1))));

%H(1) = shadedErrorBar(x, y, {@mean, @(x) 2*std(x)  }, '-r', 0);


set(gca,'linew',3)
ylabel('Distance')
ytickformat('%.1f')
xlabel('Starting position')
title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
xlim([0,400-lenseq])
%set(gca, 'YTick', round(1:4))

set(gca,'FontSize',36)
ax = gca;
grid on
box on