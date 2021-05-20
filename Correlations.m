%% correlation coefficient yeast and random sequences


paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');
%seqnum = 20;
fulllen = 400;
%from already generated files
% distances = load('distances100sequences.mat');
% 
% distances = cell2mat(struct2cell(distances)); %;zeros(seqnum,350);

%% select specific groups/organisms
vec = [data.Seq.group];

str = 'Yeast';
idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
% idx2 = strcmp('Mouse',{data.Seq.group})
%idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth


% %% all sequences
% str = 'All';
% datagroup = data.Seq; %% all sequences

% %% random sequences
% 
%  str = 'Random';
% % for j =1 :50
% %     datagroup(j).Seq = [random_dna(400)];
% % end
% % 
% %load 500 previously generated random sequences
% data500 = load('Random500Sequences400length.mat');
% 
% datagroup = data500.datagroup;


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

Int = zeros(1,fulllen-lenseq-lenseq);
minimumarea = zeros (1,seqnum);
for j = 1:seqnum
    for i = 1:fulllen-lenseq-lenseq
        rA = [i:lenseq+i];
        Int(i) = trapz(rA, distances(j,rA));
        
    end
    [value, index] = min(Int);
    minimumarea(j) = index;
end


figure 
minimumarea = minimumarea +36;
histogram(minimumarea,10)

title('Minimum area')


coeff = zeros(seqnum,seqnum);

% find average position of the 73 minimum positions
avmin = zeros(1,seqnum);
for j=1:seqnum
    x =smooth((distances(j,:)));

    [xs,index] = sort(x);
    
%     [B,TF] = rmoutliers(index(1:146),'median');
%     avmin(j) = mean(B);
%     sd_vct(j) = std(B);
    
    avmin(j) = mean((index(1:lenseq)));
    sd_vct(j) = std((index(1:lenseq)));
    
    
end

% figure 
% 
% scatter([1:seqnum],avmin)  
% %     y = avmin;
% %     x = 1:seqnum;
% %    hold on
% % errorbar(x,y,sd_vct, 'bp')
% xlabel('Sequence number')
% ylabel('Average minimal position')
% ylim([0,400-lenseq])
% title(str)
% 
% set(gca,'linew',3)
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on


figure
histogram(avmin)
title(str)

ylabel('Frequency')
%ytickformat('%.1f')
xlabel('Average minimum values')
%title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])

%set(gca, 'YTick', round(1:4))

set(gca,'FontSize',36)
ax = gca;
box on


% pd = fitdist(avmin','Normal')
% 
% %% frequency
% figure
% avminnorm = avmin - mean(avmin);
% 
% [pxx,f] = periodogram(avmin,[],[],10);
% 
% plot(f,pxx)
% ax = gca;
% ax.XLim = [0 10];
% xlabel('Frequency (cycles/week)')
% ylabel('Magnitude')



% for j = 1:seqnum
%     for k = 1:seqnum
%         
%         
%         %A = corrcoef(distances(j,:),distances(k,:));
%         coeff(j,k) = A(1,2); % take one entry from a symmetric matrix
%     end
% end


% commenting figures
% figure
% 
% b= imagesc(coeff)
% cb = colorbar
%set(b,'AlphaData',(coeff~=0)) % set to white 0 values


% 
% maxValue = max(max(coeff));
% minValue = min(min(coeff(coeff ~= 0)));
% caxis([minValue, maxValue]);
% set(gca,'linew',3)
% 
% 
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on
% 
