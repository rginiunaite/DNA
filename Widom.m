%% correlation coefficient yeast and random sequences


paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
fileName = 'Widom.txt';
% my sample text contains ---> apple, baseball, car, donut, & elephant in single column.
FID = fopen(fileName);
data = textscan(FID,'%s');
fclose(FID);
stringData = string(data{:});

%stringData = 'ATGACTAAATGCTTGCATCACAATACTTGAAGTTGACAATATTATTTAAGGACCTATTGTTTTTTCCAATAGGTGGTTAGCAATCGTCTTACTTTCTAACTTTTCTTACCTTTTACATTTCAGCAATATATATATATATTTCAAGGATATACCATTCTAATGTCTGCCCCTAAGAAGATCGTCGTTTTGCCAGGTGACCACGTTGGTCAAGAAATCACAGCCGAAGCCATTAAGGTTCTTAAAGCTATTTCTGATGTTCGTTCCAATGTCAAGTTCGATTTCGAAAATCATTTAATTGGTGGTGCTGCTATCGATGCTACAGGTGTCCCACTTCCAGATGAGGCGCTGGAAGCCTCCAAGAAGGTTGATGCCGTTTTGTTAGGTGCTGTGGGTGGTCCTA';


str = 'Widom';






%% for all
seqnum = 1;


seq = convertStringsToChars(stringData);

fulllen = length(seq);

lenseq =73;


%% norm
% 
% distances = zeros(seqnum,fulllen-lenseq);
% 
% 
% k=1; %count;
% for j = 1:seqnum
% 
%  
%         [shapes, stiff] = constructSeqParms(seq, paramset);
% 
%         abs_coord = frames(shapes); % relative to absolute coordinates
% 
%         for i=1:fulllen-lenseq 
%             distances(k,i) = norm(abs_coord(i).rc -abs_coord(i+lenseq).rc);
% 
%         end
%         k=k+1;
%    
% 
% end

%% arc length

distances = zeros(seqnum,fulllen-lenseq-1);


k=1; %count;
for j = 1:seqnum

 
        [shapes, stiff] = constructSeqParms(seq, paramset);

        abs_coord = frames(shapes); % relative to absolute coordinates

        for i=1:fulllen-lenseq -1
            sumall = 0;
            for nr = 1:lenseq
                sumall = sumall + norm(abs_coord(i+nr).rc -abs_coord(i+nr+1).rc);
            end
            distances(k,i) = sumall;
        end
   

end


%% compare areas under the curve

% 
% 
%     for i = 1:fulllen-lenseq-lenseq
%         rA = [i:lenseq+i];
%         Int(i) = trapz(rA, distances(j,rA));
%         
%     end
% 
% [value, index] = min(Int);
% 
 figure
plot(distances(1,:));




set(gca,'linew',3)
ylabel('Distance')
ytickformat('%.1f')
xlabel('Starting position')
title([str])%,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
xlim([0,fulllen-lenseq])

%set(gca, 'YTick', round(1:4))

set(gca,'FontSize',36)
ax = gca;
grid on
box on


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

figure 

scatter([1:seqnum],avmin)  
%     y = avmin;
%     x = 1:seqnum;
%    hold on
% errorbar(x,y,sd_vct, 'bp')
xlabel('Sequence number')
ylabel('Average minimal position')
ylim([0,fulllen-lenseq])
title(str)
% 
% set(gca,'linew',3)
% 
% set(gca,'FontSize',36)
% ax = gca;
% grid on
% box on


figure
histogram(avmin,20)
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
