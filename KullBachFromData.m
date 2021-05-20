% Kullbach-Leibler divergences
str = 'Yeast';
paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');


numberseq = 73;

for i = 1:numberseq
    filename = sprintf('DataKL/KullDiv%i.mat',i);
    
    sepdata(i) = load(filename);
end

AverageMatrix = zeros(35,35);

for i=1:numberseq
    AverageMatrix = AverageMatrix + sepdata(i).KullDiv;
end

AverageMatrix = AverageMatrix/50;

% set diagonal values to zero, they are already approximately 0
for i=1:35
    AverageMatrix(i,i)=0;
end
% 
% % specific curve
% for j=1:10
%     figure
%     b= imagesc(sepdata(j).KullDiv)
%     set(b,'AlphaData',(sepdata(j).KullDiv~=0)) % set to white 0 values
% 
%     cb = colorbar
%     maxValue = max(max(sepdata(j).KullDiv));
%     minValue = min(min(sepdata(j).KullDiv(sepdata(j).KullDiv ~= 0)));
%     caxis([minValue, maxValue]);
%     set(gca,'linew',3)
%     figure
%     plot(mean(sepdata(j).KullDiv))
% 
%     set(gca,'linew',3)
%     ylabel('Kullback-Leibler')
%     ytickformat('%.1f')
%     xlabel('Starting position')
%     title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])
%     
% end

filename = sprintf('DataKL/AverageMatrix74.mat');
    
dataAv= load(filename);

AverageMatrix = dataAv.AverageMatrix;

AverageMatrix = AverageMatrix/numberseq;

for i=1:35
    AverageMatrix(i,i)=0;
end
figure
% 
b= imagesc(AverageMatrix)
set(b,'AlphaData',(AverageMatrix~=0)) % set to white 0 values

cb = colorbar
maxValue = max(max(AverageMatrix));
minValue = min(min(AverageMatrix(AverageMatrix ~= 0)));
caxis([minValue, maxValue]);
set(gca,'linew',3)

title('Average')

set(gca,'FontSize',36)
ax = gca;
grid on
box on


%% with error bars

% figure
% plot(mean(distances))

figure
plot((mean(AverageMatrix)),'LineWidth',3)

y = mean(AverageMatrix);
x = 1:length(y);

sd_vct = std(AverageMatrix);
hold on
errorbar(x,y,sd_vct, 'bp')

set(gca,'linew',3)
ylabel('Kullback-Leibler')
ytickformat('%.1f')
xlabel('Starting position')
%title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])

%set(gca, 'YTick', round(1:4))

set(gca,'FontSize',36)
ax = gca;
grid on
box on


%% just means

figure
plot(mean(AverageMatrix))

set(gca,'linew',3)
ylabel('Kullback-Leibler')
ytickformat('%.1f')
xlabel('Starting position')
%title([str,' (', num2str(k-1),  '), sequence length ' , num2str(lenseq) ])

%set(gca, 'YTick', round(1:4))

set(gca,'FontSize',36)
ax = gca;
grid on
box on
