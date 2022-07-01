    %% Plot CpG and not CpG, MpN and not MpN, HpK and not HpK
    
    % cpg islands
    nbins=50;
    
    i1=2;
    i2=2;
    
    for ich=[5,6,7,9,10,11,13,14,15,16,17,18,21,22,23,24]
    filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iCpG147.txt',ich);
    data = load(filename);
    data(data == 0) = [];
    CpGenergies = nonzeros(data/147); 
    
    filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iMpN147.txt',ich);
    data = load(filename);
    MpNenergies = nonzeros(data/147); 
    
    filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iHpK147.txt',ich);
    data = load(filename);
    HpKenergies = nonzeros(data/147); 
    
    %filename = sprintf('data/energiesSeqnotcpgislands.txt');
    filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotCpG147.txt',ich);
    data = load(filename);
    data(data == 0) = [];
    notCpGenergies = nonzeros(data/147); 
    
    filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotMpN147.txt',ich);
    data = load(filename);
    notMpNenergies = nonzeros(data/147); 
    
    filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotHpK147.txt',ich);
    data = load(filename);
    notHpKenergies = nonzeros(data/147); 
    
    
    length(CpGenergies)
    
    length(MpNenergies)
    
    cpgmean = mean(CpGenergies)
    cpgstd = std(CpGenergies)
    
    
    mpnmean = mean(MpNenergies)
    mpnstd = std(MpNenergies)
    
    
    hpkmean = mean(HpKenergies)
    hpkstd = std(HpKenergies)
    
    notcpgmean = mean(notCpGenergies)
    notcpgstd = std(notCpGenergies)
    
    notmpnmean = mean(notMpNenergies)
    notmpnstd = std(notMpNenergies)
    
    
    nothpkmean = mean(notHpKenergies)
    nothpkstd = std(notHpKenergies)
    
    
    %diff = MpNenergies-CpGenergies;
    
        
        %figure
        
        %histogram(diff,nbins)
        %xlabel('Difference in energy (MpN-CpG)')
        %ylabel('Frequency')
        
figure 
cpgener = histogram(CpGenergies)

cpgener.BinWidth = 0.1;
%title('20 CpGs')
%ylim([0,400])
%title('Human Ch1 CpG islands')


hold on
mpnener = histogram(notCpGenergies)
mpnener.BinWidth = 0.1;
%title('20 CpGs')
%ylim([0,400])
%title('Human Ch1 CpG islands')
xlabel('Energy, kT')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on

legend([cpgener,mpnener],'CpG','not CpG')
ft = 'Times';
fsz = 36;  
set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
xtickformat('%,.1f')

% plot all

figure


cpg = scatter(15, cpgmean,'k','LineWidth',6);

hold on 
h=errorbar(15,cpgmean, cpgstd,'k.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;
hold on

mpn = scatter(20, mpnmean,'r','LineWidth',6);

hold on 
h=errorbar(20,mpnmean, mpnstd,'r.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;

hpk = scatter(25, hpkmean,'b','LineWidth',6);

hold on 
h=errorbar(25,hpkmean, hpkstd,'b.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;
notcpg = scatter(45, notcpgmean,'k','LineWidth',6);

hold on 
h=errorbar(45,notcpgmean, notcpgstd,'k.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;
notmpn = scatter(50, notmpnmean,'r','LineWidth',6);

hold on 
h=errorbar(50,notmpnmean, notmpnstd,'r.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;

nothpk = scatter(55, nothpkmean,'b','LineWidth',6);

hold on 
h=errorbar(55,nothpkmean, nothpkstd,'b.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;

ax=gca;
box on
grid on
xlabel('Sequence location')
ylabel('Energy, kT/bp')
xlim([10,60])
xticks([20 50])
xticklabels({'CGIs','Not CGIs'})

legend([cpg mpn hpk],{'CpG','MpN', 'HpK'})
ft = 'Times';
fsz = 36;  
set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)
ytickformat('%,.1f')
ylim([4.5,6.5])
end
    
