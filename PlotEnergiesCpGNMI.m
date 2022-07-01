%% Compare ALLenergies CpG NMI
addpath  /home/rasa/Documents/DNA/cgDNA+_clean/cgDNA+_clean/optimisation/data
addpath  /home/rasa/Documents/DNA/cgDNA+_clean/cgDNA+_clean/SkirmantasData/EnergiesCh37
%filename = sprintf('ALLenergiesCpGintersectNMI146.txt');
filename = sprintf('energiesCh1CpGintersectNMI147.txt');
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data1 = load(filename);

averagesCpGNMI= mean(data1/147);
stdsCpGNMI = std(data1/147);


%filename = sprintf('ALLenergiesNMINOTintersectCpG146.txt');
filename = sprintf('energiesCh1NMINOTintersectCpG147.txt');
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data2 = load(filename);

%data2 = data2(1:2000);

averagesNMINOTCpG= mean(data2/147)
stdsNMINOTCpG = std(data2/147)


%filename = sprintf('ALLenergiesCpGNotintersectNMI146.txt');
filename = sprintf('energiesCh1CpGNOTintersectNMI147.txt');
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data3 = load(filename);

averagesCpGnotNMI= mean(data3/147)
stdsCpGnotNMI = std(data3/147)

%filename = sprintf('ALLenergiesneither146.txt');
filename = sprintf('energiesCh1neither147.txt');
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data4 = load(filename);

%data4 = data4(1:2000);

averagesneither= mean(data4/147)
stdsneither = std(data4/147)

%filename = sprintf('ALLenergiesMethCpGNOTintersectNMI146.txt');
filename = sprintf('energiesCh1MpNNOTintersectNMI147.txt');
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data5 = load(filename);

averagesmethCpGnotNMI= mean(data5/147)
stdsmethCpGnotNMI = std(data5/147)


%filename = sprintf('ALLenergiesMethneither146.txt');
filename = sprintf('energiesCh1neithermeth147.txt');
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data6 = load(filename);

data6 = nonzeros(data6);

averagesmethneither= mean(data6/147)
stdsmethneither = std(data6/147)

    figure
    h1 = histogram(data1(1:100),10,'FaceColor','b') % for random sequences
    xlabel('Energies')
    ylabel('Frequency')
%     
    set(gca,'FontSize',36)
    ax = gca;
    box on
    
    hold on
        h2 = histogram(data2(1:100),10,'FaceColor','r') 
      
        h3 = histogram(data3(1:100),10,'FaceColor','y') 
        
        h4 = histogram(data4(1:100),10,'FaceColor','m') 
    
    
    legend('CpG and NMI','NMI not CpG','CpG not NMI','Not CpG and Not NMI')
    
    


figure


cpg = scatter(15, averagesCpGNMI,'k','LineWidth',6);

hold on 
h=errorbar(15,averagesCpGNMI, stdsCpGNMI,'k.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;
hold on

mpn = scatter(20, averagesNMINOTCpG,'r','LineWidth',6);
hold on 
h=errorbar(20,averagesNMINOTCpG, stdsNMINOTCpG,'r.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;

%hpk = scatter(25, averagesCpGnotNMI,'b','LineWidth',6);
%hold on 
%h=errorbar(25,averagesCpGnotNMI, stdsCpGnotNMI,'b.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;


%notcpg = scatter(30, averagesneither,'g','LineWidth',6);
%hold on 
%h=errorbar(30,averagesneither, stdsneither,'g.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;

ytickformat('%.1f')

methcpgnotNMI = scatter(25, averagesmethCpGnotNMI,'b','LineWidth',6);
hold on 
h=errorbar(25,averagesmethCpGnotNMI, stdsmethCpGnotNMI,'b.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;

notcpg = scatter(30, averagesmethneither,'c','LineWidth',6);
hold on 
h=errorbar(30,averagesmethneither, stdsmethneither,'g.','linewidth',6)
set(gca,'FontSize',36)
h.CapSize = 12;


ax=gca;
box on
grid on
xlabel('Sequence location')
ylabel('Energy, kT/bp')
xlim([10,35])
xticks([15 20 25 30])
xticklabels({'CpG and NMI','NMI not CpG','CpG not NMI','Not CpG and Not NMI'})%, 'Meth CpG not NMI', 'MethNeither'})

%legend([cpg mpn hpk],{'CpG','MpN', 'HpK'})
    
