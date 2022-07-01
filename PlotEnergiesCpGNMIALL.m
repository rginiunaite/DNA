%% Compare ALLenergies CpG NMI
addpath  /home/rasa/Documents/DNA/cgDNA+_clean/cgDNA+_clean/optimisation/data
addpath  /home/rasa/Documents/DNA/cgDNA+_clean/cgDNA+_clean/SkirmantasData/EnergiesCh37


i1=1;
i2=1;
k=1;
for ich = [1,5,6,9,10,13,14,15,17,18,21,22,23,24]%i1:i2

ich
%filename = sprintf('ALLenergiesCpGintersectNMI146.txt');
filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iCpGintersectNMI147.txt',ich);
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data1 = load(filename);
if ich==1
a1=data1;
else
a1 = [a1,data1];
end
averagesCpGNMI= mean(a1/147);
stdsCpGNMI = std(a1/147);

lena1(k) = length(data1);


%filename = sprintf('ALLenergiesNMINOTintersectCpG146.txt');
filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iNMINOTintersectCpG147.txt',ich);
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data2 = load(filename);
if ich==1
a2=data2;
else
a2 = [a2,data2];
end

lena2(k) = length(data2);
%data2 = data2(1:2000);

averagesNMINOTCpG= mean(a2/147)
stdsNMINOTCpG = std(a2/147)


%filename = sprintf('ALLenergiesCpGNotintersectNMI146.txt');
%filename = sprintf('energiesCh1CpGNOTintersectNMI147.txt');
filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iCpGNOTintersectNMI147.txt',ich);
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data3 = load(filename);
if ich==1
a3=data3;
else
a3 = [a3,data3];
end

lena3(k) = length(data3);

averagesCpGnotNMI= mean(a3/147)
stdsCpGnotNMI = std(a3/147)

%filename = sprintf('ALLenergiesneither146.txt');
%filename = sprintf('energiesCh1neither147.txt');
filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%ineither147.txt',ich);
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data4 = load(filename);
if ich==1
a4=data4;
else
a4 = [a4,data4];
end


lena4(k) = length(data4);
%data4 = data4(1:2000);

averagesneither= mean(a4/147)
stdsneither = std(a4/147)

%filename = sprintf('ALLenergiesMethCpGNOTintersectNMI146.txt');
filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iMpNNOTintersectNMI147.txt',ich);
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data5 = load(filename);
if ich==1
a5 = data5;
else
a5 = [a5,data5];
end

lena5(k) = length(data5);
averagesmethCpGnotNMI= mean(a5/147)
stdsmethCpGnotNMI = std(a5/147)


%filename = sprintf('ALLenergiesMethneither146.txt');
filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%ineithermeth147.txt',ich);
%filename = sprintf('data/ALLenergiesWidomEvery10.txt'); % HumanSeq1
%filename = sprintf('data/ALLenergiesHumanEvery10.txt');
data6 = load(filename);
if ich==1
a6=data6;
else
a6 = [a6,data6];
end

a6 = nonzeros(a6)';

lena6(k) = length(nonzeros(data6));

averagesmethneither= mean(a6/147)
stdsmethneither = std(a6/147)
k=k+1;
end
    %figure
    %h1 = histogram(data1(1:100),10,'FaceColor','b') % for random sequences
    %xlabel('Energies')
    %ylabel('Frequency')
%     
    %set(gca,'FontSize',36)
    %ax = gca;
    %box on
    
    %hold on
     %   h2 = histogram(data2(1:100),10,'FaceColor','r') 
      
     %   h3 = histogram(data3(1:100),10,'FaceColor','y') 
        
     %   h4 = histogram(data4(1:100),10,'FaceColor','m') 
    
    
    %legend('CpG and NMI','NMI not CpG','CpG not NMI','Not CpG and Not NMI')
    
    


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
xticklabels({'CGI and NMI','NMI not CGI','CGI not NMI','Not CGI and not NMI'})%, 'Meth CpG not NMI', 'MethNeither'})

%legend([cpg mpn hpk],{'CpG','MpN', 'HpK'})
    
