% Persistence length versus energy   


filename = sprintf('data/energiesLinkerNEWALL.txt');
dataenergiesLinker = load(filename);

filename = sprintf('data/energiesNucleosomesNEWALL.txt');
dataenergiesNucleosomes = load(filename);

filename = sprintf('data/mineigenvaluesLinkerALL.txt');
datamineigLinker = load(filename);

filename = sprintf('data/mineigenvaluesNucleosomesALL.txt');
datamineigNucleosomes = load(filename);

filename = sprintf('data/maxeigenvaluesLinkerALL.txt');
datamaxeigLinker = load(filename);

filename = sprintf('data/maxeigenvaluesNucleosomesALL.txt');
datamaxeigNucleosomes = load(filename);

filename = sprintf('data/meaneigenvaluesLinkerALL.txt');
datameaneigLinker = load(filename);

filename = sprintf('data/meaneigenvaluesNucleosomesALL.txt');
datameaneigNucleosomes = load(filename);

filename = sprintf('../PLdataNewParam/LinkerPl/Plt0Dyn1.txt');
    dataLinkerPL = load(filename);
    
filename = sprintf('../PLdataNewParam/NucleosomesPl/Plt0Dyn1.txt');
    dataNucleosomesPL = load(filename);
    
    
filename = sprintf('data/energiesLinkerALL.txt');
dataenergiesLinker = load(filename);

filename = sprintf('data/energiesNucleosomesALL.txt');
dataenergiesNucleosomes = load(filename);

   
filename = sprintf('data/arcLengthsLinkerALL.txt');
dataarcLinker = load(filename);

filename = sprintf('data/arcLengthsNucleosomesALL.txt');
dataarcNucleosomes = load(filename);


filename = sprintf('data/distancestartendLinkerALL.txt');
datadistLinker = load(filename);

filename = sprintf('data/distancestartendNucleosomesALL.txt');
datadistNucleosomes = load(filename);


mean(dataenergiesLinker)
mean(dataenergiesNucleosomes)


xL=datamineigLinker;
yL=datamaxeigLinker;

xN = datamineigNucleosomes;
yN = datamaxeigNucleosomes;

meanL = datameaneigLinker;
meanN = datameaneigNucleosomes;
corrcoef(xL, yL)
corrcoef(xN,yN)

mean(mean(xL))
mean(mean(xN))

mean(mean(yL))
mean(mean(yN))

mean(meanL)
mean(meanN)



%std(xL)
%std(xN)

corrcoef(dataLinkerPL,dataenergiesLinker )
corrcoef(dataNucleosomesPL,dataenergiesNucleosomes )






% energies vs PL nucleosomes
 figure
      y = dataenergiesNucleosomes;
   x=dataNucleosomesPL;
   scatter(x,y,'filled')
    P = polyfit(x,y,1);
    yfit = P(1)*x+P(2);
    hold on;
    plot(x,yfit,'r-.');
    ylabel('Energies')
    xlabel('Persistence length, bp')
    
ft = 'Times';
fsz = 36;  

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

%     
    set(gca,'FontSize',36)
    ax = gca;
    box on
    
    style = hgexport('factorystyle');
style.Color = 'gray'

style = hgexport('factorystyle');
style.Color = 'gray'
hgexport(gcf,'test.eps',style)
%     

% energies vs PL linker
figure
  y = dataenergiesLinker;
  x=dataLinkerPL;
  scatter(x,y,'filled')
   P = polyfit(x,y,1);
   yfit = P(1)*x+P(2);
   hold on;
   plot(x,yfit,'r-.','LineWidth',2);
    grid on
   ylabel('Energies')
   xlabel('Persistence length, bp')
    
   set(gca,'FontSize',36)
   ax = gca;
   box on
   
    
    style = hgexport('factorystyle');
style.Color = 'gray'

style = hgexport('factorystyle');
style.Color = 'gray'
hgexport(gcf,'test.eps',style)
    
    
    % nucleosomes PL
 figure
    %histogram(dataenergiesNucleosomes,50,'FaceColor','b') % nucleosomes
    nucl = histogram(dataNucleosomesPL,'FaceColor','k') % nucleosomes
    nucl.BinWidth = 0.8;
    

    xlabel('Energy')
    ylabel('Frequency')
%     
    set(gca,'FontSize',36)
    ax = gca;
    box on
    ft = 'Times';
fsz = 36;  

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

    
hold on 
% linker PL 
    %histogram(dataenergiesLinker,50,'FaceColor','r') % Linker
    linker = histogram(dataLinkerPL,'FaceColor',[0.5 0.5 0.5]) % Linker 
    linker.BinWidth = 0.8;  
    %xlabel('Energy')
    xlabel('Persistence length, bp')
    ylabel('Frequency')

	ylim([0,155])
    set(gca,'FontSize',36)
    ax = gca;
    box on
    grid on
    legend('Nucleosomal','Linker','FontName', ft)
    
    
style = hgexport('factorystyle');
style.Color = 'gray'

style = hgexport('factorystyle');
style.Color = 'gray'
hgexport(gcf,'test.eps',style)


ft = 'Times';
fsz = 36;  

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

ax=gca;
grid on
%box on
ax.BoxStyle = 'full';
ax.LineWidth = 2;


set(gca,'TickLabelInterpreter','latex')






