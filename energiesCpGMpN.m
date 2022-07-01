%% plot CpG data

nr =6;

meanCpG = zeros(1,nr);
stdCpG = zeros(1,nr);
meanMpN = zeros(1,nr);
stdMpN = zeros(1,nr);
meanHpK = zeros(1,nr);
stdHpK = zeros(1,nr);
    for j=1:nr
        i = (j-1)*10;
        %filename = sprintf('data/energiesSeqCpG%i.txt',i);
    	filename = sprintf('data/DiffNrCpGNEW/energiesSeq1000CpG%iLength147.txt',i);
    	data = load(filename);
    	CpGenergies = data/147; 
    
    	%filename = sprintf('data/energiesSeqMpN%i.txt',i);
    	filename = sprintf('data/DiffNrCpGNEW/energiesSeq1000MpN%iLength147.txt',i);
    	data = load(filename);
    	MpNenergies = data/147; 
    
        %filename = sprintf('data/energiesSeqHpK%i.txt',i);
    	filename = sprintf('data/DiffNrCpGNEW/energiesSeq1000HpK%iLength147.txt',i);
    	data = load(filename);
    	HpKenergies = data/147; 
    
    	meanCpG(j) = mean(CpGenergies);
    	stdCpG(j) = std(CpGenergies);
    	meanMpN(j) = mean(MpNenergies);
    	stdMpN(j) = std(MpNenergies);
    	meanHpK(j) = mean(HpKenergies);
    	stdHpK(j) = std(HpKenergies);
    
    end
    
    
    figure 
    
    x= [0,10,20,30,40,50];
    cpg = scatter(x,meanCpG,'k','LineWidth',6)

    hold on
    h = errorbar (x,meanCpG, stdCpG,'k.','linewidth',6)
 	h.CapSize = 12;
    hold on
    mpn = scatter(x,meanMpN,'r','LineWidth',6)
    h = errorbar (x,meanMpN, stdMpN,'r.','linewidth',6)
    	h.CapSize = 12;
    hpk = scatter(x,meanHpK,'b','LineWidth',6)
    h = errorbar (x, meanHpK,stdHpK,'b.','linewidth',6)

	h.CapSize = 12;
    xlabel('Number of CpGs')
    ylabel('Energy, kT/bp')
    ytickformat('%,.1f')
    set(gca,'FontSize',36)
    ax = gca;
    box on
    grid on
    legend([cpg,mpn,hpk],'CpG','MpN','HpK')
    
         ft = 'Times'
fsz = 36;    

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

    


