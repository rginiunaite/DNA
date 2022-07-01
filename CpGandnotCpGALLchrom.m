% plot CpG and not CpG
i1=1;
i2=24;

meanenergiesCpG = zeros(1,24);
meanenergiesnotCpG = zeros(1,24);
for ich=i1:i2



filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iCpG147.txt',ich);

datagroupCpG = load(filename);

filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotCpG147.txt',ich);

datagroupnotCpG = load(filename);

dimersALL = zeros(1,length(datagroupCpG));

%% pick CpG sequences with the count of CpGs greater than 10. 20/147

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpG147.mat',ich);

datagroup = load(filename);
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	dimersALL(k) = Dimers.CG;
	k=k+1;
end
end

ind = (dimersALL>0); %% only the ones with sufficiently large number of CpGs


%% pick notCpG sequences with the count of CpGs less than 10. 

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%inotCpG147.mat',ich);

datagroup = load(filename);
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	dimersALL(k) = Dimers.CG;
	k=k+1;
end
end

%indnot = (dimersALL<10);  %% only the ones with sufficiently small number of CpGs
indnot=(dimersALL<10000);



meanenergiesCpG(ich) = mean(nonzeros(datagroupCpG(ind)));
meanenergiesnotCpG(ich) = mean(nonzeros(datagroupnotCpG(indnot)));


if ich==1

allenergiesCpG=datagroupCpG(ind);
allenergiesnotCpG=datagroupnotCpG(indnot);

else

allenergiesCpG = [allenergiesCpG, datagroupCpG(ind)];

allenergiesnotCpG=[allenergiesnotCpG, datagroupnotCpG(indnot)];

end


end

    figure
    a1 = histogram(nonzeros(allenergiesCpG/147),50,'FaceColor','b') % for random sequences
    
    
    a1.BinWidth = 0.05;
    xlabel('Energies')
    ylabel('Frequency')
%     
    set(gca,'FontSize',36)
 
    
    
        hold on
    a2 = histogram(nonzeros(allenergiesnotCpG/147),45,'FaceColor','r') % for random sequences
    a2.BinWidth = 0.05;
    xlabel('Energy, kT/bp')
    ylabel('Frequency')
  	legend([a1,a2],'CGIs','not CGIs')
%     

grid on
    set(gca,'FontSize',36)
ax=gca;
ax.BoxStyle = 'full';
ax.LineWidth = 3;
         
     ft = 'Times'
fsz = 36;    

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)



figure 

plot(meanenergiesCpG,'r')

hold on

plot(meanenergiesnotCpG,'b')
