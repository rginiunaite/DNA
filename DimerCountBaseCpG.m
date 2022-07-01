% explore sequences with different number of CpG dinucleotides

clear all
i1=1;
i2=24;
ka=1;
for ich =i1:i2%[1,5,6,7,9,10,11,13,14,15,16,17,18,21,22,23,24]%i1:i2
ich
filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpG147.mat',ich);

datagroup = load(filename);

dimersALL=0;
dimersALLnot=0;
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	
	Basecount = basecount(seq);
	
	baseA(i) = Basecount.A;
	baseT(i) = Basecount.T;
	baseC(i) = Basecount.C;
	baseG(i) = Basecount.G;
	
	dimer(1,i) = Dimers.AA;
	dimer(2,i) = Dimers.AC;
	dimer(3,i) = Dimers.AG;
	dimer(4,i) = Dimers.AT;
	dimer(5,i) = Dimers.CA;
	dimer(6,i) = Dimers.CC;
	dimer(7,i) = Dimers.CG;
	dimer(8,i) = Dimers.CT;
	dimer(9,i) = Dimers.GA;
	dimer(10,i) = Dimers.GC;
	dimer(11,i) = Dimers.GG;
	dimer(12,i) = Dimers.GT;
	dimer(13,i) = Dimers.TA;
	dimer(14,i) = Dimers.TC;
	dimer(15,i) = Dimers.TG;
	dimer(16,i) = Dimers.TT;
	k=k+1;
	end
end


start = length(datagroup.allsubsequences);

%% add notCpG 

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%inotCpG147.mat',ich);

datagroup = load(filename);

dimersALL=0;
dimersALLnot=0;
k=1;
for i=start + 1:start + length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i-start);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	
	Basecount = basecount(seq);
	
	baseA(i) = Basecount.A;
	baseT(i) = Basecount.T;
	baseC(i) = Basecount.C;
	baseG(i) = Basecount.G;
	
	dimer(1,i) = Dimers.AA;
	dimer(2,i) = Dimers.AC;
	dimer(3,i) = Dimers.AG;
	dimer(4,i) = Dimers.AT;
	dimer(5,i) = Dimers.CA;
	dimer(6,i) = Dimers.CC;
	dimer(7,i) = Dimers.CG;
	dimer(8,i) = Dimers.CT;
	dimer(9,i) = Dimers.GA;
	dimer(10,i) = Dimers.GC;
	dimer(11,i) = Dimers.GG;
	dimer(12,i) = Dimers.GT;
	dimer(13,i) = Dimers.TA;
	dimer(14,i) = Dimers.TC;
	dimer(15,i) = Dimers.TG;
	dimer(16,i) = Dimers.TT;
	k=k+1;
	end
end

indices1 = dimer(7,:)<15 & dimer(7,:)>=5; % indices with a specific number of CpG

for k=1:16
	meanvalch(ka,k)= mean(dimer(k,indices1));
	stdvalch(ka,k) = std(dimer(k,indices1));
end

ka=ka+1;
end



figure


bar(meanvalch')
% from one chromosome
%bar([mean(dimer(1,indices1)),mean(dimer(2,indices1)),mean(dimer(3,indices1)),mean(dimer(4,indices1)),mean(dimer(5,indices1)),mean(dimer(6,indices1)),mean(dimer(7,indices1)),mean(dimer(8,indices1)),mean(dimer(9,indices1)),mean(dimer(10,indices1)),mean(dimer(11,indices1)),mean(dimer(12,indices1)),mean(dimer(13,indices1)),mean(dimer(14,indices1)),mean(dimer(15,indices1)),mean(dimer(16,indices1))])

ax=gca;
box on
grid on
%xlabel('Sequence location')
%title('Lowest PL')
ylabel('Av. number')
%xlim([10,60])
%xticks([1 2 3 4])
%xticklabels({'A','T','C','G'})

xticks([1:16])
xticklabels({'AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT'})

set(gca,'FontSize',16)

ft = 'Times';
fsz = 16;
% < 15
title('25 \leq CpG','FontSize',16,'FontName',ft)
  
%ylim([0,17])
set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

figure


bar(stdvalch')
% from one chromosome
%bar([mean(dimer(1,indices1)),mean(dimer(2,indices1)),mean(dimer(3,indices1)),mean(dimer(4,indices1)),mean(dimer(5,indices1)),mean(dimer(6,indices1)),mean(dimer(7,indices1)),mean(dimer(8,indices1)),mean(dimer(9,indices1)),mean(dimer(10,indices1)),mean(dimer(11,indices1)),mean(dimer(12,indices1)),mean(dimer(13,indices1)),mean(dimer(14,indices1)),mean(dimer(15,indices1)),mean(dimer(16,indices1))])

ax=gca;
box on
grid on
%xlabel('Sequence location')
%title('Lowest PL')
ylabel('Av. number')
%xlim([10,60])
%xticks([1 2 3 4])
%xticklabels({'A','T','C','G'})

xticks([1:16])
xticklabels({'AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT'})

set(gca,'FontSize',16)

ft = 'Times';
fsz = 16;

title('5 \leq CpG < 15, standard deviation','FontSize',16,'FontName',ft)
  
%ylim([0,17])
set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)





