%%% explore sequences with different number of CpG dinucleotides for random sequences

clear all

i1=1;
i2=24;



for ich=i1:i2


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


figure

bar([mean(dimer(1,:)),mean(dimer(2,:)),mean(dimer(3,:)),mean(dimer(4,:)),mean(dimer(5,:)),mean(dimer(6,:)),mean(dimer(7,:)),mean(dimer(8,:)),mean(dimer(9,:)),mean(dimer(10,:)),mean(dimer(11,:)),mean(dimer(12,:)),mean(dimer(13,:)),mean(dimer(14,:)),mean(dimer(15,:)),mean(dimer(16,:))])

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

title('10 CpGs','FontSize',16,'FontName',ft)
ylim([0,17])  

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

end


