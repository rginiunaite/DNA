% I will only look at chromosome 1 and will check whether the number of CpG dinucleotides correlate with the energy

clear all
i1=1;
i2=24;%1,3,4,8,20
ka=1;
one_seq=false;
for ich =[2,5,6,7,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24]%i1:i2
dimersALL = 0;
dimersALLnot = 0;

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpG147.mat',ich);

datagroup = load(filename);

dimersALL=0;
dimersALLnot=0;
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	dimersALL(k) = Dimers.CG;
	k=k+1;
	end
end


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iCpG147.txt',ich);

energiesCpG = load(filename);


%cor = corrcoef(dimersALL(1:end-1),nonzeros(energiesCpG)); % the last one is just zeros, so I remove it.


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%inotCpG147.mat',ich);

datagroup = load(filename);
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	dimersALLnot(k) = Dimers.CG;
	k=k+1;
	end
end


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotCpG147.txt',ich);

energiesnotCpG = load(filename);

%cor = corrcoef(dimersALLnot(1:end-1),nonzeros(energiesnotCpG)); % the last one is just zeros, so I remove it.


x=dimersALL(1:end-1);
x=[x,dimersALLnot(1:end-1)];
y1=nonzeros(energiesCpG/147);
y1=[y1',nonzeros(energiesnotCpG/147)'];

k1=1;
k2=1;
k3=1;
k4=1;
k5=1;
k6=1;
k7=1;
values = zeros(length(x),4);

for i=1:length(x)

if x(i)<5
	values(k1,1)=y1(i);
	k1=k1+1;
elseif x(i)>=5 && x(i)<15
	values(k2,2)=y1(i);
	k2=k2+1;
elseif x(i)>=15 && x(i) < 25
	values(k3,3)=y1(i);
	k3=k3+1;
elseif x(i)>=25 
	values(k4,4)=y1(i);
	k4=k4+1;
end
end


for i=1:4
	meanval1(i,ka) = mean(nonzeros(values(:,i)));
	stdval1(i,ka) = std(nonzeros(values(:,i)));
	
	
end

% check how many sequences with different number of CpGs
chromnum(1,ka)= k1;
chromnum(2,ka)= k2;
chromnum(3,ka)= k3;
chromnum(4,ka)= k4;





%% MpN

dimersALL = 0;
dimersALLnot = 0;
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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iMpN147.txt',ich);

energiesCpG = load(filename);

%cor = corrcoef(dimersALL(1:end-1),nonzeros(energiesCpG)); % the last one is just zeros, so I remove it.


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%inotCpG147.mat',ich);

datagroup = load(filename);
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	dimersALLnot(k) = Dimers.CG;
	k=k+1;
	end
end


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotMpN147.txt',ich);

energiesnotCpG = load(filename);

%cor = corrcoef(dimersALLnot(1:end-1),nonzeros(energiesnotCpG)); % the last one is just zeros, so I remove it.

x=dimersALL(1:end-1);
x=[x,dimersALLnot(1:end-1)];
y1=nonzeros(energiesCpG/147);
y1=[y1',nonzeros(energiesnotCpG/147)'];
k1=1;
k2=1;
k3=1;
k4=1;
k5=1;
k6=1;
k7=1;

values = zeros(length(x),4);
for i=1:length(x)

if x(i)<5
	values(k1,1)=y1(i);
	k1=k1+1;
elseif x(i)>=5 && x(i)<15
	values(k2,2)=y1(i);
	k2=k2+1;
elseif x(i)>=15 && x(i) < 25
	values(k3,3)=y1(i);
	k3=k3+1;
elseif x(i)>=25 
	values(k4,4)=y1(i);
	k4=k4+1;
	
end
end

for i=1:4
	meanval2(i,ka) = mean(nonzeros(values(:,i)));
	stdval2(i,ka) = std(nonzeros(values(:,i)));
end



%% HpK

dimersALL = 0;
dimersALLnot = 0;
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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%iHpK147.txt',ich);

energiesCpG = load(filename);

%cor = corrcoef(dimersALL(1:end-1),nonzeros(energiesCpG)); % the last one is just zeros, so I remove it.


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%inotCpG147.mat',ich);

datagroup = load(filename);
k=1;
for i=1:length(datagroup.allsubsequences)

	seq=datagroup.allsubsequences(i);
	 if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
	Dimers = dimercount(seq);
	dimersALLnot(k) = Dimers.CG;
	k=k+1;
	end
end





filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%inotHpK147.txt',ich);

energiesnotCpG = load(filename);

x=dimersALL(1:end-1);
x=[x,dimersALLnot(1:end-1)];
y1=nonzeros(energiesCpG/147);
y1=[y1',nonzeros(energiesnotCpG/147)'];

k1=1;
k2=1;
k3=1;
k4=1;
k5=1;
k6=1;
k7=1;

values = zeros(length(x),4);

for i=1:length(x)

if x(i)<5
	values(k1,1)=y1(i);
	k1=k1+1;
elseif x(i)>=5 && x(i)<15
	values(k2,2)=y1(i);
	k2=k2+1;
elseif x(i)>=15 && x(i) < 25
	values(k3,3)=y1(i);
	k3=k3+1;
elseif x(i)>=25 
	values(k4,4)=y1(i);
	k4=k4+1;
	
end
end

for i=1:4
	meanval3(i,ka) = mean(nonzeros(values(:,i)));
	stdval3(i,ka) = std(nonzeros(values(:,i)));
end


ka=ka+1; 
end

if one_seq == true
    figure 
    
    %xval= [0,5,10,15,20,25,30];
    xval= [0,10,20,30];
    cpg = scatter(xval,(meanval1'),'k','LineWidth',6)

    hold on
    h = errorbar (xval,(meanval1'), (stdval1'),'k.','linewidth',6)
 	h.CapSize = 12;

    hold on 
    

    mpn = scatter(xval,(meanval2'),'r','LineWidth',6)

    hold on
    h = errorbar (xval,(meanval2'), (stdval2'),'r.','linewidth',6)
 	h.CapSize = 12;

    hold on 
    

    hpk = scatter(xval,(meanval3'),'b','LineWidth',6)

    hold on
    h = errorbar (xval,(meanval3'), (stdval3'),'b.','linewidth',6)
 	h.CapSize = 12;
else
    figure 
    
    %xval= [0,5,10,15,20,25,30];
    xval= [0,10,20,30];
    cpg = scatter(xval,mean(meanval1'),'k','LineWidth',6)

    hold on
    h = errorbar (xval,mean(meanval1'), mean(stdval1'),'k.','linewidth',6)
 	h.CapSize = 12;

    hold on 
    

    mpn = scatter(xval,mean(meanval2'),'r','LineWidth',6)

    hold on
    h = errorbar (xval,mean(meanval2'), mean(stdval2'),'r.','linewidth',6)
 	h.CapSize = 12;

    hold on 
    

    hpk = scatter(xval,mean(meanval3'),'b','LineWidth',6)

    hold on
    h = errorbar (xval,mean(meanval3'), mean(stdval3'),'b.','linewidth',6)
 	h.CapSize = 12;
end

	h.CapSize = 12;
    xlabel('Number of CpGs')
    ylabel('Energy, kT/bp')
    %title('Chromosome %i',ich)
    ytickformat('%,.1f')

    legend([cpg,mpn,hpk],'CpG','MpN','HpK')
    
	xlim([0,40])
	
	ylim([4.5,7.0])
grid on
    set(gca,'FontSize',36)
    ax = gca;
    box on
    grid on
         
     ft = 'Times'
fsz = 36;    

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)



   
total = sum(chromnum')



