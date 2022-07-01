% I will only look at chromosome 1 and will check whether the number of CpG dinucleotides correlate with the energy

clear all
filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh5CpG147.mat');

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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh5CpG147.txt');

energiesCpG = load(filename);

cor = corrcoef(dimersALL(1:end-1),nonzeros(energiesCpG)); % the last one is just zeros, so I remove it.


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh5notCpG147.mat');

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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh5notCpG147.txt');

energiesnotCpG = load(filename);

cor = corrcoef(dimersALLnot(1:end-1),nonzeros(energiesnotCpG)); % the last one is just zeros, so I remove it.


figure

a1= scatter(dimersALL(1:end-1),nonzeros(energiesCpG/147),'filled','k');
x=dimersALL(1:end-1);
x=[x,dimersALLnot(1:end-1)];
y1=nonzeros(energiesCpG/147);
y1=[y1',nonzeros(energiesnotCpG/147)'];
hold on
a21 = scatter(dimersALLnot(1:end-1),nonzeros(energiesnotCpG/147),'filled','k');



%x2=dimersALLnot(1:end-1);

%y2=nonzeros(energiesnotCpG/147);


P = polyfit(x,y1,1);
yfit = polyval(P,x);
hold on;
plot(x,yfit,'k-.','LineWidth',3);

%P = polyfit(x2,y2,1);
%yfit = polyval(P,x2);
%hold on;
%plot(x2,yfit,'k-.','LineWidth',3);

grid on
%legend('CpG','not CpG')

xlabel('Number of CpGs')
ylabel('Energy, kT/bp')
%title('Not CpG islands')
%ytickformat('%.1f')
%     
    set(gca,'FontSize',36)
    ax = gca;
    box on
         
     ft = 'Times'
fsz = 36;    

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)



%% MpN

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh5CpG147.mat');

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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh5MpN147.txt');

energiesCpG = load(filename);

cor = corrcoef(dimersALL(1:end-1),nonzeros(energiesCpG)); % the last one is just zeros, so I remove it.


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh5notCpG147.mat');

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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh5notMpN147.txt');

energiesnotCpG = load(filename);

%cor = corrcoef(dimersALLnot(1:end-1),nonzeros(energiesnotCpG)); % the last one is just zeros, so I remove it.


hold on

a1= scatter(dimersALL(1:end-1),nonzeros(energiesCpG/147),'filled','r');
x=dimersALL(1:end-1);
x=[x,dimersALLnot(1:length(nonzeros(energiesnotCpG/147)))];
y1=nonzeros(energiesCpG/147);
y1=[y1',nonzeros(energiesnotCpG/147)'];
hold on
a22 = scatter(dimersALLnot(1:length(nonzeros(energiesnotCpG/147))),nonzeros(energiesnotCpG/147),'filled','r');

x2=dimersALLnot(length(nonzeros(energiesnotCpG/147)));

y2=nonzeros(energiesnotCpG/147);


P = polyfit(x,y1,1);
yfit = polyval(P,x);
hold on;
plot(x,yfit,'r-.','LineWidth',3);

%% HpK


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh5CpG147.mat');

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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh5HpK147.txt');

energiesCpG = load(filename);

cor = corrcoef(dimersALL(1:end-1),nonzeros(energiesCpG)); % the last one is just zeros, so I remove it.


filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh5notCpG147.mat');

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


filename = sprintf('../SkirmantasData/EnergiesCh37/energiesCh5notHpK147.txt');

energiesnotCpG = load(filename);

%cor = corrcoef(dimersALLnot(1:end-1),nonzeros(energiesnotCpG)); % the last one is just zeros, so I remove it.


hold on

a1= scatter(dimersALL(1:end-1),nonzeros(energiesCpG/147),'filled','b');
x=dimersALL(1:end-1);
x=[x,dimersALLnot(1:length(nonzeros(energiesnotCpG/147)))];
y1=nonzeros(energiesCpG/147);
y1=[y1',nonzeros(energiesnotCpG/147)'];
hold on
a23 = scatter(dimersALLnot(1:length(nonzeros(energiesnotCpG/147))),nonzeros(energiesnotCpG/147),'filled','b');

x2=dimersALLnot(1:end-1);

y2=nonzeros(energiesnotCpG/147);




P = polyfit(x,y1,1);
yfit = polyval(P,x);
hold on;
plot(x,yfit,'b-.','LineWidth',3);
legend([a21,a22,a23],'CpG','MpN','HpK')

ytickformat('%,.1f')

