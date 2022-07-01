%% count the number of sequences with different number of CpG dinucleotides


i1=1;
i2=24;
k=1;
for ich = i1:i2

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

%x=dimersALL(1:end-1);
x = dimersALLnot(1:end-1);
%x=[x,dimersALLnot(1:end-1)];

k1=1;
k2=1;
k3=1;
k4=1;
k5=1;
k6=1;
k7=1;
for i=1:length(x)

if x(i)<5

	k1=k1+1;
elseif x(i)>=5 && x(i)<15

	k2=k2+1;
elseif x(i)>=15 && x(i) < 25

	k3=k3+1;
elseif x(i)>=25 

	k4=k4+1;
end
end

% check how many sequences with different number of CpGs
chromnum(1,ich)= k1;
chromnum(2,ich)= k2;
chromnum(3,ich)= k3;
chromnum(4,ich)= k4;



end


