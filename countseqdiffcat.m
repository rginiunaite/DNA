%% count the number of sequences in each category


i1=1;
i2=24;
k=1;
for ich = i1:i2

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpGintersectNMI147.mat',ich);

datagroup = load(filename);

first(ich,1)=length(datagroup.allsubsequences);

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%iNMINOTintersectCpG147.mat',ich);

datagroup = load(filename);

first(ich,2)=length(datagroup.allsubsequences);

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpGNOTintersectNMI147.mat',ich);

datagroup = load(filename);

first(ich,3)=length(datagroup.allsubsequences);

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%ineither147.mat',ich);

datagroup = load(filename);

first(ich,4)=length(datagroup.allsubsequences);


end


