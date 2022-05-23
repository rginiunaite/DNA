%% all from mat to txt files

for ich=1:24

name = sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpGNOTintersectNMI147.mat',ich);


Data=load(name);
DataField = fieldnames(Data);
%writematrix(Data.(DataField{1}),'SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpG147.txt',ich)
newname =  sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%iCpGNOTintersectNMI147.txt',ich);

writematrix(Data.(DataField{1}),newname)


end
