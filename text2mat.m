%% all from txt to mat files


for ich=1:24

name = sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%iMpNNOTintersectNMI147.txt',ich);

newname =  sprintf('SkirmantasData/AllCentresCh37/allsubsequencesCh%iMpNNOTintersectNMI147.mat',ich);
data=importdata(name);
save(newname, 'data');


end
