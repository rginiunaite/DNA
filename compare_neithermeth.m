%compare energies: CpG sequences
function u = compare_neithermeth(i1, i2)

for ich=i1:i2

filename = sprintf('../SkirmantasData/AllCentresCh37/allsubsequencesCh%ineithermeth147.mat',ich);

%filename = sprintf('SkirmantasData/AllChromosomesSimulations/coordinatesChYCpG.mat');

datagroup = load(filename);

seqnum = length(datagroup.data);

energies = zeros(1,length(datagroup.data));
k=1;
for j = 1:length(datagroup.data)

    seq = datagroup.data(j); % for random sequences
    seq = char(seq);
    %    seq = seq(3:end-2);
    %seq = seq(28:end-27);
    if(strlength(seq)>1)
                                if(seq~='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')
     energies(k) = seq2nucleosome_grad_descent_new2(seq);
     k=k+1;
   end
   end

end


%writematrix(energies, "energiesSeq10000CpG0.txt");


name = sprintf('../SkirmantasData/EnergiesCh37/energiesCh%ineithermeth147.txt',ich);
writematrix(energies, convertCharsToStrings(name));

end
