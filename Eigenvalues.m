%% eigenvalues

paramset = load('ParameterSets/Di-hemi-meth-hmeth.mat');

% widom
%sequence = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT';

% widom centre 73

%sequence = 'GTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTC';

%sequence = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';

%sequence = 'ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA';

%sequence =  [random_dna(147)];
fileID = fopen('Sequences/Sequences1000CpG0.txt','r');

D = textscan(fileID,'%s');
fclose(fileID);

a = D{1,1};

for i=1:1000

    sequence = a{i};

    [shapesW, stiffW] = constructSeqParms(sequence, paramset);
    
    eigW = eig(stiffW);

    %minW = min(eigW)

    maxW = maxk(eigW,5);
    
    CpG0eigenvalues(:,i) = maxW;

    %meanW = mean(eigW)

end

save('CpG0eigenvalues.mat','CpG0eigenvalues')




% figure 
% 
% plot(eigW)
% 
% title('A')