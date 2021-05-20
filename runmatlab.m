%% compute persistence lengths
function [pl_app,pl_dyn] = runmatlab(seqnr) 

tang1 = "example_bash_new/Output_t0";
tang2 = "example_bash_new/Output_t0_intr";

addpath('ttc_output')

filet0 = fopen(tang1,'r');
sizeA = [2 Inf];
A1 = fscanf(filet0,'%d %f',sizeA);

t0 = A1(2,:);

nbp = length(t0);

filet0intr = fopen(tang2,'r');
A2 = fscanf(filet0intr,'%d %f',sizeA);
t0intr = A2(2,:);


[pl_app,pl_dyn,~] = Compute_pl(nbp,t0,t0intr);

%filename='matlaboutput/plapp.mat';
%save(filename,'pl_app','-append');

textfile = sprintf('matlaboutput/YeastPl/PlApp%i.txt',seqnr);
dlmwrite(textfile,pl_app,'-append','newline','pc')

%fileID = fopen('plapp.txt','w');
%fprintf(fileID, '%d\n',pl_app);
%fclose(fileID);

%filename='matlaboutput/pldyn.mat';
%save(filename,'pl_dyn','-append');
textfile2 = sprintf('matlaboutput/YeastPl/PlDyn%i.txt',seqnr);

dlmwrite(textfile2,pl_dyn,'-append','newline','pc')


end





