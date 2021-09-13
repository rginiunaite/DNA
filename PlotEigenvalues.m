%% plot eigenvalues with different number of CpG islands


for i =1:5
        j = (i-1)*20;	
	eigenvalues = load('EigenvaluesData/CpG0eigenvalues.mat');
	a = eigenvalues.(1);
	mean(mean (a))

end
