%% load all chromosomes data and plot

chr1=1;
chr2=1;

for i=chr1:chr2

	xxxcpg = load('chr%isplinecpgxxx.mat',i);
	yyycpg = load('chr%isplinecpgyyy.mat',i);
	xxxneither = load('chr%isplineneitherxxx.mat',i);
	yyyneither = load('chr%isplineneitheryyy.mat',i);

end



