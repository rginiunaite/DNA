%my_file = fopen('../SkirmantasData/Harwood2019/Ch1Harwooddata.csv');
my_file = fopen('../SkirmantasData/Harwood2019/Ch16Harwooddata.csv');
%my_file = fopen('../SkirmantasData/Harwood2019/GSM3312888_iPS_chr1-22_138-161bp.csv');

b = fscanf(my_file,'chr16 %d %d \n', [2,Inf]);
fclose(my_file);
clear my_file
b = b';
N = 10*length(b);


ba = zeros(N,1);
for i = 1:length(b)
   ba(b(i,1)+1:b(i,1)+10) = b(i,2);
   %ba(b(i,1)+1:b(i,1)) = b(i,2);
end

ba(ba>100) = 100;

%nmi_file = fopen('../SkirmantasData/NMILongPaper.csv');
nmi_file = fopen('../SkirmantasData/Harwood2019/Ch16NMILongPaper.csv');
nma = fscanf(nmi_file,'chr16,%d,%d,%d \n', [3,Inf]);
fclose(nmi_file);
clear nmi_file
nma = nma';

%cgi_file = fopen('../SkirmantasData/dataCpGLong.csv');
cgi_file = fopen('../SkirmantasData/Harwood2019/Ch16CGILong.csv');
cga = fscanf(cgi_file,'chr16,%d,%d,CpG:_%d\n', [3,Inf]);
fclose(cgi_file);
clear cga_file
cga = cga';
cga = cga(:,1:2);


cgai = zeros(N,1);
for i = 1:length(cga)
   cgai(cga(i,1):cga(i,2)) = 1;
end  

nmai = zeros(N,1);
for i = 1:length(nma)
   nmai(nma(i,1):nma(i,2)) = 1;
end
nmai = nmai(1:N);
cgai = cgai(1:N);

cgn = cgai.*nmai;

%ii = 1:n;
%jj(:,1) = find((cgn(2:end) - cgn(1:end-1))>0)+1;
%jj(:,2) = find((cgn(2:end) - cgn(1:end-1))<0);
%cNMI = jj;

neither = (~cgai).*(~nmai);
%njj(:,1) = find((cgn(2:end) - cgn(1:end-1))>0)+1;
%njj(:,2) = find((cgn(2:end) - cgn(1:end-1))<0);


onlycpg = (cgai).*(~nmai);
onlynmi = (~cgai).*(nmai);

%vidurkiai
meancpgnmi = mean(ba(logical(cgn)))
meanneither = mean(ba(logical(neither)))

meanonlycpg = mean(ba(logical(onlycpg)))
meanonlynmi = mean(ba(logical(onlynmi)))
std(ba(logical(cgn)))
std(ba(logical(neither)))
std(ba(logical(cgn)))
std(ba(logical(neither)))

bacpg = ba(logical(cgn));
neitherba = ba(logical(neither));


%save('Ch3cpg_nucl_scores.mat','bacpg')
%save('Ch3neither_nucl_scores.mat','neitherba')

