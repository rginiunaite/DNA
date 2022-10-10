%my_file = fopen('../SkirmantasData/Harwood2019/Ch1Harwooddata.csv');
%my_file = fopen('../SkirmantasData/Harwood2019/Ch16Harwooddata.csv');

%my_file = fopen('../SkirmantasData/Yazdi2015/chr1_yatsi_H1_dyad_withNOS.csv');

chrnr = 21;

filename = sprintf('../SkirmantasData/Yazdi2015/chr%iYazdi.csv',chrnr);
my_file = fopen(filename);
%my_file = fopen('../SkirmantasData/Shep2015/GM.nucpos.csv');

%my_file = fopen('../SkirmantasData/Harwood2019/GSM3312888_iPS_chr1-22_138-161bp.csv');



b = fscanf(my_file,'chr%d %d %d %f \n', [3,Inf]);% shep
size(b)

fclose(my_file);
clear my_file
b = b';
N = b(end,2);


ba = zeros(N,1);
for i = 1:length(b)
   ba(b(i,2)) = b(i,3);
   
end

%ba(ba>100) = 100;
filename = sprintf('../SkirmantasData/Ch%iNMILongPaper.csv',chrnr);
nmi_file = fopen(filename);
%nmi_file = fopen('../SkirmantasData/Harwood2019/NMILongPaper.csv');



%filename = sprintf('chr%i,%d,%d,%d \n',chrnr);


%nma = fscanf(nmi_file,'chr22,%d,%d,%d \n', [3,Inf]);
nma = fscanf(nmi_file,'chr%d,%d,%d,%d \n', [3,Inf]);
fclose(nmi_file);
clear nmi_file
nma = nma';
filename =sprintf('../SkirmantasData/Ch%iCGILong.csv',chrnr);
cgi_file = fopen(filename);

%filename = sprintf('chr22,%d,%d,CpG:_%d\n',chrnr);

%cgi_file = fopen('../SkirmantasData/Harwood2019/Ch16CGILong.csv');
cga = fscanf(cgi_file,'chr%d,%d,%d,CpG:_%d\n', [3,Inf]);
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

nonzeros_meancpgnmi = mean(nonzeros(ba(logical(cgn))))
nonzeros_meanneither = mean(nonzeros(ba(logical(neither))))

nonzeros_meanonlycpg = mean(nonzeros(ba(logical(onlycpg))))
nonzeros_meanonlynmi = mean(nonzeros(ba(logical(onlynmi))))

meancpgnmi = mean((ba(logical(cgn))))
meanneither = mean((ba(logical(neither))))

meanonlycpg = mean(ba(logical(onlycpg)))
meanonlynmi = mean(ba(logical(onlynmi)))


%std(ba(logical(cgn)))
%std(ba(logical(neither)))
%std(ba(logical(cgn)))
%std(ba(logical(neither)))

bacpg = ba(logical(cgn));
neitherba = ba(logical(neither));


%save('Ch3cpg_nucl_scores.mat','bacpg')
%save('Ch3neither_nucl_scores.mat','neitherba')

