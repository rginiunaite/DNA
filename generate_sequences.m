
% R is value specifying the minimum CpGobserved/CpGexpected ratio 
% in each window needed to mark a base. 
% Choices are a value between 0 and 1. Default is 0.6. This ratio is defined as:
% CPGobs/CpGexp = (NumCpGs*Length)/(NumGs*NumCs)

N = 200;
SeqN = 100;

filename = 'SequencesCpG20Trial.txt';  % better to use fullfile(path,name) 


for count = 1:SeqN
S = randseq(N);
R =(dimercount(S).CG*length(S))/(basecount(S).C*basecount(S).G);

CpGnum = 20; 
CpGcount = dimercount(S).CG;
while CpGcount < CpGnum
    
  i = round(rand*N);
  
  if ((i<N-1)&&(i>1)&&(dimercount(S(i-1:i+2)).CG == 0))     %(~strcmpi(S(i:i+1),'CG')))
      
      S(i:i+1) = 'CG';
      CpGcount = CpGcount + 1;
      
  end    

end    

Check = dimercount(S).CG;  %testing

    fid = fopen(filename,'a');
    fprintf(fid,'%s\n',S);          % Write the char array, interpret newline as new line
    fclose(fid);  

end
%CpGcount

%-------------------------------
% CpG islands: real and random
% Histograms of length

if 0 %not active

S1 = randseq(248956422); %249250621 ?
cpgStruct = cpgisland(S1);
L = cpgStruct.Stops - cpgStruct.Starts;
max(L)
hist(L)

sd = fastaread('sequence.fasta');
S2 = sd.Sequence;
cpgStruct = cpgisland(S2);
L = cpgStruct.Stops - cpgStruct.Starts;
max(L)
hist(L)

end



