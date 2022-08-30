%% plot nucleosome scores and energies



Table = readtable('../SkirmantasData/Yazdi2015/chr1_myseq_yatsi_dyad.csv'); % data from that specific location


energiesdata = load('data/Yatzi_Harwood/human_200k_u.mat');


figure



%scatter(Table.Var4,energiesdata.utot(Table.Var2-800001+73)) % add 73, since nucleosome score is for centres, energies for start


hold on 
plot(800001:length(energiesdata.utot)+800000,energiesdata.utot-mean(energiesdata.utot),'g') % add 73, since nucleosome score is for centres, energies for start
for i=1:length(Table.Var1)
hold on
plot([Table.Var2(i)-73,Table.Var3(i)+74],[Table.Var4(i)/max(Table.Var4),Table.Var4(i)/max(Table.Var4)],'b','LineWidth',5)
end
hold on
scatter(Table.Var2,energiesdata.utot(Table.Var2-800001+73)) % add 73, since nucleosome score is for centres, energies for start

% Compute FFT for energies

X = energiesdata.utot;
Fs = length(X);            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(X);             % Length of signal
t = (0:L-1)*T;
Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
f = Fs*(0:(L/2))/L;
plot(f,P1) 
ylim([0,0.01])
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')



% Compute FFT for nucleosome scores
aint = interp1(Table.Var2-800001,Table.Var4,1:10:length(X),'spline'); %% interpolate data

X1 = aint;
Fs = length(aint);            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(aint);             % Length of signal
t = (0:L-1)*T;
Y = fft(X1);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%ylim([0,0.01])



