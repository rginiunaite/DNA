% periodicity, frequency, autocorrelation, maybe fft

paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');
%seqnum = 20;

%from already generated files
% distances = load('distances100sequences.mat');
% 
% distances = cell2mat(struct2cell(distances)); %;zeros(seqnum,350);

%% select specific groups/organisms
vec = [data.Seq.group];
% 
str = 'Yeast';
idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
% idx2 = strcmp('Mouse',{data.Seq.group})
%idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth


% %% all sequences
% str = 'All';
% datagroup = data.Seq; %% all sequences

% %% random sequences
% 
%  str = 'Random';
% for j =1 :50
%     datagroup(j).Seq = [random_dna(400)];
% end

% %load 500 previously generated random sequences
% data500 = load('Random500Sequences400length.mat');
% 
% datagroup = data500.datagroup;

%% for all
seqnum = length(datagroup);

lenseq =73;

distances = zeros(seqnum,400-lenseq);


k=1; %count;
for j = 1:1

    seq = datagroup(j).S; % from nucleosomes data
    %seq = datagroup(j).Seq; % for random sequences
    
    
    if length(seq)==400
        
        [shapes, stiff] = constructSeqParms(seq, paramset);
        abs_coord = frames(shapes); % relative to absolute coordinates

        for i=1:400-lenseq 
            distances(k,i) = norm(abs_coord(i).rc -abs_coord(i+lenseq).rc);

        end
        k=k+1;
    end

end





for i =1:1
   figure
    x = distances(i,:);
    plot((x))
    %[ind(i),minim(i)]= min(x);
    [pks, locs] = min((x))
    
  %  GMModel = fitgmdist(x',2)
    title([str,num2str(i)])
end

cycles = diff(locs);

meancucle = mean(cycles);



% y = fft(x);
% 
% iy = ifft(y);
% 
% x2 = real(iy);
% 
% norm(x-x2);
% 
% dt = 1;
% et = 400-lenseq;
% t = 1:dt:et;
% y = distances(1,:);
% 
% figure 
% 
% subplot(2,1,1)
% plot(t,y); grid on
% 
% xlabel('sequence number')
% ylabel('Amplitude')
% 
% Y = fft(y);
% n = size(y,2)/2;
% amp_spec = abs(Y)/n;
% 
% subplot(2,1,2)
% freq = (0:79)/(2*n*dt);
% 
% plot(freq,amp_spec(1:80)); grid on
% 
% xlabel('Frequency, Hz')
% ylabel('Amplitude')
% ylim([0,1])



















% % figure
% % plot (distances(1,:));
% % 
% % f = fft(distances(1,:));
% % figure
% % plot(f)
% % ax = gca;
% % %ax.XLim = [0 10];
% % xlabel('Frequency (cycles/week)')
% % ylabel('FFT')
% 
% 
% figure
% ac=xcorr(distances(1,:),distances(1,:));
% [~,locs]=findpeaks(ac);
% mean(diff(locs)*0.1)
% 
% 
% 
% figure
% [autocor,lags] = xcorr(tempnorm,400-lenseq,'coeff');
% 
% plot(lags/fs,autocor)
% xlabel('Lag (days)')
% ylabel('Autocorrelation')
% 
% 
% [pksh,lcsh] = findpeaks(autocor);
% short = mean(diff(lcsh))/fs
% 
% [pklg,lclg] = findpeaks(autocor, ...
%     'MinPeakDistance',10,'MinPeakheight',0.1); %;ceil(short)*fs
% mean(diff(lclg))
% long = mean(diff(lclg))/fs
% 
% 
% hold on
% pks = plot(lags(lcsh)/fs,pksh,'or', ...
%     lags(lclg)/fs,pklg+0.05,'vk');
% hold off
% legend(pks,[repmat('Period: ',[2 1]) num2str([short;long],0)])

