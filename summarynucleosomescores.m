% nucleosome scores visualisation/summary


cpgdata = load('data/Yatzi_Harwood/chYCpGscores.csv'); % Yatzi
notcpgdata = load('data/Yatzi_Harwood/ch1YazdinotCpGscores.csv');      % Yatzi

%cpgdata = load('data/Yatzi_Harwood/ch10HarwoodCpGscores.csv'); % Harwood
%notcpgdata = load('data/Yatzi_Harwood/ch10HarwoodnotCpGscores.csv');    % Harwood

cpgmean = mean(cpgdata)
notcpgmean = mean(notcpgdata)

cpgstd = std(cpgdata)
notcpgstd = std(notcpgdata)
