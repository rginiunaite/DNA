%% read and compare CpG and NMIs

fileID = fopen('SkirmantasData/NMILongPaper.bed','r');

D = textscan(fileID,'%s');
fclose(fileID);

a = D{1,1};

Ch1numberNMI = 4912;
coordinatesNMI = zeros(Ch1numberNMI,2);
k=1;
for i =1:2:Ch1numberNMI*2

	coordinatesNMI(k,:)= [str2num(cell2mat(a(2*i))),str2num(cell2mat(a(2*i+1)))];
	k=k+1;
end


%% CpGs

Ch1numberCpG = 2462;

fileID = fopen('SkirmantasData/dataCpGLong.csv','r');

D = textscan(fileID,'%s');
fclose(fileID);

a = D{1,1};

coordinatesCpG = zeros(Ch1numberCpG,2);
k=1;
for i =1:2:Ch1numberCpG*2

	coordinatesCpG(k,:)= [str2num(cell2mat(a(2*i))),str2num(cell2mat(a(2*i+1)))];
	k=k+1;
end





%figure 

%for i=1:100
%	hold on
%	plot([coordinatesNMI(i,1),coordinatesNMI(i,2)],[1,1],'r','LineWidth',10);
%	plot([coordinatesCpG(i,1),coordinatesCpG(i,2)],[1.1,1.1],'b','LineWidth',10);	
%end

%ylim([-5,5]);


%% find where CpG intersect NMIs (also NMI free CpGs)
coordinatesCpGintersectNMI = zeros(Ch1numberCpG,2);
coordinatesCpGNOTintersectNMI = zeros(Ch1numberCpG,2);
k=1;
c=1;
for j =1:Ch1numberCpG
	interval1 = fixed.Interval(coordinatesCpG(j,1),coordinatesCpG(j,2));
		found = 0; % check whether intersection was found
		for i =1:Ch1numberNMI
			interval2 = fixed.Interval(coordinatesNMI(i,1),coordinatesNMI(i,2));
			if overlaps(interval1, interval2) == 1
			
			
			intersection = intersect(interval1,interval2);
			coordinatesCpGintersectNMI(k,1) = intersection.LeftEnd;
			coordinatesCpGintersectNMI(k,2) = intersection.RightEnd;
			found =1;
			k=k+1;
			break	
					
			end		 
		end
		if found == 0 %% if no intersection found
			coordinatesCpGNOTintersectNMI(c,1) = interval1.LeftEnd;
			coordinatesCpGNOTintersectNMI(c,2) = interval1.RightEnd;
			c=c+1;
		end
end

save('coordinatesCpGintersectNMI.mat','coordinatesCpGintersectNMI') 
Data = load('coordinatesCpGintersectNMI.mat');                     
DataField = fieldnames(Data);
dlmwrite('coordinatesCpGintersectNMI.txt', Data.(DataField{1}));

save('coordinatesCpGNOTintersectNMI.mat','coordinatesCpGNOTintersectNMI') 
Data = load('coordinatesCpGNOTintersectNMI.mat');                     
DataField = fieldnames(Data);
dlmwrite('coordinatesCpGNOTintersectNMI.txt', Data.(DataField{1}));



%% find intervals with NMIs but not CpGs


coordinatesNMINOTintersectCpG = zeros(Ch1numberNMI,2);
k=1;
c=1;
for j =1:Ch1numberNMI
	interval1 = fixed.Interval(coordinatesNMI(j,1),coordinatesNMI(j,2));
		found = 0; % check whether intersection was found
		for i =1:Ch1numberNMI
			interval2 = fixed.Interval(coordinatesCpG(i,1),coordinatesCpG(i,2));
			if overlaps(interval1, interval2) == 1
			intersection = intersect(interval1,interval2);
			found =1;
			k=k+1;
			break	
					
			end		 
		end
		if found == 0 %% if no intersection found
			coordinatesNMINOTintersectCpG(c,1) = interval1.LeftEnd;
			coordinatesNMINOTintersectCpG(c,2) = interval1.RightEnd;
			c=c+1;
		end
end


save('coordinatesNMINOTintersectCpG.mat','coordinatesNMINOTintersectCpG') 
Data = load('coordinatesNMINOTintersectCpG.mat');                     
DataField = fieldnames(Data);
dlmwrite('coordinatesNMINOTintersectCpG.txt', Data.(DataField{1}));
