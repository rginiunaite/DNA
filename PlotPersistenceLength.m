%% plot persistence length
str = 'Rats and mice';
persistence = load('RatsMiceFull500Every25.mat');
%persistence = load('RandomFull100.mat');

persapp= persistence.apparent;

persdyn = persistence.dynamic;


%x = [1,50,100,150,200,250,300];
x=1:25:327;
%figure

% plot(x,mean(persapp))
% hold on
% plot(x,mean(persdyn))


figure
plot1 = plot(x,(mean(persapp)),'b','LineWidth',3)

y = mean(persapp);

sd_vct = std(persapp);
hold on
%errorbar(x,y,sd_vct/sqrt(length(persapp)), 'b','LineWidth',3)
shadedErrorBar(x,y,sd_vct/sqrt(length(persapp(:,1))));

%hold on
plot2 = plot(x,(mean(persdyn)),'r','LineWidth',3);

y = mean(persdyn);

sd_vct = std(persdyn);
%hold on
%errorbar(x,y,sd_vct/sqrt(length(persdyn)), 'r','LineWidth',3')
shadedErrorBar(x,y,sd_vct/sqrt(length(persapp(:,1))));

set(gca,'linew',3)
ylabel('Persistence length')
%ytickformat('%.1f')
xlabel('Starting position')
title([str])
xlim([0,327])
%set(gca, 'YTick', round(1:4))

legend([plot1, plot2],'Apparent','Dynamic')

set(gca,'FontSize',36)
ax = gca;
grid on
box on