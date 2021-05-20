%% plot persistence length
str = 'Rat and mouse';

for i = 1:10 % HUman 7, Drosophila 21, Rat and mouse 10
    
    filename = sprintf('RatMousePl/PlApp%i.txt',i);
    data = load(filename);
    plApp(:,i) = data; % columns different sequences, rows, different positions in a sequence
    
    filename = sprintf('RatMousePl/PlDyn%i.txt',i);
    data = load(filename);
    plDyn(:,i) = data; % columns different sequences, rows, different positions in a sequence
end



%x = [1,50,100,150,200,250,300];
x=1:25:327;
%figure

figure
plot1 = plot(x,(mean(plApp')),'b','LineWidth',3)

y = mean(plApp');

sd_vct = std(plApp');
hold on
%errorbar(x,y,sd_vct/sqrt(length(persapp)), 'b','LineWidth',3)
shadedErrorBar(x,y,sd_vct/sqrt(length(plApp(:,1))));

%hold on
plot2 = plot(x,(mean(plDyn')),'r','LineWidth',3);

y = mean(plDyn');

sd_vct = std(plDyn');
%hold on
%errorbar(x,y,sd_vct/sqrt(length(persdyn)), 'r','LineWidth',3')
shadedErrorBar(x,y,sd_vct/sqrt(length(plApp(:,1))));

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