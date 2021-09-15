%% plot eigenvalues with different number of CpG islands

avvalue = zeros(1,5);
avstd = zeros(1,5);
      	
% eigenvalues = load('EigenvaluesData/CpG0eigenvalues.mat');
% a = eigenvalues.CpG0eigenvalues(1,:);
% avvalue(1) = mean(a);
% avstd(1) = std(a);
% 
% 
% eigenvalues = load('EigenvaluesData/CpG20eigenvalues.mat');
% a = eigenvalues.CpG0eigenvalues(1,:);
% avvalue(2)=mean(a);
% avstd(2) = std(a);
% 
% eigenvalues = load('EigenvaluesData/CpG40eigenvalues.mat');
% a = eigenvalues.CpG0eigenvalues(1,:);
% avvalue(3) = mean(a);
% avstd(3) = std(a);
% 
% eigenvalues = load('EigenvaluesData/CpG60eigenvalues.mat');
% a = eigenvalues.CpG0eigenvalues(1,:);
% avvalue(4) = mean(a);
% avstd(4) = std(a);

eigenvalues = load('EigenvaluesData/CpG80eigenvalues.mat');
a = eigenvalues.CpG0eigenvalues(1,:);
avvalue(5) = mean(a);
avstd(5) = std(a);

figure 
x= [0,20,40,60,80]
plot(x, avvalue,'LineWidth',2)

hold on 
errorbar(x,avvalue, avstd,'k.','linewidth',2)
set(gca,'FontSize',36)
ax=gca;
box on
grid on
xlabel('Number of CpGs')
ylabel('Maximum eigenvalue')




