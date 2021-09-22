%% compare PL lengths



AppCpGs = load('PLdata/CpG10000output/PlCpG0t0App1.txt');
AppMpNs = load('PLdata/MpN10000output/PlCpG0t0App1.txt');


DynCpGs = load('PLdata/CpG10000output/PlCpG0t0Dyn1.txt');
DynMpNs = load('PLdata/MpN10000output/PlCpG0t0Dyn1.txt');

l = length(AppCpGs);

StatCpGs = zeros(1,l);
StatMpNs = zeros(1,l);

for i=1:l
    StatCpGs(i) = 1/(1/AppCpGs(i) - 1/DynCpGs(i));
    StatMpNs(i) = 1/(1/AppMpNs(i) - 1/DynMpNs(i));
end


%% diff dynamic - apparent 

diffCpGDynminusApp = DynCpGs - AppCpGs;

diffMpNDynminusApp = DynMpNs - AppMpNs;

figure

histogram(diffCpGDynminusApp)
title('CpG difference (Dyn-App)')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on

figure

histogram(diffMpNDynminusApp)
title('MpN difference (Dyn-App)')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on



%% plot CpGs

figure 
histogram(DynCpGs)
title('CpG dynamic PL')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on


figure 
histogram(AppCpGs)
title('CpG apparent PL')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on

figure 
histogram(StatCpGs)
title('CpG static PL')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on


%% plot MpNs

figure 
histogram(DynMpNs)
title('MpN dynamic PL')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on


figure 
histogram(AppMpNs)
title('MpN apparent PL')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on

figure 
histogram(StatMpNs)
title('MpN static PL')
xlabel('Persistence length')
ylabel('Frequency')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on

