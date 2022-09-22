cpgload = load('../SkirmantasData/Harwood2019/Ch3cpg_nucl_scores.mat') 

cpgdata= cpgload.bacpg;


neitherload = load('../SkirmantasData/Harwood2019/Ch3neither_nucl_scores.mat')

neitherdata = neitherload.neitherba;


fig=figure; set(fig,'visible','off');
h = histogram(cpgdata,'Normalization','probability');
edges = h.BinEdges;  
values = h.Values; 
close(fig);

xx= 0:1:edges(end-1) ;
yy = spline(edges(1:end-1), values(1:end),xx)  ;

xxx=0:0.1:edges(end-1) ;
yyy = spline(xx, yy,xxx);  

figure              
a=plot(xxx,yyy,'k','LineWidth',2)           

fig=figure; set(fig,'visible','off');
h = histogram(neitherdata,'Normalization','probability');
edges = h.BinEdges;  
values = h.Values; 
close(fig);
  
xx= 0:1:edges(end-1) ;
yy = spline(edges(1:end-1), values(1:end),xx);  

xxx=0:0.1:edges(end-1) ;
yyy = spline(xx, yy,xxx);  

hold on           
b=plot(xxx,yyy,'g','LineWidth',2)  ;
xlim([0,6]) 

 xlabel('Nucleosome mid-point scores')
 ylabel('Normalised frequency')
    
        
    set(gca,'FontSize',36)
ax=gca;
ax.BoxStyle = 'full';
ax.LineWidth = 3;
         
     ft = 'Times'
fsz = 36;    

set(findall(gcf,'type','text'), 'FontSize', fsz, 'Color', 'k','FontName', ft)

ytickformat('%.1f')

legend([a,b],'CGI and NMI','Not CGI and not NMI')
    


