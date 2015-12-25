%Prelim figures

%Glacier bed
%helheim_fronts

bed=runmean(flowline(1:344,4),5);
surf=runmean(flowline(1:344,5),5);

flow=flowline(1:344,1);
[C,IA,IC]=unique(flow);
flow=flow(IA);
bed=bed(IA);
surf=surf(IA);

bed_fronts=interp1(flow,bed,fronts(:,2));
bed_rifts=interp1(flow,bed,rifts(:,2));

figure;
plot(flow/10^3,surf,'k','linewidth',2); hold on;
plot([flowline(344,1)/10^3 flowline(344,1)/10^3],[bed(end) surf(end)],'k','linewidth',2)
plot(flow/10^3,bed,'color',[0.55 0.27 0.07] ,'linewidth',2); hold on;
plot([flow(end)/10^3 flow(end)/10^3+5],[0 0],'b','linewidth',2);
plot([flow(end)/10^3 flow(end)/10^3+5],[bed(end) bed(end)],'--','color',[0.55 0.27 0.07],'linewidth',2);
xlim([0 flow(end)/10^3+5])
ylabel('Elevation (m)','fontsize',16,'fontname','arial');
xlabel('Distance along flowline (km)','fontsize',16,'fontname','arial');
ylim([-1200 1400]);
set(gca,'fontsize',12);
plot([36 42 42 36 36],[-1000 -1000 400 400 -1000],'k--')
%plot(fronts(:,2)/10^3,bed_fronts,'bo');

figure;
plot(flow/10^3,surf,'k','linewidth',2); hold on;
plot([flowline(344,1)/10^3 flowline(344,1)/10^3],[bed(end) surf(end)],'k','linewidth',2)
plot(flow/10^3,bed,'color',[0.55 0.27 0.07] ,'linewidth',2); hold on;
plot([flow(end)/10^3 flow(end)/10^3+5],[0 0],'b','linewidth',2);
plot([flow(end)/10^3 flow(end)/10^3+5],[bed(end) bed(end)],'--','color',[0.55 0.27 0.07],'linewidth',2);
xlim([32 42]);
ylim([-1000 400]);
ylabel('Elevation (m)','fontsize',16,'fontname','arial');
xlabel('Distance along flowline (km)','fontsize',16,'fontname','arial');
set(gca,'fontsize',12);
h(1)=plot(fronts(:,2)/10^3,bed_fronts,'ko','markersize',10,'markerfacecolor','k');
h(2)=plot(rifts(:,2)/10^3,bed_rifts,'ro','markersize',10,'markerfacecolor','r');


cd ~/Data/Mosaics/Helheim/
%dimensions=[2.94*10^5 3.15*10^5 -2.57*10^6 -2.586*10^6];
%I=geotiffread('mosaicHelheim.2013-205.148.33882_1-20mgeo.tif',dimensions);
%figure;
%imshow(I.z,'xdata',I.x,'ydata',I.y); hold on;
%set(gca,'ydir','normal');
%plot([dimensions(2)-5700,dimensions(2)-700],[dimensions(4)+1400 dimensions(4)+1400],'k','linewidth',2)
%text(dimensions(2)-4000,dimensions(4)+800,'5 km','fontsize',24)
%plot(flowline(:,2),flowline(:,3),'k','linewidth',3);



