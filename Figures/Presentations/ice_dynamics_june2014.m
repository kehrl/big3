%Let's look at only 2012-2013.
only_2012=0;
if only_2012
    ind1 = find(floor(time)==2012);
    ind2 = find(floor(fronts(:,1))==2012);
else
    ind1=find(floor(time) >= 2000);
    ind2=find(floor(fronts(:,1)) >= 2000);
end

%Make plots?
plot_vel_term1=1;
plot_vel_term2=0;
plot_flowline=0;


%Plot velocity vs. terminus position
fronts_interp2=interp1(fronts(ind2,1),fronts(ind2,2),time(ind1));
for i=1:length(ind1)
    [dist(i),lag_5km(i)]=min(abs(flowline(:,1)-(fronts_interp2(i)-5*10^3)));
    [~,eul_5km(i)]=min(abs(flowline(:,1)-(flowline(end,1)-5*10^3)));
end
if plot_vel_term1
figure;
set(gcf,'PaperUnits','centimeters');
    xSize = 17; ySize = 9.5;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[750 4 xSize*50 ySize*50])
    
    x = 3/xSize;
    y = 1/ySize; 
    width = 1-x;
    height = 1-y;
    axes('position',[x y width height])

subplot(2,1,1);
plot(fronts(ind2,1),(fronts(ind2,2)-mean(fronts(ind2,2)))/10^3,'-','color',[0.8 0.8 0.8],'linewidth',1.5); hold on;
plot(fronts(ind2,1),(fronts(ind2,2)-mean(fronts(ind2,2)))/10^3,'ko','markerfacecolor','k','markersize',2.5);
%plot(rifts(:,1),(rifts(:,2)-flowline(344,1))/10^3,'ro','markerfacecolor','r','markersize',4);
xlim([2000 2014])
ylim([-2 5])
%legend('Terminus','Rift');
ylabel('Terminus position (km)','fontsize',9,'fontname','arial');
set(gca,'fontsize',9,'fontname','arial');
set(gca,'xticklabel',{});
xlabel('Year','fontsize',9,'fontname','arial');
set(gca,'ticklength',[0.025 0.025])
set(gca,'ytick',-2:2:4);
text(2000.5,3.7,'a','fontsize',9','fontname','arial','fontweight','bold');

h=subplot(2,1,2);
nonnan = find(~isnan(velocity(eul_5km(1),ind1)));
plot(time(ind1(nonnan)),velocity(eul_5km(1),ind1(nonnan))/10^3,'-','color',[0.8 0.8 0.8],'linewidth',1.5); hold on;
plot(time(ind1(nonnan)),velocity(eul_5km(1),ind1(nonnan))/10^3,'ko','markerfacecolor','k','markersize',2.5);
xlim([2000 2014])
xlabel('Year','fontsize',9,'fontname','arial');
ylabel('Glacier velocity (km/yr)','fontsize',9,'fontname','arial');
set(gca,'fontsize',9,'fontname','arial');
p=get(h,'position')
set(h,'position',[p(1) p(2)+0.10 p(3) p(4)])
set(gca,'ticklength',[0.025 0.025])
ylim([5.2 8.2])
set(gca,'ytick',[6:1:8])
text(2000.5,7.7,'b','fontsize',9','fontname','arial','fontweight','bold');
end

if plot_vel_term2
n=1; clear summer
for i=1:length(ind1)
    time1=time(ind1(i))-floor(time(ind1(i)));
    if time1 > 152/365 && time1 < 274/365
        summer(n)=i;
        n=n+1;
    end
end
figure;
set(gcf,'PaperUnits','centimeters');
    xSize = 9.3; ySize = 5;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[750 4 xSize*50 ySize*50])
    
    x = 3/xSize;
    y = 1/ySize; 
    width = 1-x;
    height = 1-y;
    axes('position',[x y width height])
subplot(1,2,1)
fronts_interp=interp1(fronts(ind2,1),fronts(ind2,2)-flowline(344,1),time(ind1));
for j = 1:length(ind1)
    h=plot(fronts_interp(j)/10^3,velocity(eul_5km(j),ind1(j))/10^3,'ko','markerfacecolor','b','markersize',4); hold on;
end
for j = 1:length(summer)
    g=plot(fronts_interp(summer(j))/10^3,velocity(eul_5km(summer(j)),ind1(summer(j)))/10^3,'ko','markerfacecolor','r','markersize',4);
end
p=polyfit(fronts_interp/10^3,velocity(eul_5km(i),ind1)/10^3,1);
plot(fronts_interp/10^3,fronts_interp/10^3*p(1)+p(2),'k','linewidth',1);
%text(-0,8.2,'y=-0.65x+6.7','fontsize',9)
text(-0,8.2,'R^2=0.50','fontsize',9)
xlabel('Terminus position (km)','fontsize',9,'fontname','arial');
ylabel('Eulerian velocity (km/yr)','fontsize',9,'fontname','arial');
text(-1.6,6.25,'a','fontsize',9,'fontweight','bold','fontname','arial');
ylim([6 8.5])
axis square
set(gca,'ticklength',[0.05 0.05])

subplot(1,2,2);
for j = 1:length(ind1)
    h=plot(fronts_interp(j)/10^3,velocity(lag_5km(j),ind1(j))/10^3,'ko','markerfacecolor','b','markersize',4); hold on;
end
for j = 1:length(summer)
    g=plot(fronts_interp(summer(j))/10^3,velocity(lag_5km(summer(j)),ind1(summer(j)))/10^3,'ko','markerfacecolor','r','markersize',5);
end
xlabel('Terminus position (km)','fontsize',9,'fontname','arial');
ylabel('Lagrangian velocity (km/yr)','fontsize',9,'fontname','arial');
ylim([6 8.5])
%legend([h(1) g(1)],'Oct.-May','June-Sept.','location','northeast');
for i=1:length(ind1)
    temp(i)=velocity(lag_5km(i),ind1(i));
end
p=polyfit(fronts_interp/10^3,temp/10^3,1);
%clear temp
plot(fronts_interp/10^3,fronts_interp/10^3*p(1)+p(2),'k','linewidth',1);
%text(-1,8.2,'y=-0.36x+6.6','fontsize',9)
text(-0,8.2,'R^2=0.26','fontsize',9)
text(-1.8,6.25,'b','fontsize',9,'fontweight','bold','fontname','arial');
axis square
set(gca,'ticklength',[0.05 0.05])

end

if plot_flowline
    figure;
    set(gcf,'PaperUnits','centimeters');
    xSize = 9; ySize = 10;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[750 4 xSize*50 ySize*50])
    
    x = 3/xSize;
    y = 1/ySize; 
    width = 1-x;
    height = 1-y;
    axes('position',[x y width height])
    
    bed=runmean(flowline(1:344,4),5);
    surf=runmean(flowline(1:344,5),5);
    bed(345:355)=bed(344);
    surf(345:355)=surf(344);
    
    flow=flowline(:,1);
    [C,IA,IC]=unique(flow);
    flow=flow(IA);
    bed=bed(IA);
    surf=surf(IA);
    [~,front1]=min(abs(min(fronts(:,2))-flowline(:,1)));
    [~,front2]=min(abs(max(fronts(:,2))-flowline(:,1)));
    
    for i = 1:length(flowline)
        nonnan = find(~isnan(velocity(i,:)));
        ave_vel(i) = mean(velocity(i,nonnan));
        std_vel(i) = std(velocity(i,nonnan));
    end
    
    s=subplot(2,1,1); 
    plot((flowline(:,1)-flowline(344,1))/10^3,velocity/10^3,'color',[0.7 0.7 0.7],'linewidth',1); hold on;
    plot((flowline(:,1)-flowline(344,1))/10^3,ave_vel/10^3,'k','linewidth',1.5);
    xlim([-32 5])
    ylim([2 11])
    set(gca,'xticklabel',{})
    ylabel('Velocity (km/yr)','fontsize',9,'fontname','arial');
    text(-30.5,10,'a','fontsize',9,'fontname','arial','fontweight','bold');
    set(gca,'ticklength',[0.025 0.025])
    set(gca,'ticklength',[0.025 0.025])
    p=get(s,'position');
    set(s,'position',[p(1)+0.05 p(2) p(3) p(4)])
    
    s=subplot(2,1,2);
    plot((flow(1:344)-flowline(344,1))/10^3,bed(1:344),'color',[0.55 0.27 0.07],'linewidth',1.5); hold on;
    plot((flow(1:front2)-flowline(344,1))/10^3,surf(1:front2),'k','linewidth',1.5);
    plot([(flow(front2)-flowline(344,1))/10^3 5],[0 0],'b','linewidth',1.5);
    %plot([0 0],[bed(end) surf(end)],'k','linewidth',1.5);
    h(2)=plot((flowline(min(lag_5km):max(lag_5km),1)-flowline(344,1))/10^3,surf(min(lag_5km):max(lag_5km)),'r','linewidth',1.5);
    h(1)=plot((flowline(eul_5km(1),1)-flowline(344,1))/10^3,surf(eul_5km(1)),'ko','markerfacecolor','r','markersize',4);
    plot([0 2.4],[bed(end) bed(end)],'--','color',[0.55 0.27 0.07],'linewidth',1.5)
    text(1.9,bed(end),'?','color',[0.55 0.27 0.07],'fontsize',14);
    plot([3.4 5],[bed(end) bed(end)],'--','color',[0.55 0.27 0.07],'linewidth',1.5)
    %plot((flow(front1:front2)-flowline(344,1))/10^3,bed(front1:front2),'k--','linewidth',1.5)
    plot([flow(front1)-flowline(344,1) flow(front1)-flowline(344,1)]/10^3,[bed(front1) surf(front1)],'k--','linewidth',1.5);
    plot([flow(front2)-flowline(344,1) flow(front2)-flowline(344,1)]/10^3,[bed(front2) surf(front2)],'k--','linewidth',1.5);
    xlim([-32 5])
    set(gca,'ytick',[-1000:500:1000]);
    ylim([-1200 1200])
    ylabel('Elevation (masl)','fontsize',9,'fontname','arial');
    xlabel('Distance along flowline (km)','fontsize',9,'fontname','arial');
    set(gca,'fontsize',9);
    text(-30.5,800,'b','fontsize',9,'fontname','arial','fontweight','bold');
    legend(h,'Eulerian','Lagrangian','location','northeast');
    legendshrink(0.4);
    legend('boxoff');
    plot([-8.5 4.3 4.3 -8.5 -8.5],[400 400 1100 1100 400],'k');
    set(gca,'ticklength',[0.025 0.025])
    p=get(s,'position');
    set(s,'position',[p(1)+0.05 p(2)+0.1 p(3) p(4)])
    
end

%Terminus position and velocity change
%Case 1: 2010
%cd ~/Data/Velocity/TSX/Helheim/Outputs/
%a2010=load('vel_track-16848.mat');
%r2010=load('vel_track-17850.mat');

%Case 2: 2013
%r2013_1=load('vel_track-30542.mat');
%a2013_1=load('vel_track-31544.mat');
%r2013_2=load('vel_track-32713.mat');