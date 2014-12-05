% Figures for my poster for McCarthy summer school.

if ~exist('fronts','var')
    helheim_fronts
end
%if ~exist('velocity','var')
%    helheim_vel_flowline
%end

plot_timeline=0;
plot_crossection1=0;
plot_crossection2=1;

% Terminus and velocity timeline
if plot_timeline
    figure;
    set(gcf,'PaperUnits','centimeters');
    xSize = 40; ySize = 20;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[750 4 xSize*50 ySize*50])
    
    x = 3/xSize;
    y = 1/ySize; 
    width = 1-x;
    height = 1-y;
    axes('position',[x y width height])
    
    s(1)=subplot(2,1,1);
    plot(fronts(:,1),(fronts(:,2)-fronts(end,2)+10^3)/10^3,'--','linewidth',1.5,'color',[0.6 0.6 0.6]); hold on;
    plot(fronts(:,1),(fronts(:,2)-fronts(end,2)+10^3)/10^3,'k.','markersize',20)
    ylabel('Terminus (km)','fontsize',24,'fontname','arial')
    ylim([0 7.5]);
    xlim([2000.5 2014]);
    set(gca,'xtick',2001:2:2013)
    set(gca,'xticklabel',{})
    set(gca,'fontsize',20,'fontname','arial')
    set(gca,'ytick',0:2:6)
    
    s(2)=subplot(2,1,2);
    fronts_interp2=interp1(fronts(:,1),fronts(:,2),time);
    for i=1:length(time)
        [dist(i),lag_5km(i)]=min(abs(flowline(:,1)-(fronts_interp2(i)-5*10^3)));
        [~,eul_5km(i)]=min(abs(flowline(:,1)-(flowline(end,1)-5*10^3)));
    end
    nonnan = find(~isnan(velocity(eul_5km(1),:)));
    plot(time(nonnan),velocity(eul_5km(1),nonnan)/1000,'--','linewidth',1.5,'color',[0.6 0.6 0.6]); hold on;
    plot(time,velocity(eul_5km(1),:)/1000,'k.','markersize',20)
    ylabel('Velocity (km/yr)','fontsize',24,'fontname','arial')
    ylim([5 8.5])
    xlim([2000.5 2014]);
    set(gca,'fontsize',20,'fontname','arial')
    set(gca,'ytick',5:1:8)
    set(gca,'xtick',2001:2:2013)
    xlhand = get(gca,'xlabel')
    
    p=get(s(2),'position');
    set(s(2),'position',[p(1) p(2)+0.1 p(3) p(4)])
end

if plot_crossection1
    figure;
    set(gcf,'PaperUnits','centimeters');
    xSize = 25; ySize = 15;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[750 4 xSize*50 ySize*50])
    
    x = 3/xSize;
    y = 1/ySize; 
    width = 1-x;
    height = 1-y;
    axes('position',[x y width height])
    
    bed_fronts=interp1(flowline(:,1),flowline(:,4),fronts(:,2));
    patch([flowline(end,1)/10^3; flowline(end,1)/10^3+10; flowline(end,1)/10^3+10; flowline(end,1)/10^3],[0; 0; flowline(end,4); flowline(end,4)],[173/255 216/255 230/255],'edgecolor','none'); hold on;
    patch([flowline(end,1)/10^3; flowline(end,1)/10^3+10; flowline(end,1)/10^3+10; flowline(end,1)/10^3],[-1000; -1000; flowline(end,4); flowline(end,4)],[222/255 184/255 135/255],'edgecolor','none')
    patch([flowline(:,1)/10^3; flowline(end,1)/10^3; flowline(1,1)/10^3; flowline(1,1)/10^3],[flowline(:,4); -1000; -1000; flowline(1,4)],[222/255 184/255 135/255],'edgecolor','none')
    plot(flowline(:,1)/10^3,flowline(:,4),'k','linewidth',1.5);
    plot(flowline(:,1)/10^3,flowline(:,5),'k','linewidth',1.5);
    plot([flowline(end,1)/10^3 flowline(end,1)/10^3],[flowline(end,4) flowline(end,5)],'k','linewidth',1.5)
    plot([flowline(end,1)/10^3 flowline(end,1)/10^3+10],[flowline(end,4) flowline(end,4)],'k--','linewidth',1.5)
    xlim([0 50])
    xlabel('Distance along flowline (km)','fontsize',26,'fontname','arial')
    ylabel('Elevation (m asl)','fontsize',26,'fontname','arial')
    set(gca,'fontsize',20,'fontname','arial')
    
    plot([35 45 45 35 35],[-800 -800 200 200 -800],'k--','linewidth',1.5)
    
    p=get(gca,'position');
    set(gca,'position',[p(1)+0.05 p(2)+0.1 p(3)-0.1 p(4)-0.2])
end

if plot_crossection2
    figure;
    set(gcf,'PaperUnits','centimeters');
    xSize = 15; ySize = 25;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[750 4 xSize*50 ySize*50])
    
    x = 3/xSize;
    y = 1/ySize; 
    width = 1-x;
    height = 1-y;
    axes('position',[x y width height])
    
    bed_fronts=interp1(flowline(:,1),flowline(:,4),fronts(:,2));
    bed_surfs=interp1(flowline(:,1),flowline(:,5),fronts(:,2));
    patch([flowline(end,1)/10^3; flowline(end,1)/10^3+10; flowline(end,1)/10^3+10; flowline(end,1)/10^3],[0; 0; flowline(end,4); flowline(end,4)],[173/255 216/255 230/255],'edgecolor','none'); hold on;
    patch([flowline(end,1)/10^3; flowline(end,1)/10^3+10; flowline(end,1)/10^3+10; flowline(end,1)/10^3],[-1000; -1000; flowline(end,4); flowline(end,4)],[222/255 184/255 135/255],'edgecolor','none')
    patch([flowline(:,1)/10^3; flowline(end,1)/10^3; flowline(1,1)/10^3; flowline(1,1)/10^3],[flowline(:,4); -1000; -1000; flowline(1,4)],[222/255 184/255 135/255],'edgecolor','none')
    
    for i = 1:length(fronts)
        plot([fronts(i,2)/10^3 fronts(i,2)/10^3],[bed_fronts(i),bed_surfs(i)],'linewidth',1.5,'color','r')
    end
    plot(flowline(:,1)/10^3,flowline(:,4),'k','linewidth',1.5);
    plot(flowline(:,1)/10^3,flowline(:,5),'k','linewidth',1.5);
    plot([flowline(end,1)/10^3 flowline(end,1)/10^3],[flowline(end,4) flowline(end,5)],'k','linewidth',1.5)
    plot([flowline(end,1)/10^3 flowline(end,1)/10^3+10],[flowline(end,4) flowline(end,4)],'k--','linewidth',1.5)
    ylim([-800 300])
    set(gca,'ytick',[-500:500:500])
    xlim([35 45])
    xlabel('Distance along flowline (km)','fontsize',24,'fontname','arial')
    ylabel('Elevation (m asl)','fontsize',24,'fontname','arial')
    set(gca,'fontsize',20,'fontname','arial')
    
    p=get(gca,'position');
    set(gca,'position',[p(1)+0.05 p(2)+0.1 p(3)-0.1 p(4)-0.2]) 
end