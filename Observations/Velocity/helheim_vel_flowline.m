%This script loads and plots velocities along the Helheim flowline from
%2000-present (RADARSAT and TerraSAR-X).

%Load the flowline
junk = importdata('~/Dropbox/Code/Solver_files/Flowline/Helheim/Inputs/flowline.dat',' ',1);
flowline = junk.data; clear junk

clear time velocity xcoord ycoord

TSX=1;
RADARSAT=1;

%Make plots?
plot_bw=0;
plot_color=0;

%Load all TSX velocities
if TSX
cd ~/Data/Velocity/TSX/Helheim/Outputs/
D=dir('~/Data/Velocity/TSX/Helheim/Outputs/vel*.mat');
velocity=zeros(length(flowline),length(D));
N=length(D);
for i = 1:length(D)
    load(D(i).name)
    time(i) = date;
    
    %ind = min_dist_grid(x,y,coords);
    ind = min_dist_grid(x,y,flowline(:,2:3));
    for j = 1:length(ind)
        velocity(j,i)=v(ind(j,2),ind(j,1));
        xcoord(j,i)= x(ind(j,1));
        ycoord(j,i)=y(ind(j,2));
    end
end
else N=0;
end

if RADARSAT
%Load all RADARSAT velocities
cd ~/Data/Velocity/RADARSAT/Helheim/Outputs/
D=dir('~/Data/Velocity/RADARSAT/Helheim/Outputs/vel*.mat');
M=length(D);
velocity=[velocity zeros(length(flowline),length(D))];
for i = 1:length(D)
    load(D(i).name)
    time(i+N) = date;
    ind = min_dist_grid(x,y,flowline(:,2:3));
    %ind=min_dist_grid(x,y,coords);
    for j = 1:length(ind)
        velocity(j,i+N)=v(ind(j,2),ind(j,1));
        xcoord(j,i+N)= x(ind(j,1));
        ycoord(j,i+N)=y(ind(j,2));
    end
end
else M=0;
end

%Sort by date
[time,ind]=sort(time);
velocity=velocity(:,ind);
xcoord=xcoord(:,ind);
ycoord=ycoord(:,ind);

if plot_bw
    figure;
    fraction(1)=0;
    for i = 2:N+M
        fraction(i)=fraction(i-1)+(time(i)-time(i-1))/(time(end)-time(1));
    end
    colors=(fraction*0.8); clear fraction
    plot(flowline(:,3)/10^3,velocity(:,1),'color',[colors(1) colors(1) colors(1)],'linewidth',1); hold on;
    plot(flowline(:,3)/10^3,velocity(:,end),'color',[colors(end) colors(end) colors(end)],'linewidth',1);
    for i = 1:M+N
        plot(flowline(:,3)/10^3,velocity(:,i),'color',[colors(i) colors(i) colors(i)],'linewidth',1); 
    end
    legend(num2str(time(1)),num2str(time(end)));
    xlabel('Distance from terminus (km)','fontsize',12,'fontname','arial');
    ylabel('Velocity (m/yr)','fontsize',12,'fontname','arial');
end

if plot_color
    figure;
    ind=[1 2 3 4 5 10 14 19 49];
    colors=[1 0 0;
            1 0.5 0;
            1 1 0;
            0 1 0;
            0 1 1;
            0 0.5 1;
            0 0 1;
            0.5 0 0.5;
            0 0 0];
    
    for i = 1:length(colors)
        plot(flowline(:,1)/10^3,velocity(:,ind(i)),'color',colors(i,:),'linewidth',1.5); hold on;
    end
    for i = 1:length(colors)
        entries(i,1:6) = sprintf('%.1f',time(ind(i)));
    end
    legend(entries)
    %set(gca,'xdir','reverse')
    xlabel('x (km)','fontsize',12,'fontname','arial');
    ylabel('Velocity (m/yr)','fontsize',12,'fontname','arial');
end
