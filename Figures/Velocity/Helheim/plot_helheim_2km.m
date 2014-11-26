data = load('~/Data/Velocity/Random/Helheim.3.mat');

vel=[data.vw' data.vtsx']';
time=[data.years data.dates']';

ind=find(vel<1000);
vel(ind)=[];
time(ind)=[];
ind=find(isnan(vel));
vel(ind)=[];
time(ind)=[];
[time,ind]=sort(time);
vel=vel(ind);

figure;
plot(time,vel,'o--','color','k','linewidth',1.5,'markerfacecolor','k');
ylabel('Velocity (m/yr)','fontsize',12);
xlabel('Year');
