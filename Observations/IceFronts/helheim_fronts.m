%This file collects the Helheim ice front and rift positions from
%~2010-2013. I have yet to digitize rift positions for 2010-2011, so all
%comparison of rift position and terminus position should be limited to
%2012-present.
%
%LMK, UW, 5/1/2014

cd ~/Code/IceFronts/
system('python helheim_compile.py');

junk = importdata('~/Code/Helheim/Modeling/SolverFiles/Flowline/Helheim/Inputs/flowline.dat',' ',1);
flowline = junk.data; clear junk

%Ice front positions
cd ~/Data/Shape_files/Ice_fronts/Helheim/
files=dir('*.dat');

m=(flowline(end-2,3)-flowline(end-1,3))/(flowline(end-2,2)-flowline(end-1,2));
b=flowline(end-2,3)-m*flowline(end,2);
values=500:500:10000;
for i = 1:length(values)
    extrap(i,2)=flowline(end,2)+values(i);
    extrap(i,3)=extrap(i,2)*m+b;
    extrap(i,1)=flowline(end,1)+sqrt((flowline(end,2)-extrap(i,2))^2+(flowline(end,3)-extrap(i,3))^2);
end
flowcoords=[flowline(:,1:3);extrap];
for i = 1:length(files)
    front(i).date=str2num(files(i).name(1:4))+str2num(files(i).name(6:8))/365.25;
    front(i).coords=load(files(i).name);
    
    %Find intersection of ice front with the flowline and then compute the
    %distance along the flowline to that point
    intersect=InterX(flowcoords(:,2:3)',front(i).coords');
    [~,ind]=min(abs(intersect(1)-flowcoords(:,2)));
    if flowcoords(ind,2) > intersect(1) %make sure that we have the smaller value
        ind=ind-1;
    end
    dist=flowcoords(ind,1)+sqrt((intersect(1)-flowcoords(ind,2))^2+(intersect(2)-flowcoords(ind,3))^2);
    fronts(i,1)=front(i).date;
    fronts(i,2)=dist;
end

%Rift positions
cd ~/Data/Shape_files/Rifts/Helheim/
files=dir('*.dat');

for i = 1:length(files)
    rift(i).date=str2num(files(i).name(1:4))+str2num(files(i).name(6:8))/365.25;
    rift(i).coords=load(files(i).name);
    
    %Find intersection of rifts with the flowline and then compute the
    %distance along the flowline to that point
    intersect=InterX(flowcoords(:,2:3)',rift(i).coords');
    [~,ind]=min(abs(intersect(1)-flowcoords(:,2)));
    if flowcoords(ind,2) > intersect(1) %make sure that we have the smaller value
        ind=ind-1;
    end
    dist=flowcoords(ind,1)+sqrt((intersect(1)-flowcoords(ind,2))^2+(intersect(2)-flowcoords(ind,3))^2);    
    
    rifts(i,1)=rift(i).date;
    rifts(i,2)=dist;
end