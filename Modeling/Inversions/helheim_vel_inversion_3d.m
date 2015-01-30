%This code takes the global velocity map and fills in gaps, so that we can
%use it in inversions.

addpath('~/Code/Util/Mfiles')

%Large velocity map
cd ~/Data/Velocity/Random/Greenland/track-07to10/
[x,y,vx,vy,ex,ey,v]=read_velocity('mosaicOffsets',0);

%Individual velocity map
cd ~/Data/Velocity/TSX/Helheim/Outputs/
map = load('vel_track-31544');


%Chop map to desired size
[junk,xmin]=min(abs(x-200000));
[junk,xmax]=min(abs(x-321000));
[junk,ymin]=min(abs(y--2620000));
[junk,ymax]=min(abs(y--2510000));

x=x(xmin:xmax);
y=y(ymin:ymax);
v=v(ymin:ymax,xmin:xmax);
vx=vx(ymin:ymax,xmin:xmax);
vy=vy(ymin:ymax,xmin:xmax);
ex=ex(ymin:ymax,xmin:xmax);
ey=ey(ymin:ymax,xmin:xmax);

vx_fill=vx;
vy_fill=vy;
for i=1:length(y)
    nans=find(isnan(vx(i,:)));
    nonnans=find(~isnan(vx(i,:)));
    for j=1:length(nans)
       [~,ind]=min(abs(nans(j)-nonnans));
       vx_fill(i,nans(j))=vx(i,nonnans(ind));
       vy_fill(i,nans(j))=vy(i,nonnans(ind));
    end
end

%Final grid
xr=x(1):100:x(end);
yr=y(1):100:y(end);

Fx = griddedInterpolant(vx_fill,'linear');
vxr1=Fx({1:1/5:size(vx_fill,1),1:1/5:size(vx_fill,2)});
Fy = griddedInterpolant(vy_fill,'linear');
vyr1=Fy({1:1/5:size(vy_fill,1),1:1/5:size(vy_fill,2)});

[junk,xind1]=min(abs(xr-map.x(1)));
[~,xind2]=min(abs(xr-map.x(441)));
[~,yind1]=min(abs(yr-map.y(11)));
[~,yind2]=min(abs(yr-map.y(end)));

vxr2=vxr1;
vyr2=vyr1;
vmag2=zeros(size(vyr1));
vxr2(yind1:yind2,xind1:xind2)=map.vx(11:end,1:441);
vyr2(yind1:yind2,xind1:xind2)=map.vy(11:end,1:441);

for i=1:length(xr)
    for j=1:length(yr)
        if isnan(vxr2(j,i))
            vxr2(j,i)=vxr1(j,i);
            vyr2(j,i)=vyr1(j,i);
        end
    end
end

for i=1:length(xr)
    for j=1:length(yr)
        vmag2(j,i)=sqrt(vxr2(j,i)^2+vyr2(j,i)^2);
    end
end

%Convert to x,y,v form for saving
vmag_vec=reshape(vmag2,[length(xr)*length(yr),1]);
vx_vec=reshape(vxr2,[length(xr)*length(yr),1]);
vy_vec=reshape(vyr2,[length(xr)*length(yr),1]);
[xgrid,ygrid]=meshgrid(xr,yr);
x_vec=reshape(xgrid,[length(xr)*length(yr),1]);
y_vec=reshape(ygrid,[length(xr)*length(yr),1]);

cd ~/Code/Helheim/Modeling/SolverFiles/3D/Inputs/
fid = fopen('UDEM.xy','w');
fprintf(fid,'%d\n%d\n',length(xr),length(yr));
fclose(fid);
dlmwrite('UDEM.xy',[x_vec y_vec vx_vec],'-append','delimiter',' ','precision','%.1f')

fid = fopen('VDEM.xy','w');
fprintf(fid,'%d\n%d\n',length(xr),length(yr));
fclose(fid);
dlmwrite('VDEM.xy',[x_vec y_vec vy_vec],'-append','delimiter',' ','precision','%.1f')

fid = fopen('VMAG.xy','w');
fprintf(fid,'%d\n%d\n',length(xr),length(yr));
fclose(fid);
dlmwrite('VMAG.xy',[x_vec y_vec vmag_vec],'-append','delimiter',' ','precision','%.1f')


