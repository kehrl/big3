%Load data
data=importdata('~/Data/Climate/IceTemperature/Helheim/xyzTAhelheim.txt',',',1);
data=data.data;
data(:,4)=data(:,4);

%Load mesh extent
mesh=load('~/Code/Modeling/SolverFiles/3D/Helheim/Inputs/mesh_extent.dat');

%Find temperature measurements inside our mesh area
rectangle=[min(mesh(:,1))-5000 min(mesh(:,2))-5000; min(mesh(:,1))-1000 max(mesh(:,2))+5000; 
    max(mesh(:,1))+5000 max(mesh(:,2))+5000; max(mesh(:,1))+5000 min(mesh(:,2))-5000]; 
inside=inpolygon(data(:,1),data(:,2),rectangle(:,1),rectangle(:,2));
data=data(inside,:);

if length(unique(data(:,1))) > length(unique(data(:,2)))
    [~,ind_locs]=unique(data(:,1));
else
    [~,ind_locs]=unique(data(:,2));
end

for i=1:length(ind_locs)
    junk=find(data(:,1)==data(ind_locs(i),1));
    column=find(data(junk,2)==data(ind_locs(i),2));
    column=junk(column);
    [~,ind_zs]=max(data(column,3));
    Ts(i)=data(column(ind_zs),4);
    [~,ind_zb]=min(data(column,3));
    Tb(i)=data(column(ind_zb),4);
    Tave(i)=mean(data(column,4));
    [~,ind_min]=min(data(column,4));
    Zmin(i)=data(column(ind_min),3);
    Tmin(i)=data(column(ind_min),4);
end

fid = fopen('~/Code/Modeling/SolverFiles/3D/Helheim/Inputs/icetemperature.xyz','w');
fprintf(fid,'%d\n',length(data));
for i = 1:length(data)
    fprintf(fid,'%f %f %f %f\n',data(i,1),data(i,2),data(i,3),data(i,4));
end