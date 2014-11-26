cd ~/Data/Velocity/TSX/Helheim/

D=dir('track*');

for i=1:length(D)
    clear x y vx vy ex ey v info date file_out
    file_out_mat=sprintf('Outputs/vel_%s.mat',D(i).name);
    file_out_txt=sprintf('Outputs/vel_%s.txt',D(i).name);
    if exist(file_out_txt,'file')
        fprintf('Already imported \n');
    else  
        cd(D(i).name)
        [x,y,vx,vy,ex,ey,v]= read_velocity('mosaicOffsets',0);
        
        info=fileread('mosaicOffsets.meta');
        date=jdate2doy(str2num(info(36:47)));
        
        cd ../
        save(file_out_mat,'x','y','vx','vy','ex','ey','v','info','date')
        
        n=1;
        for i = 1:length(x)
            for j = 1:length(y)
                x1(n)=x(i);
                y1(n)=y(j);
                vx1(n)=vx(j,i);
                ex1(n)=ex(j,i);
                vy1(n)=vy(j,i);
                ey1(n)=ey(j,i);
                v1(n)=v(j,i);
                n=n+1;
            end
        end
        dlmwrite(file_out_txt,[x1' y1' v1' vx1' ex1' vy1' ey1'],'delimiter',' ','precision','%.6f');
    end
end

%Combine all the velocities into one matrix
cd ~/Data/Velocity/TSX/Helheim/Outputs/
D=dir('*.mat');
for i=1:length(D)
  load(D(i).name)
  if i==1
      xmin=min(x);
      xmax=max(x);
      ymin=min(y);
      ymax=max(y);
  end
  if min(x) < xmin
      xmin=min(x);
  end
  if max(x) > xmax
      xmax=max(x);
  end
  if min(y) < ymin
      ymin=min(y);
  end   
    if max(y) > ymax
      ymax=max(y);
  end   
end
xall=xmin:100:xmax;
yall=ymin:100:ymax;

dates=zeros(length(D));
velocities=zeros(length(yall),length(xall),length(D));
velocities(:,:,:)=NaN;
for i=1:length(D)
    load(D(i).name)
    dates(i)=date;
    [~,xmin]=min(abs(x(1)-xall));
    [~,xmax]=min(abs(x(end)-xall));
    [~,ymin]=min(abs(y(1)-yall));
    [~,ymax]=min(abs(y(end)-yall));
    velocities(ymin:ymax,xmin:xmax,i)=v;
end

%Compute some statistics
meanvelocities=zeros(length(yall),length(xall));
meanvelocities(:,:)=NaN;
numbers=zeros(length(yall),length(xall));
numbers(:,:)=NaN;
stdvelocities=zeros(length(yall),length(xall));
stdvelocities(:,:)=NaN;
minvelocities=zeros(length(yall),length(xall));
minvelocities(:,:)=NaN;
maxvelocities=zeros(length(yall),length(xall));
maxvelocities(:,:)=NaN;
for i=1:length(xall)
    for j=1:length(yall)
        nonnan=find(~isnan(velocities(j,i,:)));
        numbers(j,i)=length(nonnan);
        if length(nonnan) > 0
            meanvelocities(j,i)=mean(velocities(j,i,nonnan));
            stdvelocities(j,i)=std(velocities(j,i,nonnan));
            minvelocities(j,i)=min(velocities(j,i,nonnan));
            maxvelocities(j,i)=max(velocities(j,i,nonnan));
        end
    end
end