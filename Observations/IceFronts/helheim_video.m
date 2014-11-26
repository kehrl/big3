%Video of Helheim front positions from TerraSAR-X
%
%LMK, UW, 4/22/2014

dimensions=[3.03*10^5 3.15*10^5 -2.572*10^6 -2.583*10^6];

%For now just 2012 radar images but we can update this.
cd ~/Data/Mosaics/Helheim/
files=[dir('*.2009*1-20m*'); dir('*.2010*1-20m*'); dir('*.2011*1-20m*'); dir('*.2012*1-20m*'); dir('*.2013*1-20m*')];

%Get data
for i=1:length(files)
    time(i,1)=str2num(files(i).name(15:18));
    time(i,2)=str2num(files(i).name(20:22));
    time(i,3)=time(i,1)+time(i,2)/365.25;
    I(i) = geotiffread(files(i).name,dimensions);
    %Normalize data so that they are closer in scales
    %minvalue=min(min(I(i).z));
    %scale=230-double(minvalue);

    %I(i).z=(I(i).z-minvalue+15)*(230/scale);
end

limit=11/365.25; %if time difference between images exceeds this, we want to
%hold the frame so that timing is right in the video

n=1;
len=length(find(time(:,1) >= 2012));
cd ~/Desktop/test
clear writerObj
%writerObj = VideoWriter('helheim.avi','uncompressed avi');
%writerObj.FrameRate = 4;
%writerObj.Quality=100;
%open(writerObj);
for i=1:length(time)
    if time(i,1) >= 2012
    hold off;
    imshow(I(i).z,'xdata',I(i).x,'ydata',I(i).y); hold on;
    set(gca,'ydir','normal');
    text(dimensions(1)+500,dimensions(3)-800,datestr(doy2date(time(i,2),time(i,1))),'backgroundcolor','white','fontsize',24)
    plot([dimensions(2)-1700,dimensions(2)-700],[dimensions(4)+900 dimensions(4)+900],'k','linewidth',2)
    text(dimensions(2)-1700,dimensions(4)+500,'1 km','fontsize',24)
    %
    axis equal
    pos = get(gca,'position');
    set(gca,'position',[0 0.01 1 0.99])
    eval(sprintf('print -dpng test%.3d.png',n))
    n=n+1;
    %Take into account breaks between radar images
    ctime=time(i,3);
    if i ~=length(time(:,1))
    while time(i+1,3)-ctime > limit
        hold off;
        fig=imshow(I(i).z,'xdata',I(i).x,'ydata',I(i).y); hold on;
        set(gca,'ydir','normal');
        text(dimensions(1)+500,dimensions(3)-800,datestr(doy2date(time(i,2),time(i,1))),'backgroundcolor','white','fontsize',24)
        plot([dimensions(2)-1700,dimensions(2)-700],[dimensions(4)+900 dimensions(4)+900],'k','linewidth',2)
        text(dimensions(2)-1700,dimensions(4)+500,'1 km','fontsize',24)
        axis equal
        pos = get(gca,'position');
        set(gca,'position',[0 0.01 1 0.99])
        eval(sprintf('print -dpng test%.3d.png',n))
        eval(sprintf('print -dpng test%.3d.png',n))
        n=n+1;
        ctime=ctime+limit;
    end
    end
    %writeVideo(writerObj,frame);
    end
end
%close(writerObj);

%Save movie
cd ~/Desktop
%writerObj = VideoWriter('helheim.avi')
%movie2avi(Vid,'helheim.avi','fps',3,'quality',100)