cd ~/Data/Velocity/RADARSAT/Helheim/

D=dir('wint*');

for i=1:length(D)
    clear x y vx vy ex ey v info date file_out
    file_out=sprintf('Outputs/vel_%s.mat',D(i).name);
    if exist(file_out,'file')
        fprintf('Already imported \n');
    else  
        cd(D(i).name)
        [x,y,vx,vy,ex,ey,v]= read_velocity('mosaicOffsets',0);
        
        info=fileread('mosaicOffsets.meta');
        date=str2num(info(16:21));
        
        cd ../
        save(file_out,'x','y','vx','vy','ex','ey','v','info','date')
    end
end