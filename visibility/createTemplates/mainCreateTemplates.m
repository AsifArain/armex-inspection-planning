close all; clear all; clc;



fov = [45,90,180,300,360]%270%45%90%[270]; %300
range = [5,10]%20%[30];
fovfaces = [360];



for i = 1:numel(fov)
    for j = 1:numel(range)
        for k = 1%:numel(fovfaces)            
                
            FoV = fov(i)
            SensorRange.cell = range(j)
            FoVfaces.num = fovfaces(k)
            
            %profile on;
            
            [ tDomain ] = fDomain( FoV,FoVfaces,SensorRange )
            %[ tDomain ] = fCopyDomain( FoV,FoVfaces,SensorRange )
            
            
            
            %profile_name = sprintf('profiles/range%02d_fov%03d_fovfaces%03d',...
            %                SensorRange.cell,FoV,FoVfaces.num)
            
            %mkdir(profile_name)            
            %profile off;
            %profsave(profile('info'),profile_name);
            
            %pause
            
        end
    end
end


