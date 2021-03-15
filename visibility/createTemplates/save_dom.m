close all; clear all; clc;



fov = 90;
range = 15;
fovfaces = 360;

f = fovfaces;




FoV = fov
SensorRange.cell = range
FoVfaces.num = fovfaces


filename = sprintf('../sD_oD_VObs_Tomography_range%02d_fov%03d_fovfaces%03d.mat',...
                    SensorRange.cell,FoV,f);
load(filename);


sD_all = sD;
oD_all = oD;
VObs_all = VObs;


sD = [];
oD = [];
VObs = [];



for j = 1:f
    
    
    dir__=sprintf('../sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/orientation_%03ddeg/',SensorRange.cell,FoV,j);
    
    if ~exist(dir__,'dir')
        mkdir(dir__);
    end
    
    sD = sD_all{j};
    oD = oD_all{j};
    VObs = VObs_all{j};
    
    filename = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
                        SensorRange.cell,FoV,j);
    save([dir__,filename],'sD','oD','VObs','-V7.3');

end



