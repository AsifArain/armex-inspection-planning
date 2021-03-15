clear all; close all; clc;

disp('#____________________________________________________________________#')
disp('#                                                                    #')
disp('#                                                                    #')
disp('#                            SPP-REM                                 #')
disp('#       ------------------------------------------------------       #')
disp('#           Sensor Planning for Robot Emission Monitoring            #')
disp('#       ------------------------------------------------------       #')
disp('#                      One-step Monitoring                           #')
disp('#                      Two-step Monitoring                           #')
disp('#                    Human-expert Monitoring                         #')
disp('#                                                                    #')
disp('#                     ---- EVALUATION ----                           #')
disp('#                                                                    #')
disp('#                     Author: Asif Arain                             #')
disp('#                                                                    #')
disp('#____________________________________________________________________#')



exp_nums = 4:10;


tsp_matrix = zeros(numel(exp_nums),1);
ees_matrix = zeros(numel(exp_nums),1);

for i = 1:numel(exp_nums)
        
        
    dir_tsp = sprintf('results/%s-%02d/%s/tsp/',...
        'prismaforum5',exp_nums(i),'1t-armex');

    dir_ees = sprintf('results/%s-%02d/%s/ees/',...
        'prismaforum5',exp_nums(i),'1t-armex');


    tsp_this = load([dir_tsp,'planned_dist_global.dat']);
    ees_this = load([dir_ees,'planned_dist_global.dat']);

    tsp_matrix(i) = tsp_this;
    ees_matrix(i) = ees_this;
        
        
end



tsp_mean = mean(tsp_matrix)
ees_mean = mean(ees_matrix)

