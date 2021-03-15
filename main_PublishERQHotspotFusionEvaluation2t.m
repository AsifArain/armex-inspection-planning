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



exp_nums = 6:10;

w_nearness  = zeros(numel(exp_nums),1);
w_jsd       = zeros(numel(exp_nums),1);
w_tp        = zeros(numel(exp_nums),1);
w_fp        = zeros(numel(exp_nums),1);
w_precision = zeros(numel(exp_nums),1);
w_recall    = zeros(numel(exp_nums),1);
w_fscore    = zeros(numel(exp_nums),1);

x_nearness  = zeros(numel(exp_nums),1);
x_jsd       = zeros(numel(exp_nums),1);
x_tp        = zeros(numel(exp_nums),1);
x_fp        = zeros(numel(exp_nums),1);
x_precision = zeros(numel(exp_nums),1);
x_recall    = zeros(numel(exp_nums),1);
x_fscore    = zeros(numel(exp_nums),1);

x_confs     = zeros(numel(exp_nums),1);
w_confs     = zeros(numel(exp_nums),1);

x_dist      = zeros(numel(exp_nums),1);
w_dist      = zeros(numel(exp_nums),1);

for i = 1:numel(exp_nums)
    
    dir_x = sprintf('results/%s-%02d/%s/gdm-withouthotspotfusion/evaluation/',...
        'prismaforum5',exp_nums(i),'2t-armex');
    
    dir_w = sprintf('results/%s-%02d/%s/gdm-withhotspotfusion/evaluation/',...
        'prismaforum5',exp_nums(i),'2t-armex');
        
    x2_this = load([dir_x,'evaluation2.mat']);
    x3_this = load([dir_x,'evaluation3.mat']);
    
    w2_this = load([dir_w,'evaluation2.mat']);
    w3_this = load([dir_w,'evaluation3.mat']);
        
    x_nearness(i)   = x2_this.nearness;
    x_jsd(i)        = x2_this.jsd;
    x_tp(i)         = x3_this.num_of_true_positives;
    x_fp(i)         = x3_this.num_of_false_positives;
    x_precision(i)  = x3_this.precision;
    x_recall(i)     = x3_this.recall;
    x_fscore(i)     = x3_this.f_measure;
        
    w_nearness(i)   = w2_this.nearness;
    w_jsd(i)        = w2_this.jsd;
    w_tp(i)         = w3_this.num_of_true_positives;
    w_fp(i)         = w3_this.num_of_false_positives;
    w_precision(i)  = w3_this.precision;
    w_recall(i)     = w3_this.recall;
    w_fscore(i)     = w3_this.f_measure;
        
end


for i = 1:numel(exp_nums)
            
    dir_x = sprintf('results/%s-%02d/%s/gdm-withouthotspotfusion/',...
        'prismaforum5',exp_nums(i),'2t-armex');
    
    dir_w = sprintf('results/%s-%02d/%s/gdm-withhotspotfusion/',...
        'prismaforum5',exp_nums(i),'2t-armex');
    
    x_conf_this = load([dir_x,'executed_confs_final.txt']);
    w_conf_this = load([dir_w,'executed_confs_final.txt']);
    
    x_dist_this = load([dir_x,'planned_dist_global.dat']);
    w_dist_this = load([dir_w,'planned_dist_global.dat']);
    
    x_confs(i)   = size(x_conf_this,1);
    w_confs(i)   = size(w_conf_this,1);
    
    x_dist(i)    = x_dist_this;
    w_dist(i)    = w_dist_this;
        
end


