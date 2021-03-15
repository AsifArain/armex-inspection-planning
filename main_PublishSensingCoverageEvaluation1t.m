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



exp_nums = 4:10
% covrgs   = 0.95:0.01:1.00
covrgs   = 1.0:-0.01:0.95%0.90


conf_matrix = zeros(numel(exp_nums),numel(covrgs));
dist_matrix = zeros(numel(exp_nums),numel(covrgs));

f1___matrix = zeros(numel(exp_nums),numel(covrgs));
fp___matrix = zeros(numel(exp_nums),numel(covrgs));
tp___matrix = zeros(numel(exp_nums),numel(covrgs));
p____matrix = zeros(numel(exp_nums),numel(covrgs));
r____matrix = zeros(numel(exp_nums),numel(covrgs));


for i = 1:numel(exp_nums)
    for j = 1:numel(covrgs)
        
        dir_ = sprintf('results/sensing_coverage_comparison/%s-%02d/%s/coverage%03d/',...
            'prismaforum5',exp_nums(i),'1t-armex',covrgs(j)*100);
        
        confs_file = load([dir_,'executed_confs_final.txt']);
        conf_this = size(confs_file,1);
        
        dist_this = load([dir_,'planned_dist_global.dat']);
        
        eval_this = load([dir_,'evaluation/evaluation.mat']);
        
        eval_this
        
        conf_matrix(i,j) = conf_this;
        dist_matrix(i,j) = dist_this;
        
        f1___matrix(i,j) = eval_this.f_measure;
        tp___matrix(i,j) = eval_this.num_of_true_positives;
        fp___matrix(i,j) = eval_this.num_of_false_positives;
        p____matrix(i,j) = eval_this.precision;
        r____matrix(i,j) = eval_this.recall;
        
        
        
    end
end




conf_nums = mean(conf_matrix);
conf_percn = (mean(conf_matrix)./max(mean(conf_matrix)))*100;


dist_nums = mean(dist_matrix);
dist_percn = (mean(dist_matrix)./max(mean(dist_matrix)))*100;


%figure; plot(mean(conf_matrix))
%figure; plot((mean(conf_matrix)./max(mean(conf_matrix)))*100)
%figure; plot(mean(dist_matrix))



% color_str = [205/256, 092/256, 092/256;...
%              046/256, 139/256, 087/256;...
%              075/256, 000/256, 130/256;...
%              139/256, 069/256, 019/256];
% color_str = lines(10);
% color_str = [     0    0.4470    0.7410;...
%              0.8500    0.3250    0.0980;...
%              205/256, 092/256, 092/256;...
%              046/256, 139/256, 087/256;...
%              075/256, 000/256, 130/256;...
%              139/256, 069/256, 019/256];
color_str = lines(3);
col1 = color_str(1,:);
col2 = color_str(2,:);
col3 = color_str(3,:);



%***************************************
%
%       PLOT - SENSING CONFIGURATIONS
%
%***************************************
h = figure('name','Sensing conf'); 
hold on;
h_conf(1) = plot(conf_nums,'-o','color',col1,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col1);
hXLabel = xlabel('Minimum sensing coverage (%)');
hYLabel = ylabel('Sensing configurations');
set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'     , ...
      'YMinorTick'  , 'off'     , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'XTick'       , (1:7),    ...  
      'YTick'       , (10:2:35),   ...  
      'LineWidth'   , 1         );
%xticklabels({'100%','99%','98%','97%','96%','95%'})
xticklabels({'100','99','98','97','96','95','94','93','92','91','90'})
%ylim([55,105])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set([gca,hXLabel,hYLabel],'FontSize',8);
%xlim([0.5,6.5])
%ylim([00,50])
print_file_name = 'fig/sensing-coverage-vs-confs';
print(print_file_name, '-painters', '-dpdf', '-r1200');
matlab2tikz([print_file_name,'.tex'],...
    'height','\figureheight',...
    'width','\figurewidth',...
    'extraAxisOptions','xlabel near ticks',...
    'extraAxisOptions','ylabel near ticks');
  
%***************************************
%
%       PLOT - TRAVELLING DIST
%
%***************************************
h = figure('name','Travelling distance'); 
hold on;
h_conf(1) = plot(dist_nums,'-o','color',col1,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col1);
hXLabel = xlabel('Minimum sensing coverage (%)');
hYLabel = ylabel('Travelling distance (m)');
set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'     , ...
      'YMinorTick'  , 'off'     , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'XTick'       , (1:7),    ...  
      'YTick'       , (300:5:350),   ...  
      'LineWidth'   , 1         );
%xticklabels({'100%','99%','98%','97%','96%','95%'})
xticklabels({'100','99','98','97','96','95','94','93','92','91','90'})
%ylim([55,105])
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set([gca,hXLabel,hYLabel],'FontSize',8);
%xlim([0.5,6.5])
%ylim([00,50])
print_file_name = 'fig/sensing-coverage-vs-dist';
print(print_file_name, '-painters', '-dpdf', '-r1200');
matlab2tikz([print_file_name,'.tex'],...
    'height','\figureheight',...
    'width','\figurewidth',...
    'extraAxisOptions','xlabel near ticks',...
    'extraAxisOptions','ylabel near ticks');

%***************************************
%
%       PLOT - OPERATIONAL COST
%
%***************************************
h = figure('name','Operational Cost - Percentage'); 
hold on;
h_conf(1) = plot(conf_percn,'-o','color',col1,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col1);
h_conf(2) = plot(dist_percn,'-o','color',col2,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col2);
hLegend = legend(h_conf,...
                 'Sensing configurations',...
                 'Travelling distance',...
                 'location','NorthOutside',...
                 'Orientation','horizontal'); %'North' %vertical
legend('boxoff');    
set(hLegend,'Interpreter','latex');
hXLabel = xlabel('Minimum sensing coverage (\Omega)');
hYLabel = ylabel('Operational cost (%)');
set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'     , ...
      'YMinorTick'  , 'off'     , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'XTick'       , (1:11),    ...  
      'YTick'       , (50:05:100),   ...  
      'LineWidth'   , 1         );
%xticklabels({'100%','99%','98%','97%','96%','95%'})
%xticklabels({'100','99','98','97','96','95'})
xticklabels({'1.00','0.99','0.98','0.97','0.96','0.95','0.94','0.93','0.92','0.91','0.90'})
ylim([50,50+((100-50)*1.05)]) %ylim([50,100+(5*.25)]) % ylim([50,100])
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set([gca,hXLabel,hYLabel],'FontSize',8);
%
print_file_name = 'fig/sensing-coverage-vs-operational-cost';
print(print_file_name, '-painters', '-dpdf', '-r1200');
matlab2tikz([print_file_name,'.tex'],...
    'height','\figureheight',...
    'width','\figurewidth',...
    'extraAxisOptions','xlabel near ticks',...
    'extraAxisOptions','ylabel near ticks');

%***************************************
%
%       PLOT - RECONSTRUCTION QUALITY
%
%***************************************
h = figure('name','Reconstruction Quality'); 
hold on;
h_conf(1) = plot(mean(p____matrix),'-o','color',col1,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col1);
h_conf(2) = plot(mean(r____matrix),'-o','color',col2,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col2);
h_conf(3) = plot(mean(f1___matrix),'-o','color',col3,'LineWidth',1,'MarkerSize',4,'MarkerFaceColor',col3);
hLegend = legend(h_conf,...
                 'Precision',...
                 'Recall',...
                 'F1-score',...
                 'location','NorthOutside',...
                 'Orientation','vertical'); %'North'
legend('boxoff');    
set(hLegend,'Interpreter','latex');
hXLabel = xlabel('Minimum sensing coverage (\Omega)');
hYLabel = ylabel('Gas source estimations');
set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'     , ...
      'YMinorTick'  , 'off'     , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'XTick'       , (1:11),    ...  
      'YTick'       , (0.8:0.05:1.00),   ...  
      'LineWidth'   , 1         );
%xticklabels({'100%','99%','98%','97%','96%','95%'})
%xticklabels({'100','99','98','97','96','95'})
xticklabels({'1.00','0.99','0.98','0.97','0.96','0.95','0.94','0.93','0.92','0.91','0.90'})
ylim([0.8,1.0])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set([gca,hXLabel,hYLabel],'FontSize',8);
%xlim([0.5,6.5])
%ylim([00,50])
print_file_name = 'fig/sensing-coverage-vs-source-estimation';
print(print_file_name, '-painters', '-dpdf', '-r1200');
matlab2tikz([print_file_name,'.tex'],...
    'height','\figureheight',...
    'width','\figurewidth',...
    'extraAxisOptions','xlabel near ticks',...
    'extraAxisOptions','ylabel near ticks');