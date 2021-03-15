function [ mse,...
           nearness,...
           ple,...
           md,...
           nae,...
           nk,...
           psnr,...
           sc,...
           mssim,...
           kld,...
           jsd,...
           d_eucl,...
           d_corr] = fEvaluation( recon_filename,...
                                  para_,...
                                  dir_,...
                                  visualize )


% fEvaluation performs evaluation of reconstructions.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)


% fprintf('\n\nEvaluation:\n');
% ______________________________________________________________________ %
%                   upload maps and reconstructions                      %
% ______________________________________________________________________ %


% mean_map = dlmread([results_dir,'mean_map.dat']);
% mean_map = dlmread([dir_.recon_mean,gdm_filename]);

% --- cell_coords_x
% cell_coords_x_file_name = strrep(gdm_filename,'mean_map','cell_coords_x');
% cell_coords_y_file_name = strrep(gdm_filename,'mean_map','cell_coords_y');

% cell_coords_x = dlmread([dir_.recon_cellx,cell_coords_x_file_name]);
% cell_coords_y = dlmread([dir_.recon_celly,cell_coords_y_file_name]);


load([dir_.recon,recon_filename],'mean_map','cell_coords_x','cell_coords_y')
mean_map = mean_map';

% map_conc_______gt = dlmread([dir_.maps,'map_conc_______gt.dat']);
% map_conc_______gt = dlmread([dir_.maps,gt_filename]);

% gt_filename = sprintf('map_conc_%d.dat',para_.GroundTruthTimeStamp);
gt_filename = sprintf('map_con.dat');
map_conc_______gt = dlmread([dir_.GroundTruthLogs,gt_filename]);


% ______________________________________________________________________ %
%                           ground truth                                 %
% ______________________________________________________________________ %

% grndtruth_m = map_conc_______gt(round(cell_coords_x)+1,round(cell_coords_y));
grndtruth_m = map_conc_______gt(round(cell_coords_x),round(cell_coords_y));


% figure; pcolor(grndtruth_m)
% figure; pcolor(mean_map)

if visualize == 1
    figure; 
    subplot(1,2,1); pcolor(grndtruth_m); axis equal
    subplot(1,2,2); pcolor(mean_map);    axis equal
end

if visualize == 1
    figure; 
    subplot(1,2,1);
    image(cell_coords_x,cell_coords_y,grndtruth_m,'CDataMapping','scaled');
    %image(1:100,1:100,gt,'CDataMapping','scaled');hold on;        
    % title(sprintf('Mean Map, file %s', A_FILE));
    colormap hot;
    %colormap gray; 
    colormap(flipud(colormap))
    colorbar;
    %caxis([0 5000]);
    %xlim([0,100]); ylim([0,100]);
    %}
    hold on;
    axis equal;

    subplot(1,2,2);
    image(cell_coords_x,cell_coords_y,mean_map,'CDataMapping','scaled');
    %image(1:100,1:100,gt,'CDataMapping','scaled');hold on;        
    % title(sprintf('Mean Map, file %s', A_FILE));
    colormap hot; 
    %colormap gray; 
    colormap(flipud(colormap))
    %title([num2str(table(i,2:n_c+1)),' <',...
    %    num2str(rad2deg(normalizeAngle(deg2rad(table(i,3)-table(i,2))))),'>']);
    colorbar;
    %caxis([0 5000]);
    %xlim([0,100]); ylim([0,100]);
    %}
    hold on;
    axis equal;
end


% size(grndtruth_m)
% size(mean_map)

% Display
% figure;
% subplot(2,3,1);
% pcolor(environment_model.cell_coords_x(1:end-1),environment_model.cell_coords_y(1:end-1),grndtruth_m); 
% hold on; 
% title('Ground Truth (Original Resolution)');
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
% 
% groundTruth.meanMap.originalResolution = grndtruth_m;


% ______________________________________________________________________ %
%                           MSE                                          %
% ______________________________________________________________________ %

mse = mean( (grndtruth_m(:)-mean_map(:)).^2 );
fprintf('- MSE = %f\n',mse);
% new_MSE = MeanSquareError(grndtruth_m, mean_map)

% ______________________________________________________________________ %
%                           Nearness                                     %
% ______________________________________________________________________ %

nom = sum( (grndtruth_m(:)-mean_map(:)).^2 );
den = sum( ( grndtruth_m(:)- (sum(grndtruth_m(:))/numel(grndtruth_m)) ).^2 );

nearness = nom/den;
fprintf('- Nearness = %f\n',nearness);

% ______________________________________________________________________ %
%                         peak location error                            %
% ______________________________________________________________________ %

[~,num_gt] = max(grndtruth_m(:));
[gt_row,gt_col] = ind2sub(size(grndtruth_m),num_gt);

[~,num_mm] = max(mean_map(:)); 
[mm_row,mm_col] = ind2sub(size(mean_map),num_mm);


ple = sqrt((mm_row-gt_row)^2 + (mm_col-gt_col)^2);
fprintf('- Peak location error = %f\n',ple);

% ______________________________________________________________________ %
%                         maximum difference                             %
% ______________________________________________________________________ %

md = MaximumDifference(grndtruth_m, mean_map);
fprintf('- Maximum difference = %f\n',md);


% ______________________________________________________________________ %
%                     normalized absolute error                          %
% ______________________________________________________________________ %

nae = NormalizedAbsoluteError(grndtruth_m,mean_map);
fprintf('- Normalized absolute error = %f\n',nae);


% ______________________________________________________________________ %
%                     normalized cross correlation                       %
% ______________________________________________________________________ %

nk = NormalizedCrossCorrelation(grndtruth_m,mean_map);
fprintf('- Normalized cross correlation = %f\n',nk);


% ______________________________________________________________________ %
%                     Peak signal to noice ratio                         %
% ______________________________________________________________________ %

psnr = PeakSignaltoNoiseRatio(grndtruth_m,mean_map);
fprintf('- Peak signal to noise ration = %f\n',psnr);


% ______________________________________________________________________ %
%                          Structural content                            %
% ______________________________________________________________________ %

sc = StructuralContent(grndtruth_m,mean_map);
fprintf('- Structural content = %f\n',sc);


% ______________________________________________________________________ %
%                   Structural SIMilarity (SSIM) index                   %
% ______________________________________________________________________ %

% [mssim,~] = ssim_index(grndtruth_m,mean_map);

K = [0.05 0.05];
% window = ones(8);
window = ones(4);
L = max(grndtruth_m(:));%100;
[mssim ssim_map] = ssim_index(grndtruth_m,mean_map,K,window,L);
fprintf('- Structural SIMilarity (SSIM) index = %f\n',mssim);

% figure; imshow(max(0, ssim_map).^4)


% ______________________________________________________________________ %
%                   Histogram (probability distributions)                %
% ______________________________________________________________________ %
edges = linspace(0,max(grndtruth_m(:)),5000);

[N_gt,edges] = histcounts(grndtruth_m,edges);
pdist_gt     = N_gt/numel(grndtruth_m);

[N_recon,edges] = histcounts(mean_map,edges);
pdist_recon     = N_recon/numel(mean_map);


% ______________________________________________________________________ %
%                        Kullback-Leibler Divergence                     %
% ______________________________________________________________________ %

[ kld ] = hcompare_KL( pdist_gt',pdist_recon');
fprintf('- Kullback-Leibler Divergence = %f\n',kld);

% ______________________________________________________________________ %
%                        Jensen-Shannon divergence                       %
% ______________________________________________________________________ %

jsd = JSDiv( pdist_gt,pdist_recon );
fprintf('- Jensen-Shannon divergence = %f\n',jsd);

% ______________________________________________________________________ %
%         Pairwise distance between two observations - euclidean         %
% ______________________________________________________________________ %

d_eucl = pdist2(pdist_gt,pdist_recon,'euclidean');
fprintf('- Pairwise euclidean dist = %f\n',d_eucl);


% ______________________________________________________________________ %
%         Pairwise distance between two observations - correlation       %
% ______________________________________________________________________ %

d_corr = pdist2(pdist_gt,pdist_recon,'correlation');
fprintf('- Pairwise correlation dist = %f\n',d_corr);





% pause
% close all

return
