function [ mean_map,...
           cell_coords_x,...
           cell_coords_y ] = fReconstructionVideo( conf_num,...
                                                   para_,...
                                                   dir_)


% fReconstruction performs reconstruction from the integral concentration measurements.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)



cellsize_recon = para_.ReconTomographyCellSizeM;
        
%====================================================
% -- measurement matrixL
%====================================================
fprintf('- measurement matrix for line model.\n');

M = [];
for i = 1:conf_num
    measu_file = sprintf('measurement_conf%02d.dat',i);
    M_current = dlmread([dir_.logs,measu_file]);
    M = [M;M_current];
end

%-- final M matrix (new)
M = M(:,[3,4,5,6,8]);



%====================================================
%	LINE MODEL 
%====================================================
fprintf('---> Generating line model.\n');

recon_size = 'environment-map';
[ A,...
  ppmm,...
  cell_coords_x,...
  cell_coords_y,...
  num_cells_x,...
  num_cells_y ] = fLineModelJFR( M,...
                                 recon_size,...
                                 cellsize_recon,...
                                 para_,...
                                 dir_ );
A = full(A);


% %--- save cell_coords_x
% x_coord_file = sprintf('x_coord.dat');
% fileID = fopen([dir_.Recon,x_coord_file],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,x_coord_file],cell_coords_x,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,x_coord_file],'a'); fclose(fileID);
% 
% %--- save cell_coords_y
% y_coord_file = sprintf('y_coord.dat');
% fileID = fopen([dir_.Recon,y_coord_file],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,y_coord_file],cell_coords_y,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,y_coord_file],'a'); fclose(fileID);


%==================================================== 
%----  GAS DISTRIBUTION MAP
%====================================================

fprintf('---> Building gas distribution map.\n');
visualize = 0;
[ mean_map ] = fGDM( A,...
                     ppmm,...
                     cell_coords_x,...
                     cell_coords_y,...
                     num_cells_x,...
                     num_cells_y,...
                     visualize );
% map_recon = mean_map';
% 
% %--- save mean map all
% recon_file = sprintf('reconstruction.dat');
% fileID = fopen([dir_.Recon,recon_file],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,recon_file],map_recon,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,recon_file],'a'); fclose(fileID);
% 
% map_recon_col = mean_map;
% gdm_colormap = flipud(autumn(512)); % for main fig.
% max(map_recon_col(:))
% max_concentration = min( [max(map_recon_col(:)),para_.PublishUpperPPM] ) % main fig.
% delta = max_concentration/(size(gdm_colormap,1)-1);
% % linear_color_number = round(map_recon_col(:)./delta)+1;
% linear_color_number = round(min(map_recon_col(:),max_concentration)./delta)+1;
%  
% 
% colR = round(gdm_colormap(linear_color_number,1)*255);
% colG = round(gdm_colormap(linear_color_number,2)*255);
% colB = round(gdm_colormap(linear_color_number,3)*255);
% 
% thisfile = sprintf('reconstructionColorR.dat');
% fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,thisfile],colR,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);
% 
% thisfile = sprintf('reconstructionColorG.dat');
% fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,thisfile],colG,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);
% 
% thisfile = sprintf('reconstructionColorB.dat');
% fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,thisfile],colB,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);






%---- VISUALIZATION
%====================================================
                
% if visualize == 1
%     
%     %map_env = map_env;
%     
%     figure('name','reconstrucion on environment map'); 
%     imshow(map_env','InitialMagnification',800); hold on;
%     set(gca,'YDir','normal');
%     colormap('gray'); caxis([-1 1.25]);
%            
%     map_recon = mean_map';
%     
%     
%     % ------------- select color map for gdm --------------
%     gdm_colormap = flipud(hot(512)); % for main fig. 
%            
%     % ------ color step size --------------------------------    
%     max_concentration = max(map_recon(:)); % main fig.
%     delta = max_concentration/(size(gdm_colormap,1)-1);
%     
%     
%     
%     for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
%         for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
%             % if positive concentration
%             
%             if map_recon(i,j)>=0
%                 
%                 if map_env(i,j)>0
%                 % ---------------------- useful trash -------------------
%                 %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);
%                 
%                 % --- if concentration is greater than 1000 ppm, its still 1000 ppm
%                 %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
%                 %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
%                 %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
%                 % -------------------------------------------------------
%                 
%                 % ----------------------------------------------------
%                 % linear color code
%                 % ----------------------------------------------------
%                 
%                 linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
%                 col = gdm_colormap( linear_color_number,:);                
%                 
%                 
%                 %----- to avoid black borders of earch cell -----
%                 % {
%                 if sum(col) == 3
%                     col = col-(1e-10);
%                 end
%                 %}
%                 % ----- plot ------------
%                 plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
%                     's',...
%                     'MarkerEdgeColor',col,...
%                     'MarkerFaceColor',col,...
%                     'MarkerSize',5);
%                 end
%             end
%         end
%     end
% 
% end
% 
% 
%                 
% if visualize == 1
%     %pause
%     % ----------------    
%     % -- splite M matrix
%     start_x__all = M(:,1);
%     start_y__all = M(:,2);
%     end___x__all = M(:,3);
%     end___y__all = M(:,4);
%     plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
%     %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')
%     axis equal
%     
%     %figure; imshow(mean_map','InitialMagnification',800); hold on;
%     %set(gca,'YDir','normal')
%     %colormap('hot'); % not go back to gray scale.
%     %%plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
%     %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')    
%     %pause(1)
%     
% end

end