function fPublishHumanReconstruction( conf_sequence_num,para_ )
%=====================================================================
%   AN ASSISTANCE FOR HUMAN EXPERT TO EXPLORE THE ENVIRONMENT
%=====================================================================
% 


para_.SamplingType = 'adaptive';

%=============================
%---- directories
%=============================
[ dir_ ] = fDirHumanJFR( para_ );
dir_.Recon = dir_.Solutions;
dir_.logs = dir_.MeasurementLogs;

%-- plan type
para_.PlanType = 'human-exploration';
dir_.ROSLogsThis = sprintf('%s%s/',dir_.ROSLogs,para_.PlanType);


%-- perform reconstruction
%--------------------------------    
visualize = 1;
[ M,...
  mean_map,...
  cell_coords_x,...
  cell_coords_y ] = fReconstructionAdaptive( conf_sequence_num,...
                                             para_,...
                                             dir_,...
                                             visualize );
pause
%-- publish this reconstruction
%-----------------------------------
fPublishReconstruction( mean_map,...
                        cell_coords_x,...
                        cell_coords_y,...
                        M,...
                        para_,...
                        dir_ );

end


function fPublishReconstruction( mean_map,...
                                 cell_coords_x,...
                                 cell_coords_y,...
                                 M,...
                                 para_,...
                                 dir_ )
% 
% 



% --- environment map data
% ================================================

file__map_env       = sprintf('%s_map_%s.dat',para_.environment,para_.EnvMapToPublish);
file__origin_env    = sprintf('%s_origin_%s.dat',para_.environment,para_.EnvMapToPublish);
file__cellsize_env  = sprintf('%s_cellsize_%s.dat',para_.environment,para_.EnvMapToPublish);

map_env      = load([dir_.env,file__map_env]);
origin_env   = load([dir_.env,file__origin_env]);
cellsize_env = load([dir_.env,file__cellsize_env]);


% --- plot environment
% ================================================

figure('name','reconstrucion on environment map'); 
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal');
colormap('gray'); caxis([-1 0.15]);

plot(origin_env(1),origin_env(2),'*r');

% --- plot reconstruction
% ================================================

map_recon = mean_map';

% ------------- select color map for gdm --------------
gdm_colormap = flipud(hot(512)); % for main fig. 

% ------ color step size --------------------------------    
max_concentration = max(map_recon(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);



for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
    for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
        % if positive concentration

        if map_recon(i,j) >0 %>=0

            %if map_env(i,j)>0
            % ---------------------- useful trash -------------------
            %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);

            % --- if concentration is greater than 1000 ppm, its still 1000 ppm
            %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
            %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
            %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
            % -------------------------------------------------------

            % ----------------------------------------------------
            % linear color code
            % ----------------------------------------------------

            linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number,:);                


            %----- to avoid black borders of earch cell -----
            % {
            if sum(col) == 3
                col = col-(1e-10);
            end
            %}
            % ----- plot ------------
            
            
            %plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
            plot( (cell_coords_x(i)*(para_.ReconTomographyCellSizeM/cellsize_env) ) +0.0,...
                  (cell_coords_y(j)*(para_.ReconTomographyCellSizeM/cellsize_env) ) +0.0,...
                 's',...
                 'MarkerEdgeColor','k',... 'MarkerEdgeColor',col,...
                 'MarkerFaceColor',col,...
                 'MarkerSize',5); %5
            %end
        end
    end
end



%--- sensing positions
%----------------------------------
sensingPositions = [unique(M(:,1)),unique(M(:,2))];
plot(sensingPositions(:,1),sensingPositions(:,2),'og');
% plot(sensingPositions(end,1),sensingPositions(end,2),'ok','MarkerFaceColor','k');
plot(sensingPositions(end,1),sensingPositions(end,2),'*g');







end

