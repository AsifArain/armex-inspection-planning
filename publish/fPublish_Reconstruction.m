function [ hRecon ] = fPublish_Reconstruction( map_recon,...
                                               cell_coords_x,...
                                               cell_coords_y,...
                                               map_env,...
                                               para_,...
                                               pubPara_)
%
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GAS CONCENTRATION MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

    
% %-- reconstruction 
% load([dir_.recon,para_.ReconstructionFile],...
%     'mean_map',...
%     'cell_coords_x',...
%     'cell_coords_y');
% map_recon = mean_map';

%dlmwrite('recon_publish.dat',map_recon,'delimiter',' ');

%--------------- select color map for gdm --------------
%gdm_colormap = flipud(colormap('hot')); % for main fig. 
gdm_colormap = flipud(hot(512)); % for main fig. 
%gdm_colormap = flipud(autumn(512)); % for main fig. 
%gdm_colormap = flipud(summer(512)); % for main fig. 
%gdm_colormap = flipud(colormap('summer')); % for insets
%-------------------------------------------------------

%------ color step size --------------------------------
%max_concentration = 10000; % inset gdm
%gt = dlmread([dir_.GroundTruthLogs,'map_con.dat']);
%max_concentration = max(gt(:))
%max_concentration = max(map_recon(:)); % main fig.
%max_concentration = 500;
max_concentration = min( [max(map_recon(:)),para_.PublishUpperPPM] ) % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);
%delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
%delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm

%------------------- nonlinear color code ---------------------  
if strcmp(pubPara_.color_scale,'nonlinear')
    %b = (exp(256/1))/255;
    %a = 256 /(log(b*256));

    x = 1:1:256;
    %epsilon = 1/(exp(1)-1); % standard
    %epsilon = 1/(exp(0.05)-1); % slower
    epsilon = 1/(exp(0.10)-1); % slower
    %epsilon = 1/(exp(0.15)-1); % slower
    %epsilon = 1/(exp(0.20)-1); % slower
    %epsilon = 1/(exp(0.90)-1); % slower
    %epsilon = 1/(exp(1.5)-1); % standard
    normlogsum_penalty = log(1+(abs(x)/epsilon));
    %normlogsum_penalty = 1.05.^(1+(abs(x)/epsilon));
    normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
    normlogsum_penalty = round(256*normlogsum_penalty);
    normlogsum_penalty(1)
    normlogsum_penalty(end)
    %figure; plot(x,normlogsum_penalty); % to test

    % publish nonlinear steps
    linear_percentage    = 10:10:100;
    linear_percentage_ind  = round((linear_percentage/100)*numel(x));
    nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
    nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
    nonlinear_percentage = round(nonlinear_percentage);
    %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
    disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
    %disp(disp_stng)
    disp(disp_stng)
end
%--------------------------------------------------------------

xyzPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);
colorPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);

reconColorR = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
reconColorG = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
reconColorB = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);

hRecon = [];
k = 1;
for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
    for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
        % if positive concentration

        if map_recon(i,j)>=10 %200%0

            if map_env(i,j)>0
            %------------------- useful trash -------------------
            %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);

            % --- if concentration is greater than 1000 ppm, its still 1000 ppm
            %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
            %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
            %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
            %-----------------------------------------------------

            %----------------------------------------------------
            % linear color code
            %----------------------------------------------------
            if strcmp(pubPara_.color_scale,'linear')
            % {
            linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number,:);                
            %}
            end
            %----------------------------------------------------

            %----------------------------------------------------
            % nonlinear color code
            %----------------------------------------------------
            if strcmp(pubPara_.color_scale,'nonlinear')
            % {
            linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
            %nonlinear_color_number = round( a*log(b*linear_color_number) );
            %nonlinear_color_number = linear_color_number;
            nonlinear_color_number = normlogsum_penalty(linear_color_number);
            %if nonlinear_color_number <1
            %    nonlinear_color_number = 1;
            %end
            col = gdm_colormap( nonlinear_color_number,:);                
            %}
            end
            %----------------------------------------------------

            %----- to avoid black borders of earch cell -----
            % {
            if sum(col) == 3
                col = col-(1e-10);
            end
            %}
            %----- plot ------------
            xyzPoints(k,:) = [cell_coords_x(i),cell_coords_y(j),0];                
            colorPoints(k,:) = col;                
            reconColorR(k) = round(col(1)*255);
            reconColorG(k) = round(col(2)*255);
            reconColorB(k) = round(col(3)*255);

            hRecon(k) = plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
                's',...
                'MarkerEdgeColor',col,...
                'MarkerFaceColor',col,...
                'MarkerSize',pubPara_.cell_marker_size);
            k = k+1;

            end
        end
    end
end

%{
map_recon = map_recon';
linear_color_number = round(map_recon(:)./delta)+1;

colR = round(gdm_colormap(linear_color_number,1)*255);
colG = round(gdm_colormap(linear_color_number,2)*255);
colB = round(gdm_colormap(linear_color_number,3)*255);

thisfile = sprintf('reconstructionColorR.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,colR,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorG.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,colG,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorB.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,colB,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

%colR(find(colR(:)~=255))

%pause
%ptCloud = pointCloud(xyzPoints,'Color',colorPoints)
%pcwrite(ptCloud,'testPlyFile','PLYFormat','binary')


test_recon = map_recon;%+rand(size(map_recon));
max(test_recon(:))
%test_recon = test_recon./(max(test_recon(:)));

%size(test_recon)
%pause

thisfile = sprintf('reconstructionColorR.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,reconColorR,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorG.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,reconColorG,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorB.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,reconColorB,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('reconstruction.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,test_recon,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('x_coord.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,cell_coords_x,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);

thisfile = sprintf('y_coord.dat');
fileID = fopen(thisfile,'wt'); fclose(fileID);
dlmwrite(thisfile,cell_coords_y,'-append','delimiter',' ');
fileID = fopen(thisfile,'a'); fclose(fileID);
%}
%pause

% --- concentration scale.
conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max_concentration);
%{
text(map_row-2,2,conc_text,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom',...
        'FontWeight','demi',...
        'FontSize',font_size_conc_scale,...
        'Color',[0,0,0],...
        'BackgroundColor',[1,1,1],...
        'Interpreter','latex',...
        'EdgeColor',[0,0,0],...
        'Margin',1,...
        'FontWeight','bold');
%}
end
