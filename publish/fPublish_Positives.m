function fPlotPositives( map_env,...
                         map_recon,...
                         cell_coords_x,...
                         cell_coords_y,...
                         sc_estimated,...
                         sc_true,...
                         tf_positives_association_ind,...
                         para_,...
                         dir_ )

figure('name','Estimated gas sources - true false positives on recon map'); 
imshow(map_env','InitialMagnification',350); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])

%------------------------------------------------------------
cell_marker_size = 4.85; %7; %4.85;
color_scale = 'linear'; %
%--------------- select color map for gdm --------------
gdm_colormap = flipud(autumn(512)); % for main fig. 

%------ color step size --------------------------------
max_concentration = min( [max(map_recon(:)),para_.PublishUpperPPM] ); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);

%------------------- nonlinear color code ---------------------  
if strcmp(color_scale,'nonlinear')
    %b = (exp(256/1))/255;
    %a = 256 /(log(b*256));

    x = 1:1:256;
    epsilon = 1/(exp(0.10)-1); % slower
    normlogsum_penalty = log(1+(abs(x)/epsilon));
    normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
    normlogsum_penalty = round(256*normlogsum_penalty);
    normlogsum_penalty(1)
    normlogsum_penalty(end)
    %figure; plot(x,normlogsum_penalty); % to test

    % publish nonlinear steps
    linear_percentage = 10:10:100;
    linear_percentage_ind = round((linear_percentage/100)*numel(x));
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

k = 1;
for i = 1:numel(cell_coords_x)
    for j = 1:numel(cell_coords_y)
        % if positive concentration
        if map_recon(i,j)>=500%50%0
            %if map_env(i,j)>0
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
            if strcmp(color_scale,'linear')
            % {
            linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number,:);                
            %}
            end
            %----------------------------------------------------

            %----------------------------------------------------
            % nonlinear color code
            %----------------------------------------------------
            if strcmp(color_scale,'nonlinear')
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

            plot(cell_coords_x(i)+0.5,cell_coords_y(j)+0.5,...
                's',...
                'MarkerEdgeColor',col,...
                'MarkerFaceColor',col,...
                'MarkerSize',cell_marker_size);
            k = k+1;
            %end
        end
    end
end
pause(1)
%--------------------------
%--- estimated sources
%--------------------------
for i = 1:size(sc_estimated,1)    
    plot(sc_estimated(i,1),sc_estimated(i,2),...
         'o',...
         'color',[0.1,0.1,0.1],...
         'MarkerFaceColor','r',...
         'MarkerSize',6.0); %7
end
%--------------------------
%--- true sources
%--------------------------
for i = 1:size(sc_true,1)    
    plot(sc_true(i,1),sc_true(i,2),...
         'o',...
         'color',[0.1,0.1,0.1],...
         'MarkerFaceColor','b',...
         'MarkerSize',4.0); %7
end    
%--------------------------
%-- true positives
%--------------------------
%{
colors = linspecer(numel(tf_positives_association_ind));
for i = 1:size(sc_true,1)
    pts = sc_estimated(tf_positives_association_ind==i,:);
    if size(pts,1)>0            
        plot(pts(:,1),pts(:,2),'o','color',colors(i,:),'MarkerSize',6.5); %3
        plot(sc_true(i,1),sc_true(i,2),'o','color',colors(i,:),'MarkerSize',6.5); %3
    end
end
%}
%{
colors = linspecer(numel(tf_positives_association_ind));
for i = 1:size(sc_true,1)        
    pts = [sc_estimated(tf_positives_association_ind==i,:);sc_true(i,:)];
    %pts = [true_positives(true_positives_association_ind==i,:);sc_true(i,:)];
    plot(pts(:,1),pts(:,2),...
         'o','color',colors(i,:),'MarkerSize',6.5); %3
end
%}
%{
%colors = linspecer(numel(tf_positives_association_ind));
colors = linspecer(size(sc_true,1));
for i = 1:size(sc_true,1)
    pts = true_positives(true_positives_association_ind==i,:);

    plot(pts(:,1),pts(:,2),...
         'o',...
         'color','k',...
         'MarkerFaceColor','g',...
         'MarkerSize',4.5); %3

    plot(pts(:,1),pts(:,2),...
         'v','color',colors(i,:),'MarkerSize',6.5); %3

    plot(sc_true(i,1),sc_true(i,2),...
         'o','color',colors(i,:),'MarkerSize',4.5); %3
end
%}
%--------------------------
%-- false positives
%--------------------------
pts = sc_estimated(tf_positives_association_ind==0,:);
plot(pts(:,1),pts(:,2),...
     'x','color','k','MarkerSize',6.5); %3
%{
plot(false_positives(:,1),false_positives(:,2),...
         'o',...
         'color','k',...
         'MarkerFaceColor','g',...
         'MarkerSize',4.5); %3         
plot(false_positives(:,1),false_positives(:,2),...
     'x','color','k','MarkerSize',6.5); %3
%}
%--------------------------
%-- LEGENDS
%--------------------------
% true-positive
plot(160,48,'o','color',[0.1,0.1,0.1],'MarkerFaceColor','r','MarkerSize',4.5);
% false-positive
plot(160,43,'o','color',[0.1,0.1,0.1],'MarkerFaceColor','r','MarkerSize',4.5);
plot(160,43,'x','color','k','MarkerSize',6.5); %3
% gas sources
plot(160,38,'o','color',[0.1,0.1,0.1],'MarkerFaceColor','b','MarkerSize',4.5);

%--------------------------
% SCALE
%--------------------------
% scale_pt1 = [090,40]
% scale_pt2 = [100,40]
scale_pt1 = [(size(map_env,1)/2)-5,40]
scale_pt2 = [(size(map_env,1)/2)+5,40]

plot([scale_pt1(1),scale_pt2(1)],[scale_pt1(2)-0,scale_pt2(2)+0],'-k','LineWidth',0.5)
plot([scale_pt1(1),scale_pt1(1)],[scale_pt1(2)-1,scale_pt1(2)+1],'-k','LineWidth',0.5)
plot([scale_pt2(1),scale_pt2(1)],[scale_pt2(2)-1,scale_pt2(2)+1],'-k','LineWidth',0.5)


pause(2)
%-- save
filename = sprintf('%s%s-positives',dir_.evaluation,para_.FilePrefix);
%print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500');
pause(1)

%-- save to a common folder
dir_fig = '/media/a/wrk/planners/spp-tomography/arain-spp-rem-jfr2017/fig/';
filename = sprintf('%s%s-positives',dir_fig,para_.FilePrefix);
%print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500');
pause(1)

end