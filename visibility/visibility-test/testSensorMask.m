%%
clear all; close all;
SensorRange.cell = 15;
FoV = 90;
load(['../SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)]);

hold on;
for i = -15:15
    for j = -15:15

        plot(i-0.5,j-0.5,'S', 'color',[0.9 0.9 0.9],'MarkerFaceColor',[0.9 0.9 0.9],'MarkerSize',10);
    end
end
% xlim([-1 5])
% ylim([-1 5])
axis equal
grid on

sensorDomain = sD{1};
plot(sensorDomain(:,1)-0.5,sensorDomain(:,2)-0.5,'+r')


obstacDomain = oD{1};
plot(obstacDomain(:,1)-0.45,obstacDomain(:,2)-0.45,'xb')

pause

VOBS = VObs{1};


radius = SensorRange.cell;
for i = 1:size(VOBS,1)
    i
    
    obsCell = obstacDomain(i,:)
    h1(1) = plot(obsCell(1)-0.5,obsCell(2)-0.5,'xk');
    
    ind = find(VOBS(i,:));
    
    visibleCells = sensorDomain(ind,:);
    
    size(visibleCells)
    
    if size(visibleCells)>0
        h1(2) = plot(visibleCells(:,1)-0.5,visibleCells(:,2)-0.5,'+r');
    end
    
    
    sensor = [0,0];
    r =0; c = 0;

    % Mark configuration
    plot(sensor(1)-0.5,sensor(2)-0.5,'o')
    
    % Draw FoV
    start_angle  = 45;     % first angle
    sector_angle = FoV;                     % increment in first angle
    [xx, yy] = SecDraw(start_angle, sector_angle, radius);
    h1(3) = plot(xx+r-0.5,yy+c-0.5,'color','r'); hold on;
    
    % Visible cells
%     vis = find(V(:,i));
%     if vis
%         [vis_r,vis_c] = ind2sub(size(map),vis);
%         h1(2) = plot(vis_r-0.5,vis_c-0.5,'+r');
%     end
%     pause
%     delete(h1(1));
%     if vis
%         delete(h1(2));
%     end
    pause
    delete(h1(1))
    if visibleCells
        delete(h1(2))
    end
    
end
