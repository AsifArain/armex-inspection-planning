function [ hSymbol ] = fConfSymobls(confs,conf_num,pubPara_)

for j = 1:conf_num
        
    if     confs(j,7) == 0
        conf_color = pubPara_.colorExploration;
    elseif confs(j,7) == 1
        conf_color = pubPara_.colorExploitation;
    end

    % -- drat sumbolic field of view
    start_angle  = confs(j,3)-(confs(j,4)/2); %FoVfaces.ang(nCONF(i),1);
    end___angle  = confs(j,3)+(confs(j,4)/2); %start_angle+FoV;
    r            = pubPara_.symbolic_fov_length;
    x1           = r*cosd(start_angle);
    y1           = r*sind(start_angle);
    x2           = r*cosd(end___angle);
    y2           = r*sind(end___angle);
    hSymbol.line1(j) = plot([confs(j,1)-0.5,x1+confs(j,1)-0.5],...
                            [confs(j,2)-0.5,y1+confs(j,2)-0.5],...
                            'color',conf_color,'LineWidth',pubPara_.line_width_symbolic_fov);
    hSymbol.line2(j) = plot([confs(j,1)-0.5,x2+confs(j,1)-0.5],...
                            [confs(j,2)-0.5,y2+confs(j,2)-0.5],...
                            'color',conf_color,'LineWidth',pubPara_.line_width_symbolic_fov); 
    hSymbol.posit(j) = plot(confs(j,1)-0.5,confs(j,2)-0.5,...
                            'o','color',conf_color,'MarkerFaceColor',conf_color);  

    % Draw arrow
    hSymbol.arrow(j) = quiver(confs(j,1)-0.5,confs(j,2)-0.5,...
                              pubPara_.arrow_length*cosd(confs(j,3)),...
                              pubPara_.arrow_length*sind(confs(j,3)),...
                              'LineWidth',pubPara_.line_width_arrow,...
                              'color',conf_color,...
                              'MaxHeadSize',pubPara_.arrow_head_size,... %'Marker','o',...
                              'MarkerSize',pubPara_.marker_size_arrow);

end

end
