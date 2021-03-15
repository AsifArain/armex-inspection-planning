function [ map_env,mapCut,robotStartCellInMap ] = fMapEnv( para_,dir_,visualize )
%
%


if ...
        strcmp(para_.ExperimentTitle,'corridor01') ||...
        strcmp(para_.ExperimentTitle,'corridor02') ||...
        strcmp(para_.ExperimentTitle,'corridor03') ||...
        strcmp(para_.ExperimentTitle,'corridor04') ||...
        strcmp(para_.ExperimentTitle,'corridor05') ||...
        strcmp(para_.ExperimentTitle,'corridor06') ||...
        strcmp(para_.ExperimentTitle,'corridor07') ||...
        strcmp(para_.ExperimentTitle,'corridor08') ||...
        strcmp(para_.ExperimentTitle,'corridor09') ||...
        strcmp(para_.ExperimentTitle,'corridor10') ||...
        strcmp(para_.ExperimentTitle,'corridor11') ||...
        strcmp(para_.ExperimentTitle,'corridor12') ||...
        strcmp(para_.ExperimentTitle,'corridor13') ||...
        strcmp(para_.ExperimentTitle,'corridor14') ||...
        strcmp(para_.ExperimentTitle,'corridor15') ||...
        strcmp(para_.ExperimentTitle,'corridor16') ||...
        strcmp(para_.ExperimentTitle,'corridor17') ||...
        strcmp(para_.ExperimentTitle,'corridor18') ||...
        strcmp(para_.ExperimentTitle,'corridor19') ||...
        strcmp(para_.ExperimentTitle,'corridor20') ||...
        strcmp(para_.ExperimentTitle,'corridor21')
        
    
    
    map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
    %map_env = 1-map_env;
    
    mapCut = [0,0,0,0];    
    robotStartCellInMap = [10,10];
    
elseif strcmp(para_.environment,'forest-area')
    
    map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
    %map_env = 1-map_env;
    mapCut = [0,0,0,0];    
    robotStartCellInMap = [10,10];    
    
elseif strcmp(para_.environment,'johnsson-metal')
    
    map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
    %map_env = 1-map_env;
    mapCut = [0,0,0,0];    
    robotStartCellInMap = [10,10];
    
elseif strcmp(para_.environment,'global-castings')
    
    map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
    %map_env = 1-map_env;
    mapCut = [0,0,0,0];    
    robotStartCellInMap = [10,10];            
    
elseif strcmp(para_.environment,'sample')
    
    map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
    %map_env = 1-map_env;
    mapCut = [0,0,0,0];    
    robotStartCellInMap = [50,10];     
    
elseif strcmp(para_.environment,'lab')
    
    map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
    %map_env = 1-map_env;
    mapCut = [0,0,0,0];
    robotStartCellInMap = [12,18];   
    
end


if visualize==1
    
    figure('name','robot start position in environment map'); 
    imshow(map_env','InitialMagnification',800); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
    
    plot(robotStartCellInMap(1),robotStartCellInMap(2),'*r')
    
    %pause
end




end