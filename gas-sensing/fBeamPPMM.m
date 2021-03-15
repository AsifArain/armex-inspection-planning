function [ppm_m,...
          beam_voxel____end,...
          beam_voxel_end_xy,...
          beam_voxel___PPMM] = fBeamPPMM(beam_voxel___list,...
                                         beam_voxel_length,...
                                         beam_voxel_points,...
                                         map___environment,...
                                         map_concentration)


[map_row,map_col] = size(map___environment);
%map___environment = ones(size(map___environment));

% -- measurement of integral concentration on grid map with obstacles

% -- end voxel of the beam
% beam_voxel____end = [];
% beam_voxel_end_xy = [];
beam_voxel____end = [beam_voxel___list(1,1), beam_voxel___list(2,1)];
beam_voxel_end_xy = [beam_voxel_points(1,1), beam_voxel_points(2,1)];

ppm_m = 0;

beam_voxel___PPMM = [];

% beam_voxel___list
for i = 1:size(beam_voxel___list,2)
    
    %beam_voxel___list(1,i)
    %beam_voxel___list(2,i)
    
    % if voxel is out of map
    if beam_voxel___list(1,i) <= 0 || beam_voxel___list(1,i) > map_row ||...
       beam_voxel___list(2,i) <= 0 || beam_voxel___list(2,i) > map_col
    
        %disp('--- broken in first condition');
        break;
        
    end
    
    % if voxel is an occupied cell    
    if map___environment(beam_voxel___list(1,i),beam_voxel___list(2,i)) == 0
        
        %disp('--- broken in second condition');
        break;
    end
    
    % -- update end voxel
    beam_voxel____end = [beam_voxel___list(1,i), beam_voxel___list(2,i)];
    beam_voxel_end_xy = [beam_voxel_points(1,i), beam_voxel_points(2,i)];
    
    
    % -- update voxel list PPMM
    beam_voxel___PPMM = [beam_voxel___PPMM, [beam_voxel___list(1,i);beam_voxel___list(2,i)]];
    
    %plot(beam_voxel__end_x,beam_voxel__end_y,'*r')
    
    %plot(beam_voxel__end_x+0.5,beam_voxel__end_y+0.5,'sr')
        
    % -- update intergrated concentration
    %disp('-- this concentration')
    %map_concentration(beam_voxel___list(1,i),beam_voxel___list(2,i))
    this_ppm_m = beam_voxel_length(i)*map_concentration(beam_voxel___list(1,i),beam_voxel___list(2,i));
    
    ppm_m = ppm_m + this_ppm_m;
    
    %ppm_m = ppm_m + beam_voxel_length(i)*...
    %                map_concentration(beam_voxel___list(1,i),beam_voxel___list(2,i));
    
end



%size(beam_voxel___list,2)
%i
%ppm_m

end