function [ tDomain ] = fCopyDomain( FoV,FoVfaces,SensorRange )
% 

tDomaini = tic; % initialize pre computational time.

f = FoVfaces.num; % for simplification.

%**********************************************************************
% COMPUTE DOMAIN
%**********************************************************************
fSaveDomain( f,...
             FoV,...
             FoVfaces,...
             SensorRange );

%**********************************************************************
% Total computation time.
%**********************************************************************
tDomain = 1e-4*(round(toc(tDomaini)*1e4));

fprintf(1,'\nComputation time: %0.4f sec \n',tDomain);




end



function [ ] = fSaveDomain( f,...
                            FoV,...
                            FoVfaces,...
                            SensorRange )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative angles
%FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.alf = wrapTo360(FoVfaces.alf);
FoVfaces.bay = wrapTo360(FoVfaces.bay);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];


%***********************************************************
%                     Sensor Domain
%***********************************************************
disp('.. Sensor Domain')

        
% For each cell of the map.
sD   = cell(1,f); % sensor domain (visible).
oD   = cell(1,f); % sensor's osbtacle domain (used to check the 
VObs = cell(1,f);               % visibility of "sD" for each cell in "oD"]

for i = 1
    %sensor = [0,0];
    % For each configuration of a cell.
    for j = 1:f
        
        
        this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',SensorRange.cell,FoV,j);
        this_dir  = sprintf('../sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/orientation_%03ddeg/',SensorRange.cell,FoV,j);

        
        this = load([this_dir,this_file],'sD','oD','VObs');

        sD{j}   = this.sD;
        oD{j}   = this.oD;
        VObs{j} = this.VObs;
    end
    
end
all_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_fovfaces%03d.mat',SensorRange.cell,FoV,f);
all_dir = sprintf('../sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',SensorRange.cell,FoV);
save([all_dir,all_file],'sD','oD','VObs','-V7.3');


end

%--------------------------- End of the document ------------------------------------
