function [ para_ ] = fExperimentParams( env_id,exp_num,para_ )
%
%fExperimentParams if used to set experimental parameters


fprintf(1,['\n\n',...
    '------------------------------------------\n',...
    '  Setting experiment parameters           \n',...
    '------------------------------------------\n'])

%---------------------------
% Experiment number
%---------------------------
%{
exp_num_corridor                = 10;%1;%1:6; %[1,3,5,6,7,8]; %[2,4,9]
exp_num_corridor_short          = 1;%1;%1:6; %[1,3,5,6,7,8]; %[2,4,9]
exp_num_forest                  = 1;
exp_num_johnsson                = 1;
exp_num_castings                = 1;%1:5;
exp_num_sample                  = 3; %4;%3%4%3;%1%:2;
exp_num_lab                     = 1;
exp_num_prismaforum1            = 1;
exp_num_prismaforum2            = 1;
exp_num_prismaforum3            = 1;
exp_num_prismaforum4            = 1;
exp_num_prismaforum5            = 4%4:10 %[4,5,6,8,9,10]; %10; [4,8]
exp_num_smokebot_sample         = 2; 
exp_num_sample_reconstruction   = 3;%1; 
exp_num_sample_env_paper        = 4; %3;
exp_num_sample_env_thesis       = 1;
exp_num_lab_ministers_demo_2019 = 1;
exp_num_sample_thesis_rwl1      = 1;
exp_num_thesis_cover            = 1;
%}

if     env_id == 1
    para_.environment = 'corridor';
    %exp_nums = exp_num_corridor;
    para_.FilePrefix1 = 'corridor';
elseif env_id == 2
    para_.environment = 'forest-area';
    %exp_nums = exp_num_forest;
    para_.FilePrefix1 = 'forest';
elseif env_id == 3
    para_.environment = 'johnsson-metal';
    %exp_nums = exp_num_johnsson;
    para_.FilePrefix1 = 'johnsson';
elseif env_id == 4
    para_.environment = 'global-castings';
    %exp_nums = exp_num_castings;
    para_.FilePrefix1 = 'global-casting';
elseif env_id == 5
    para_.environment = 'sample';
    %exp_nums = exp_num_sample;
    para_.FilePrefix1 = 'sample';
elseif env_id == 6
    para_.environment = 'lab';
    %exp_nums = exp_num_lab;
elseif env_id == 7
    para_.environment = 'prismaforum1';
    %exp_nums = exp_num_prismaforum1;
elseif env_id == 8
    para_.environment = 'prismaforum2';
    %exp_nums = exp_num_prismaforum2;
elseif env_id == 9
    para_.environment = 'prismaforum3';
    %exp_nums = exp_num_prismaforum3;    
elseif env_id == 10
    para_.environment = 'prismaforum4';
    %exp_nums = exp_num_prismaforum4;
elseif env_id == 11
    para_.environment = 'prismaforum5';
    %exp_nums = exp_num_prismaforum5;
    para_.FilePrefix1 = 'pf5';
elseif env_id == 12
    para_.environment = 'smokebot_sample';
    %exp_nums = exp_num_smokebot_sample;
elseif env_id == 13
    para_.environment = 'corridor-short';
    %exp_nums = exp_num_corridor_short;
    para_.FilePrefix1 = 'corridor-short';
elseif env_id == 14
    para_.environment = 'sample-reconstruction';
    %exp_nums = exp_num_sample_reconstruction;
    para_.FilePrefix1 = 'sample-reconstruction';
elseif env_id == 15
    para_.environment = 'sample-env-paper';
    %exp_nums = exp_num_sample_env_paper;
    para_.FilePrefix1 = 'sample-env-paper';
elseif env_id == 16
    para_.environment = 'lab-ministers-demo-2019';
    %exp_nums = exp_num_lab_ministers_demo_2019;
    para_.FilePrefix1 = 'lab-ministers-demo-2019';
elseif env_id == 17
    para_.environment = 'sample-thesis-rwl1';
    %exp_nums = exp_num_sample_thesis_rwl1;
    para_.FilePrefix1 = 'sample-thesis-rwl1';
elseif env_id == 18
    para_.environment = 'sample-env-thesis';
    %exp_nums = exp_num_sample_env_thesis;
    para_.FilePrefix1 = 'sample-env-thesis';
elseif env_id == 19
    para_.environment = 'thesis-cover';
    %exp_nums = exp_num_thesis_cover;
    para_.FilePrefix1 = 'thesis-cover';
end

%//////////////////////////////////////////////////
%
%   COMMON
%
%//////////////////////////////////////////////////

para_.ExperimentTitle = sprintf('%s-%02d',para_.environment,exp_num);
%para_.EnvironmentFigFile = sprintf('%s.dat',para_.environment);
para_.EnvironmentYamlFile = sprintf('%s.yaml',para_.environment);

para_.ROSServiceExecutedConfNum = '/executed_conf_num';
% para_.ROSServiceExecutedConfNum = '/plan_execution/conf_num';





    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'corridor')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////

    
    if exp_num == 1%2

        
        %para_.ExperimentTitle                 = 'corridor01';  % experiment title
        %para_.EnvironmentFigFile              = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'corridor.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        

    end

    if exp_num == 2%4
        % {
        %para_.ExperimentTitle               = 'corridor02';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 50; %20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 3%7
        % {
        %para_.ExperimentTitle               = 'corridor03';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 4%9
        % {
        %para_.ExperimentTitle               = 'corridor04';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 5%11
        % {
        %para_.ExperimentTitle               = 'corridor05';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 6%12
        % {
        %para_.ExperimentTitle               = 'corridor06';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 7%13
        % {
        %para_.ExperimentTitle               = 'corridor07';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 500; %20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 8%14
        % {
        %para_.ExperimentTitle               = 'corridor08';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end

    if exp_num == 9%15
        % {
        %para_.ExperimentTitle               = 'corridor09';  % experiment title
        %para_.EnvironmentFigFile            = 'corridor.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile           = 'corridor.yaml';
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        %para_.GroundTruthFile               = 'gasType1_simulation120';
        para_.GroundTruthTimeStamp          = 0;
        %para_.GroundTruthZSlice             = 5.5;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';

        %}
    end
    
    if exp_num == 10
        para_.SensingRangeM                 = 15;          % sensor range in meters
        para_.FieldOfViewDEG                = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp          = 0;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';
    end
    
end
    

    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'corridor-short')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////

    
    if exp_num == 1
        para_.SensingRangeM                 = 10;          % sensor range in meters
        para_.FieldOfViewDEG                = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM          = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM       = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp          = 0;
        para_.ReconTomographyCellSizeM      = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold  = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial      = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots            = 1;

        %para_.TomographyPlanType            = 'jfr17'; %'icra16'
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 

        para_.GroundTruthType               = 'manual';
        para_.HotspotsType                  = 'computed';
    end
    
    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
elseif strcmp(para_.environment,'forest-area')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////

    
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                 = 'forest-area01';  % experiment title
        %para_.EnvironmentFigFile              = 'forest-area.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'forest-area.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 500;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 2

        % {
        %para_.ExperimentTitle                 = 'forest-area02';  % experiment title
        %para_.EnvironmentFigFile              = 'forest-area.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'forest-area.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 500;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 3

        % {
        %para_.ExperimentTitle                 = 'forest-area03';  % experiment title
        %para_.EnvironmentFigFile              = 'forest-area.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'forest-area.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 500;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 4

        % {
        %para_.ExperimentTitle                 = 'forest-area04';  % experiment title
        %para_.EnvironmentFigFile              = 'forest-area.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'forest-area.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 500;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
   
end
    

    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'johnsson-metal')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                 = 'johnsson-metal01';  % experiment title
        %para_.EnvironmentFigFile              = 'johnsson-metal.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'johnsson-metal.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 2

        % {
        %para_.ExperimentTitle                 = 'johnsson-metal02';  % experiment title
        %para_.EnvironmentFigFile              = 'johnsson-metal.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'johnsson-metal.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 150;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 3

        % {
        %para_.ExperimentTitle                 = 'johnsson-metal03';  % experiment title
        %para_.EnvironmentFigFile              = 'johnsson-metal.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'johnsson-metal.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 100;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 4

        % {
        %para_.ExperimentTitle                 = 'johnsson-metal04';  % experiment title
        %para_.EnvironmentFigFile              = 'johnsson-metal.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'johnsson-metal.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 100;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'manual';    
        %}

    end
        
    if exp_num == 5

        % {
        %para_.ExperimentTitle                 = 'johnsson-metal05';  % experiment title
        %para_.EnvironmentFigFile              = 'johnsson-metal.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'johnsson-metal.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 100;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'manual';    
        %}

    end
    
end

    
    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'global-castings')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                 = 'global-castings01';  % experiment title
        %para_.EnvironmentFigFile              = 'global-castings.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'global-castings.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 200;%20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 2

        % {
        %para_.ExperimentTitle                 = 'global-castings02';  % experiment title
        %para_.EnvironmentFigFile              = 'global-castings.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'global-castings.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 500;%20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 3

        % {
        %para_.ExperimentTitle                 = 'global-castings03';  % experiment title
        %para_.EnvironmentFigFile              = 'global-castings.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'global-castings.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 200;%20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 4

        % {
        %para_.ExperimentTitle                 = 'global-castings04';  % experiment title
        %para_.EnvironmentFigFile              = 'global-castings.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'global-castings.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 200;%20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 5

        % {
        %para_.ExperimentTitle                 = 'global-castings05';  % experiment title
        %para_.EnvironmentFigFile              = 'global-castings.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'global-castings.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 200;%20;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}
      
    end
    
end    


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'sample')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                 = 'sample01';  % experiment title
        %para_.EnvironmentFigFile              = 'sample.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'sample.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
    if exp_num == 2

        % {
        %para_.ExperimentTitle                 = 'sample02';  % experiment title
        %para_.EnvironmentFigFile              = 'sample.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'sample.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'manual';    
        %}

    end
    
    if exp_num == 3

        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    

    end
    
    if exp_num == 4

        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    

    end
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'smokebot_sample')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    
    if exp_num == 1
        
        para_.SensingRangeM                   = 20;          % sensor range in meters
        para_.FieldOfViewDEG                  = 180;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';
        
    end
        
    if exp_num == 2
        
        para_.SensingRangeM                   = 20;          % sensor range in meters
        para_.FieldOfViewDEG                  = 45;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';
        
    end
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'sample-reconstruction')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    
    if exp_num == 1
        
        para_.SensingRangeM                   = 15;   % sensor range in meters
        para_.FieldOfViewDEG                  = 90;   % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;    % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM; % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;  % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98; % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;    % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';
    end
    
    if exp_num == 2
        
        para_.SensingRangeM                   = 30;   % sensor range in meters
        para_.FieldOfViewDEG                  = 90;   % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;    % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM; % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;  % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98; % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;    % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';
    end
    
    if exp_num == 3
        
        para_.SensingRangeM                   = 30;   % sensor range in meters
        para_.FieldOfViewDEG                  = 90;   % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;    % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM; % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;  % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98; % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;    % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'manual';
    end
        
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'sample-env-paper')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
       
    if exp_num == 1

        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';   

    end    
    
    if exp_num == 2

        para_.SensingRangeM                   = 30;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    

    end    
    
    if exp_num == 3

        para_.SensingRangeM                   = 30;          % sensor range in meters
        para_.FieldOfViewDEG                  = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10; %2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    

    end
    
    if exp_num == 4

        para_.SensingRangeM                   = 30;          % sensor range in meters
        para_.FieldOfViewDEG                  = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10; %2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    

    end    
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'sample-env-thesis')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
       
    if exp_num == 1

        para_.SensingRangeM                   = 30;          % sensor range in meters
        para_.FieldOfViewDEG                  = 270;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';   

    end
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'sample-thesis-rwl1')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
       
    if exp_num == 1

        para_.SensingRangeM                   = 5;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';   

    end    
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'thesis-cover')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
       
    if exp_num == 1

        para_.SensingRangeM                   = 5;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';   

    end    
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'lab')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                 = 'lab01';  % experiment title
        %para_.EnvironmentFigFile              = 'lab.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'lab.yaml';
        para_.SensingRangeM                   = 15;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'lab-ministers-demo-2019')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    if exp_num == 1
        
        para_.SensingRangeM                   = 10;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        
    end 

    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'prismaforum1')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                 = 'prismaforum01';  % experiment title
        %para_.EnvironmentFigFile              = 'prismaforum.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile             = 'prismaforum.yaml';
        para_.SensingRangeM                   = 10;          % sensor range in meters
        para_.FieldOfViewDEG                  = 90;          % sensor field of view in degrees
        para_.EnvironmentCellSizeM            = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM         = para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp            = 0;
        para_.ReconTomographyCellSizeM        = 1.0;   % cell size of ragt gdm
        para_.EnvironmentGray2BinThreshold    = 0.98;   % threshold to convert gray scal map to binary map.
        para_.ReconDetectionArtificial        = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots              = 0;
        para_.highConcentrationThreshold_PPM  = 10;%2500; %5000; 
        para_.GroundTruthType                 = 'manual';
        para_.HotspotsType                    = 'computed';    
        %}

    end
    
end    
    

    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'prismaforum2')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
       
       
    
    if exp_num == 1

        % {
        %para_.ExperimentTitle                   = 'prismaforum201';  % experiment title
        %para_.EnvironmentFigFile               = 'prismaforum2.dat'; % indoor_tecknikhuset
        %para_.EnvironmentYamlFile              = 'prismaforum2.yaml';
        %para_.EnvironmentYamlFile               = sprintf('%s.yaml',para_.environment);
        
        
        %para_.EnvironmentGray2BinThreshold      = 0.98;   % threshold to convert gray scal map to binary map.        
        para_.SensingRangeM                     = 10;          % sensor range in meters
        para_.FieldOfViewDEG                    = 90;          % sensor field of view in degrees
        
        %para_.EnvironmentCellSizeM              = 1;           % map cell size in meters
        para_.ReconDetectionCellSizeM           = 1; %para_.EnvironmentCellSizeM;   % cell size of ragd gdm
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 1.0;   % cell size of ragt gdm
        para_.ReconstructionMapSize             = 'environment-map';
        %para_.ReconstructionMapSize            = 'measurement-beams';
        
        
        
        para_.ReconDetectionArtificial          = 0;      % use artificial generated ragd gdm
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;%2500; %5000; 
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';    
        %}

    end    
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'prismaforum3')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
    
    
    if exp_num == 1

        para_.SensingRangeM                     = 10;
        para_.FieldOfViewDEG                    = 90;
        para_.ReconDetectionCellSizeM           = 1;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 1.0;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        
    end    


end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'prismaforum4')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
       
       
    
    if exp_num == 1

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 180;
        para_.ReconDetectionCellSizeM           = 1;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 1.0;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
    end    
    
    
end


    %///////////////////////////////////////////////////////////////////////////
    %                                                                         //
if strcmp(para_.environment,'prismaforum5')
    %                                                                         //
    %///////////////////////////////////////////////////////////////////////////
              
    
    
    
    
    if exp_num == 1

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 180;
        para_.ReconDetectionCellSizeM           = 1;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 1.0;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
    end
    
    if exp_num == 2

        para_.SensingRangeM                     = 10;
        para_.FieldOfViewDEG                    = 180;
        para_.ReconDetectionCellSizeM           = 1;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 1.0;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
    end    
    
    if exp_num == 3

        para_.SensingRangeM                     = 10;
        para_.FieldOfViewDEG                    = 90;
        para_.ReconDetectionCellSizeM           = 1;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 1.0;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
    end
    
    
    if exp_num == 4

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 100;%10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'HanFan';
    end
    
    if exp_num == 5

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'HanFan';
    end
    
    
    if exp_num == 6

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'VictorHernandez';
    end
    
    if exp_num == 7

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10;
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'VictorHernandez';
    end
    
    
    if exp_num == 8

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 10; %50; %10
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'VictorHernandez';
    end
    
    if exp_num == 9

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 20; %210; %10
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'Erik';
    end
    
    if exp_num == 10

        para_.SensingRangeM                     = 15;
        para_.FieldOfViewDEG                    = 270;
        para_.ReconDetectionCellSizeM           = 0.5;
        para_.GroundTruthTimeStamp              = 0;
        para_.ReconTomographyCellSizeM          = 0.5;
        para_.ReconstructionMapSize             = 'environment-map';
        para_.ReconDetectionArtificial          = 0;
        para_.ArtificialHotspots                = 0;
        para_.highConcentrationThreshold_PPM    = 50; %10; %50; %20; %210; %10
        para_.GroundTruthType                   = 'manual';
        para_.HotspotsType                      = 'computed';
        para_.HumanExpert                       = 'Achim';
    end
    
end


%------------------------------------------

% FILE PREFIX2
para_.FilePrefix2 = sprintf('%s-%02d',para_.FilePrefix1,exp_num);



end