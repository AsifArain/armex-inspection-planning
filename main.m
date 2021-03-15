
clear all; close all; clc;

disp('#____________________________________________________________________#')
disp('#                                                                    #')
disp('#                       ---  MAIN FILE  ----                         #')
disp('#                                                                    #')
disp('#                              ARMEx                                 #')
disp('#       ------------------------------------------------------       #')
disp('#                Autonomous Remote Methane Explorer                  #')
disp('#       ------------------------------------------------------       #')
disp('#                    One-tour Mission Strategy                       #')
disp('#                    Two-tour Mission Strategy                       #')
disp('#                      Human-expert Strategy                         #')
disp('#                                                                    #')
disp('#                                                                    #')
disp('#                       Author: Asif Arain                           #')
disp('#                                                                    #')
disp('#____________________________________________________________________#')


%{
    environment_ids:
         1 = 'corridor'
         2 = 'forest-area'
         3 = 'johnsson-metal'
         4 = 'global-castings'
         5 = 'sample'
         6 = 'lab'
         7 = 'prismaforum1'
         8 = 'prismaforum2'
         9 = 'prismaforum3'
        10 = 'prismaforum4'
        11 = 'prismaforum5'
        12 = 'smokebot_sample'
        13 = 'corridor-short'
        14 = 'sample-reconstruction'
        15 = 'sample-env-paper'
        16 = 'lab-ministers-demo-2019'
        17 = 'sample-thesis-rwl1'
        18 = 'sample-env-thesis'
        19 = 'thesis-cover'
%}

%{
    mission_strategies:
        1 = '1t-armex'
        2 = '2t-armex'
        3 = '2t-armex-icra16'
        4 = 'h-expert'
%}

environment_ids = 15; %18%5%17%5%14%5%13%11%5%11;
experiment_nums = [4];
mission_strategies = 2; %[1,2,5]

%-------------------------------------------------------------------
% Adding path directories
%-------------------------------------------------------------------
fAddDirPath( );

%-------------------------------------------------------------------
% Setting parameters related to planning algorithms and strategies
%-------------------------------------------------------------------
[ para_ ] = fAlgorithmParams( );

%-------------------------------------------------------------------
% Here we go with the experiments and the strategies
%-------------------------------------------------------------------
for env_id = environment_ids
    
    for exp_num = experiment_nums %exp_nums 
        
        close all
        
        [ para_ ] = fExperimentParams( env_id,exp_num,para_ );
        
        for mission_strategy = mission_strategies
            
            if     mission_strategy == 1
                f1tARMEx( para_ )

            elseif mission_strategy == 2
                f2tARMEx( para_ )

            elseif mission_strategy == 3
                f2tARMEx_ICRA16( para_ )
                
            elseif mission_strategy == 5
                fHumanExpert( para_ )
                
            end
        end
    end
end
