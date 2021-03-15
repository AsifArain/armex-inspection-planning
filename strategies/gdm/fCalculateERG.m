function [ ERG ] = fCalculateERG(nc,C,wV,D,conf_crossAngles_G,Conf_DistOrnGainV,para_)


% nc  = this.NumOfAllowedConf;%para_.NumberOfConf_SPPH;
% C   = selectedConf;
% G   = (alpha*beta)*(1/(nc*(nc-1)))*this.conf_crossAngles_G;
% Uc  = (alpha*(1-beta))/nc * ((gamma*sum(this.wV,1))+((1-gamma)*this.Conf_DistOrnGainV') );
% Ud  = (1-alpha)/nc*(1-(this.D/max(this.D)));
% U   = (Uc'+Ud);
% 
% thisERG   = C'*(G*C+U);

G   = (para_.alpha*para_.beta)*(1/(nc*(nc-1)))*conf_crossAngles_G;

Uc  = (para_.alpha*(1-para_.beta))/nc * ...
      ((para_.gamma*sum(wV,1))+((1-para_.gamma)*Conf_DistOrnGainV') );

Ud  = (1-para_.alpha)/nc*(1-(D/max(D)));
U   = (Uc'+Ud);

ERG = C'*(G*C+U);


    

%-- current ERG
%-------------------------------
% nc  = NumOfAllowedConf;%para_.NumberOfConf_SPPH;
% C   = reselectedConf;
% 
% G   = (para_.alpha*para_.beta)*(1/(nc*(nc-1)))*conf_crossAngles_G;
% 
% Uc  = (para_.alpha*(1-para_.beta))/nc *...
%       ((para_.gamma*sum(wV,1))+((1-para_.gamma)*Conf_DistOrnGainV') );
% 
% Ud  = (1-para_.alpha)/nc*(1-(D/max(D)));
% 
% U   = (Uc'+Ud);
% 
% thisERG = C'*(G*C+U);
%         
        
end
