function [obstacle_mat]=load_occupancy_map_ASIF(varargin)
 n_params=length(varargin);
 curr_param=1;
 location='';
 offset=[0 0];
 resolution=1;
 rotation=0;
 markercolor='k';
 markersize=1.5;
 only_data=1;
 while curr_param<n_params
     param_name=varargin{curr_param};
     
     if strcmpi(param_name,'figure')==1
         curr_param=curr_param+1;
         location=varargin{curr_param};
     elseif strcmpi(param_name,'spec_file')==1
         curr_param=curr_param+1;
         load(varargin{curr_param});
     elseif strcmpi(param_name,'offset')==1
         curr_param=curr_param+1;
         offset=varargin{curr_param};
     elseif strcmpi(param_name,'resolution')==1
         curr_param=curr_param+1;
         resolution=varargin{curr_param};
     elseif strcmpi(param_name,'rotation')==1
         curr_param=curr_param+1;
         rotation=varargin{curr_param};
     elseif strcmpi(param_name,'markerfacecolor')==1
         curr_param=curr_param+1;
         markercolor=varargin{curr_param};
     elseif strcmpi(param_name,'markersize')==1
         curr_param=curr_param+1;
         markersize=varargin{curr_param};
     elseif strcmpi(param_name,'plot_map')==1
         curr_param=curr_param+1;
         only_data=varargin{curr_param};
     else
         curr_param=curr_param+1;
         
     end
 end
 
 if isempty(location)==1
     error('Specify a valid location');
 else
     I=imread(location);
     I_gray=rgb2gray(I)';
     map_input=double(I_gray);
     offset_x=offset(1);
     offset_y=offset(2);
     [rows,cols]=find(map_input==0);
     obstacle_mat=[cols(:),rows(:)];
     obstacle_mat=obstacle_mat.*resolution;
     obstacle_mat(:,1)=obstacle_mat(:,1)+offset_x;
     obstacle_mat(:,2)=obstacle_mat(:,2)+offset_y;
     RA=[cosd(rotation) -sind(rotation); sind(rotation) cosd(rotation)];
     obstacle_mat=obstacle_mat*RA;
     
     if only_data==1
         plot(obstacle_mat(:,1),obstacle_mat(:,2),'sk','markerfacecolor',markercolor,'markersize',markersize);
     end
 end
 
 
end

