function measurement_path = GetMeasurementPath(start_x,end_x,start_y,end_y,cell_coords_x,cell_coords_y)

    cell_size = cell_coords_x(2) - cell_coords_x(1);

    cell_start_x = floor((start_x-cell_coords_x(1))./cell_size)+1;
    cell_end_x = floor((end_x-cell_coords_x(1))./cell_size)+1;
    
    cell_start_y = floor((start_y-cell_coords_y(1))./cell_size)+1;
    cell_end_y = floor((end_y-cell_coords_y(1))./cell_size)+1;

    if (cell_start_x <= cell_end_x)
        x_intersections = cell_start_x+1:1:cell_end_x;
    else
        x_intersections = cell_start_x:-1:cell_end_x+1;
    end
    
    if (cell_start_y <= cell_end_y)
        y_intersections = cell_start_y+1:1:cell_end_y;
    else
        y_intersections = cell_start_y:-1:cell_end_y+1;
    end

    %intersection with cells on x axis
    intersection_points_x = zeros(length(x_intersections),3);
    if ~isempty(x_intersections)
        intersection_points_x(:,1) = cell_coords_x(x_intersections);
        intersection_points_x(:,3) = (intersection_points_x(:,1) - start_x)./(end_x - start_x);
        intersection_points_x(:,2) = start_y + (end_y - start_y)*intersection_points_x(:,3);
    end

    %intersection with cells on y axis
    intersection_points_y = zeros(length(y_intersections),3);
    if ~isempty(y_intersections)
        intersection_points_y(:,2) = cell_coords_y(y_intersections);
        intersection_points_y(:,3) = (intersection_points_y(:,2) - start_y)./(end_y - start_y);
        intersection_points_y(:,1) = start_x + (end_x - start_x)*intersection_points_y(:,3);
    end

    intersection_points = [intersection_points_x ; intersection_points_y];
    %add start and end point
    intersection_points = [start_x start_y 0; intersection_points ; end_x end_y 1];
    intersection_points = sortrows(intersection_points,3);

    %the distance from the start of the beam to the intersection of each cell is the 4th column (t) multiplied by the
    %length of the beam
    dist = sqrt((end_x - start_x)^2 + (end_y - start_y)^2);
    distances = intersection_points(:,3) * dist;
    distances = diff(distances);
    
    cells = zeros(length(distances)-1,2);
    for i = 1:size(intersection_points,1)-1
        mid_x = (intersection_points(i,1) + intersection_points(i+1,1))/2;
        mid_y = (intersection_points(i,2) + intersection_points(i+1,2))/2;
        cells(i,1) = floor((mid_x-cell_coords_x(1))./cell_size)+1; 
        cells(i,2) = floor((mid_y-cell_coords_y(1))./cell_size)+1;
    end
  
    indexes = length(cell_coords_y)*(cells(:,1)-1) + cells(:,2);
    measurement_path = zeros(1,(length(cell_coords_x))*(length(cell_coords_y)));
    measurement_path(indexes) = distances;
end