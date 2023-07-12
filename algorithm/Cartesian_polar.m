function coordinate = cartesian_polar(cartesian_coordinate,polar_coordinate )
% Conversion of Cartesian and hyperspherical coordinates
% When cartesian coordinates are input and hyperspherical coordinates are empty, they are converted to hyperspherical coordinates. 
% The opposite is similar.
% 
% cartesian_polar(cartesian_coordinate,[]) cartesian_coordinate are
% converted to coordinate,which are hyperspherical coordinates.

% Convert to hyperspherical coordinates
if ~isempty(cartesian_coordinate)
    m=length(cartesian_coordinate);
    coordinate=zeros(1,m);
    up=sum(cartesian_coordinate.^2);
    coordinate(1)=sqrt(up);
    
    for i=2:m
        up=sum(cartesian_coordinate(i:m).^2);
        coordinate(i)=atan(sqrt(up)/cartesian_coordinate(i-1));
        if isnan(coordinate(i))
            coordinate(i)=pi/2;
        end
    end
   % Convert to Cartesian coordinates
elseif ~isempty(polar_coordinate)
    m=length(polar_coordinate);
    coordinate=zeros(1,m);
    r=polar_coordinate(1);
    coordinate(1)=r*cos(polar_coordinate(2));
    
    for i=2:m-1
        coordinate(i)=coordinate(i-1)*tan(polar_coordinate(i))*cos(polar_coordinate(i+1));
    end
    coordinate(m)=coordinate(m-1)*tan(polar_coordinate(m));
    
end

end

