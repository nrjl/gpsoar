function plot3mat(xyz, options)

in_size = size(xyz);

if ndims(xyz) == 2
    
    if in_size(1) == 3
        plot3(xyz(1,:), xyz(2,:), xyz(3,:))
    elseif in_size(2) == 3
        plot3(xyz(:,1), xyz(:,2), xyz(:,3))
    end
else
    print(['Error: Incorrect dimensions, input must be 2-dimensional ', ...
        'array with 3 rows or columns.'])
end    