function h_fuse = SBXC_fuse( )

fuse_length = 25*3*25.4/1000;
r_fuse = [0, 4, 7, 9, 11, 11.5, 12, 12, 11.5, 11, 10, 9, 7, 6, 5.5, 5, ...
4.5, 4, 3.5, 3, 3, 2.5, 2, 1.5, 1, 0.5]/3*25.4/1000;

[ZZ, YY, XX] = cylinder(r_fuse); % Note that the primary axis is the x-axis

XX = 27.5*.0254 - XX*fuse_length; % Reverse axis (positive is forward)
ZZ = ZZ + 0.02;

h_fuse = surf(XX, YY, ZZ);

% [f, v, c] = surf2patch(XX, YY, ZZ); 
% v = Ceb*v' + X(10:12)*ones(1, size(v',2)); 
% fuse = patch(v(1,:), v(2,:), v(3,:), 'b');