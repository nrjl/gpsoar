function aout = vel_set(xin, vcore, Radi, theta)

zero_finder = ~(zeros(size(xin)) | xin);

xin = xin + zero_finder;

% Velocity field function
aout = abs(Radi/pi*vcore*sin(pi*xin/Radi)./(xin.*cos(theta)));
% aout = ones(size(theta));

% call = rref([Radi^7,   Radi^5,   Radi^3,   0;...
%              7*Radi^6, 5*Radi^4, 3*Radi^2, vcore;...
%              21*Radi^5,10*Radi^3,3*Radi,   0]);
% c1 = call(1,4); c2 = call(2,4); c3 = call(3,4);
% 
% aout = -(c1*(xin-Radi).^7 + c2*(xin-Radi).^5 + c3*(xin-Radi).^3)./(xin.*cos(theta));

aout = (aout.*~zero_finder) + (vcore*zero_finder);