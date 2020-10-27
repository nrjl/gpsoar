function [x, y, z] = helical_path(ri, ro, zb, zt, n_cycles, n_samples)
%--------------------------------------------------------------------------
%
% FUNCTION:		helical_path
%
% PURPOSE:		Create helical path with given parameters
%
% SYNTAX:		[x, y, z] = helical_path(ri, ro, zb, zt, n_cyc,	n_samp)
%
% INPUTS:		ri		- Inner radius
%				ro		- Outer radius
%				zb		- z-coordinate of base
%				zt		- z-coordinate of top
%				n_cyc	- Number of cycles (loops)
%				n_samp	- Number of sample points (even angle distribution)
%
% OUTPUTS:		[x, y, z]	- Output position vectors
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		September 2009
%
% MODIFIED:     September 2009
%
% See also:		wind_GP3D
%--------------------------------------------------------------------------
% Parametrised by t. Sinusoidal variation of radius.
delta_z = zt-zb;

h_radius = @(t, n, ri, ro) (ro-ri)/2*(cos(t/n)-1) + ro;
t = linspace(0, 2*pi*n_cycles, n_samples);

x = h_radius(t, n_cycles, ri, ro).*cos(t);
y = h_radius(t, n_cycles, ri, ro).*sin(t);
z = zb + delta_z/(2*pi*n_cycles)*t;
