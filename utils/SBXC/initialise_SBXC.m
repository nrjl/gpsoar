global g
global PLANE_AERO PLANE_PARAM

g = 9.81;

x_cg = -0.5*290;
y_cg = 0;
z_cg = 0;
[SBXC_aero, SBXC_param] = SBXC_def(x_cg, y_cg, z_cg);
plane_properties(SBXC_aero, SBXC_param);
PLANE_AERO = SBXC_aero;
PLANE_PARAM = SBXC_param;