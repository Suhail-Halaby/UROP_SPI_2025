function [] = SimuPilot_Initializer ()

gravity = 9.81;
aspd_min = 5;
mass = 4;
scaler = 1;

RLL2SRV_TCONST = 0.2;
RLL2SRV_MAX = deg2rad(20);

RLL2SRV_P = 2.5;
RLL2SRV_I = 0.8;
RLL2SRV_D = 0.5;


RLL_K_P = (RLL2SRV_P - RLL2SRV_I * RLL2SRV_TCONST) * RLL2SRV_TCONST - RLL2SRV_D;
RLL_K_I = RLL2SRV_I*RLL2SRV_TCONST;


PTCH2SRV_TCONST = 0.2;
PTCH2SRV_RMAX_UP = deg2rad(10);
PTCH2SRV_RMAX_DN = -deg2rad(10);
PTCH2SRV_RLL = 0;

PTCH2SRV_P = -1.0;
PTCH2SRV_I = -0.4;
PTCH2SRV_D = -0.2;


PTCH_K_P = (PTCH2SRV_P - PTCH2SRV_I * PTCH2SRV_TCONST) * PTCH2SRV_TCONST - PTCH2SRV_D;
PTCH_K_I = PTCH2SRV_I*PTCH2SRV_TCONST;


YAW2SRV_RLL = 1.5;

YAW2SRV_SLIP = 1;
YAW2SRV_INT = 0.5;
YAW2SRV_DAMP = 0.05;


K_V = 0.1;
K_H = 0.1;
K_D = 0.05; 
K_T = 0.05;



H_INIT = 0 ;
V_INIT = 15;

vars = who;
for k = 1:numel(vars)
    if evalin('base', sprintf('exist(''%s'',''var'')', vars{k}))
        warning('Base workspace already has variable %s; it will be overwritten.', vars{k});
    end
    assignin('base', vars{k}, eval(vars{k}));
end

