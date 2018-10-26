%% Heat & Mass Transfer HX Project %%
%% Intial Conditions

git = 900; % gas inlet temperature
ait = 200; % air inlet temperature
e = 0.8381; % effectivness
pt = 0.5*(10^-3); % plate thickness
gmf = 1.66; % gas mass flow rate
amf = 2.00; % air mas flow rate
gip = 160*(10^3); % gas inlet pressure;

%% Temperature calculations

cr = 0.857; 
got = git - e*(git - ait);
aot = ait + e*(git - ait)*cr;
tgm = (git + got)/2;
tam = (ait +aot)/2;

F = 0.8; % correction factor
tlm = ((ait-got)-(aot-git))/(log((ait-got)-(aot-git))); % log mean temperature
ctlm = F*tlm; % corrected log mean temperature
tlm2 = ((git-aot)-(got-ait))/(log((git-aot)-(got-ait)));
ctlm2 = F*tlm2;

% Assume operating conditions are at 1 atm since the change in pressure doesn't affect the results much%

%% Fluid properties

acp = 1079; % air cp value interpolated
gcp = 1852.6; % gas cp value interpolated
gpr = 0.8622; % Pr value for gas
apr = 0.6891; % Pr value for air
ug = 3.8425*(10^-5); % dynamic viscocity gas
ua = 2.6441*(10^-5); % dynamic viscocity air

pli = 0.36; % entrance pressure loss coefficient
plo = 0.42; % exit pressure loss coefficient

q = e*gmf*gcp*(git-ait); % total heat transfer
%% Area Properties

ba = 2.49*(10^-3); % plate spacing air
bg = 2.49*(10^-3); % plate spacing gas
da = 1.54*(10^-3); % hydraulic diameter air
dg = 1.54*(10^-3); % hydraulic diameter gas
ala = 941.6; % m2/m3 Area/Volume air
alg = 941.6; % m2/m3 Area/Volume gas
fta = 0.102*(10^-3); % fin thickness air
ftg = 0.102*(10^-3); % fin thickness gas
fat = 0.785; % fin area/total area air
%% ntu

tlm = 155.1724; % delta tm
NTU = 8.0816; % NTU
ntu_a = 2*NTU; % ntu air
ntu_g = 2*NTU; % ntu gas
%% Mass Velocity G

jh_g = 0.01979; % jh for gas
jh_a = 0.01785; % jh for air 
cf_a = 0.07254; % cf for air
cf_g = 0.08544 % cf for gas
dp_g = 9.05*(10^3); % differential pressure for gas
dp_a = 8.79*(10^3); % differential pressure for air

Ga = 11.6562; % mass velocity G for air
Gg = 10.3058; % mass velocity G for gas
%% Reynolds Number

Re_a = 530.79; % reynolds number for air
Re_g = 412.9; % reynolds number for gas
%% Fin Calculations

w_a = (ba/2)-fta; % fin height for air
w_g = (bg/2)-ftg; % fin height for gas

lf = 3.18*(10^-3); % fin offset length
per_a = 2*(lf+fta); % wetted perimeter for air
per_g = 2*(lf+ftg); % wetted perimeter for gas

areak_a = (lf*fta); % area of fin for air
areak_g = (lf*ftg); % area of fin for gas
%% Compute h's & n0's

h_a = jh_a*Ga*acp*(apr^(-2/3)); % h value for air
h_g = jh_g*Gg*gcp*(gpr^(-2/3)); % h value for gas

k = 18; % thermal conductivity of Iconell 625 alloy
m_a = ((h_a*per_a)/(k*areak_a))^(1/2); % m for efficiency for air
m_g = ((h_g*per_g)/(k*areak_g))^(1/2); % m for efficiency for gas

nf_a = tanh(m_a*w_a)/(m_a*w_a); % fin efficiency of air
nf_g = tanh(m_g*w_g)/(m_g*w_g); % fin efficiency of gas

no_a = 1-(1-nf_a)*fat; % fin efficiency actual of air
no_g = 1-(1-nf_g)*fat; % fin efficiency actual of air
%% Compute U of hot air

R_a = fta/k; % thermal resistivity of air
R_g = ftg/k; % thermal resistivity of gas

y = 1; % alpha_cold/alpha_hot ratio
U = 1/(((1/(no_a*h_a))+(R_a/no_a))+y*((1/(no_g*h_g))+(R_g/no_g))); % heat generated
%% Calculate core dimensions

cmin = gcp; % Cmin is our cp of gas
area_h = (NTU*cmin)/U; % area of hot (gas)
area_c = area_h; % area of cold (air) is equal to area of hot (gas)
A_oh = gmf/Gg; % free flow of gas (gas)
A_oc = amf/Ga; % free flow of air (air)
sig_h = alg*(dg/4); % sigma for gas
sig_c = ala*(da/4); % sigma for air
A_frh = A_oh/sig_h; % frontal area for gas
A_frc = A_oc/sig_c; % frontal area for air

%% Final Lengths

L1 = (bg*area_h)/(4*A_oh); % length
L2 = (ba*area_c)/(4*A_oc); % height
L3 = A_frh/L2; % width

%% Drop in pressures

dens_ih = 0.3010; % density in for hot
dens_ic = 0.7481; % density in for cold
dens_oh = 0.6028; % density out for hot
dens_oc = 0.3678; % density out for cold

mean_dens_h = 1/(0.5*((1/dens_ih)+(1/dens_oh))); % mean density for hot
mean_dens_c = 1/(0.5*((1/dens_ic)+(1/dens_oc))); % mean density for cold

kc = 0.36; % pressure loss coefficient 
ke = 0.42; % pressure loss coefficient 

delta_pg = ((Gg^2)/(2*dens_ih))*(1-(sig_h^2)+kc+(2*((dens_ih/dens_oh)-1))+((cf_g*L1*dens_ih/(dg/4))*(1/mean_dens_h))-(1-(sig_h^2)-ke)*((dens_ih/dens_oh)));
% delta p for gas
delta_pa = ((Ga^2)/(2*dens_ic))*(1-(sig_c^2)+kc+(2*((dens_ic/dens_oc)-1))+((cf_a*L2*dens_ic/(da/4))*(1/mean_dens_c))-(1-(sig_c^2)-ke)*((dens_ic/dens_oc)));
% delta p for air





