%%% E (Pyramidal), I (PV) & I2 (CCK) 

%-----------------------------------
% number of E- and I-cells:
NE=80;                         % total number of pyramidal neurons in the network
perc=31;                       % (empirical) percentge of pyr neurons with FAST adaptation
num_e_a=round(perc*NE/100);    % number of pyramidal neurons with FAST adaptation

num_e=NE-num_e_a;              % number of pyramidal neurons with SLOW adaptation
num_i=10;
num_i2=20;
%-----------------------------------



%-----------------------------------
%  maximal conductance of the muscarinic potasium channel (I_m)

% for slow adapters
% g_sa_3pmc=2.52;  % mS/cm2;  L3 ACC --> PMC    
% g_sa_5pmc=2.66;  % mS/cm2;  L5 ACC --> PMC
% g_sa_3amy=2.54;  % mS/cm2;  L3 ACC --> AMY 
% g_sa_5amy=2;     % mS/cm2;  L5 ACC --> AMY 

g_sa=2.66;

% for fast adapters
g_fa=20;  % mS/cm2
%-----------------------------------



%-----------------------------------
% density of synaptic connections: 
p_ee=1.0; p_ei=1.0; p_ie=1.0; p_ii=1.0;
p_i2i=1.0; p_i2e=1.0; p_ei2=1.0; p_ii2=1.0; p_i2i2=1.0; 
p_eaea=1.0; p_eea=1.0; p_eae=1.0; p_eai=1.0; p_iea=1.0; p_i2ea=1.0; p_eai2=1.0; 
%-----------------------------------



%-----------------------------------
% rise and decay time constants associated with synapses:
%tau_r_e=0.1*ones(num_e,1); tau_d_e=4.6*ones(num_e,1);     % L3 ACC --> PMC  (EPSCs for AMPA, Maya's data); single channel based on the model
tau_r_e=0.1*ones(num_e,1); tau_d_e=4.2*ones(num_e,1);     % L5 ACC --> PMC 
%tau_r_e=0.1*ones(num_e,1); tau_d_e=4.9*ones(num_e,1);     % L3 ACC --> AMY 
%tau_r_e=0.1*ones(num_e,1); tau_d_e=4.2*ones(num_e,1);     % L5 ACC --> AMY  

%tau_r_e_a=0.1*ones(num_e_a,1); tau_d_e_a=4.6*ones(num_e_a,1);     % L3 ACC --> PMC  (EPSCs for AMPA, Maya's data); single channel based on the model
tau_r_e_a=0.1*ones(num_e_a,1); tau_d_e_a=4.2*ones(num_e_a,1);     % L5 ACC --> PMC 
%tau_r_e_a=0.1*ones(num_e_a,1); tau_d_e_a=4.9*ones(num_e_a,1);     % L3 ACC --> AMY 
%tau_r_e_a=0.1*ones(num_e_a,1); tau_d_e_a=4.2*ones(num_e_a,1);     % L5 ACC --> AMY 

tau_r_i=0.3*ones(num_i,1); tau_d_i=8*ones(num_i,1);         % fast IPSCs (PV), Medalla et al. 2017 (decay times= 2-14 ms); take the mean value  
tau_r_i2=0.3*ones(num_i2,1); tau_d_i2=42*ones(num_i2,1);    % slow IPSCs (CCK), Medalla et al. 2017 (decay times >14 to 70 ms); take the mean value       
%-----------------------------------



%-----------------------------------
% this section has been added by Sara: membrane capacitance for pyr. neurons (taken from Maya's data).
% membrane capacitance = membrane time constant / membrane resistance
% SurfArea = 559.06*10^(-8);              % [cm2]    L3 ACC --> PMC
% Ce = (30.2*10^(-3)/193.1)/SurfArea;     % [cm2]    L3 ACC --> PMC

SurfArea = 931.92 *10^(-8);             % [cm2]    L5 ACC --> PMC
Ce = (46*10^(-3)/236.4)/SurfArea;       % [uF/cm2] L5 ACC --> PMC

% SurfArea = 759.17 *10^(-8);             % [cm2]    L3 ACC --> AMY
% Ce = (28.1*10^(-3)/127.1)/SurfArea;     % [uF/cm2] L3 ACC --> AMY

% SurfArea = 921.96 *10^(-8);             % [cm2]    L5 ACC --> AMY 
% Ce = (37.2*10^(-3)/235.5)/SurfArea;     % [uF/cm2] L5 ACC --> AMY 
%-----------------------------------



%-----------------------------------
% synaptic reversal potentials: 
v_rev_e=0; v_rev_i=-68.2;    % min RMP for E - 5 & min RMP for I and I2 - 5
%-----------------------------------



%-----------------------------------
% time simulated (in ms):
t_final=500;  
%-----------------------------------



%-----------------------------------
% time step used in the midpoint method:
dt=0.02; 
%-----------------------------------



%-----------------------------------
% strength of synaptic connections
% g_hat_ie=1;                              % L3 ACC --> PMC 
% g_hat_i2e=g_hat_ie*6/12;                 % L3 ACC --> PMC 
% g_hat_ee=(g_hat_ie+g_hat_i2e)*0.086;     % L3 ACC --> PMC (E:I ratios based on total conductance event stats)
 
g_hat_ie=20/12;                          % L5 ACC --> PMC
g_hat_i2e=g_hat_ie*13/20;                % L5 ACC --> PMC
g_hat_ee=(g_hat_ie+g_hat_i2e)*0.111;     % L5 ACC --> PMC

% g_hat_ie=19/12;                          % L3 ACC --> AMY 
% g_hat_i2e=g_hat_ie*10/19;                % L3 ACC --> AMY 
% g_hat_ee=(g_hat_ie+g_hat_i2e)*0.086;     % L3 ACC --> AMY 

% g_hat_ie=12/12;                          % L5 ACC --> AMY 
% g_hat_i2e=g_hat_ie*31/12;              % L5 ACC --> AMY
% g_hat_ee=(g_hat_ie+g_hat_i2e)*0.181;     % L5 ACC --> AMY

g_hat_ei=g_hat_ee*10/4;           % Povysheva et al. 2008; Fig. 5 (monkey) 
g_hat_ei2=g_hat_ei/3;             % Bartos and Elgueta 2012; Fig. 1 (rat, hippocampus) 

g_hat_ii=g_hat_ei*5/4.5;          % Rotaru et al. 2015; Fig. 7 (monkey); https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4359991/

g_hat_ii2=g_hat_ii*1.37/2.16;     % Bezaire et al. 2016; table 1 (rodent)
g_hat_i2i2=g_hat_ii/3;            % cck->cck stronger than pv->pv; from Bezaire et al. 2016, table 1 (also some info in appendix 1 - tables 25 & 51: 7 times)

g_hat_i2i=1.75*g_hat_ii2;         % Karson et al. 2009 (rat); table 2: 7:4 ratio, 7/4=1.75
%-----------------------------------



%-----------------------------------
% external drive to the E-cells:

% deterministic drive:
Amp=21;
I_e=@(t) Amp*ones(num_e,1); 
I_e_a=@(t) Amp*ones(num_e_a,1); 

% maximum conductance, decay time, and frequency of Poisson train of excitatory input pulses:
g_stoch_e=0; tau_d_stoch_e=3; f_stoch_e=0;
%g_stoch_e=0.1; f_stoch_e=20; tau_d_stoch_e=3; 
%g_stoch_e=0.3; f_stoch_e=20; tau_d_stoch_e=3;
%-----------------------------------



%-----------------------------------
% external drive to the I-cells:

% deterministic drive:
I_i=@(t) zeros(num_i,1);

% maximum conductance, decay time, and frequency of Poisson train of excitatory input pulses:
g_stoch_i=0; tau_d_stoch_i=3; f_stoch_i=0;
%-----------------------------------



%-----------------------------------
% external drive to the I2-cells:

% deterministic drive:
I_i2=@(t) zeros(num_i2,1);

% maximum conductance, decay time, and frequency of Poisson train of excitatory input pulses:
g_stoch_i2=0; tau_d_stoch_i2=3; f_stoch_i2=0;
%-----------------------------------





