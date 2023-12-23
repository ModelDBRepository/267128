%%%% E (Pyramidal), initial values

%%% L3 ACC --> PMC
% Vrest_l3_pmc = -65.8-0.9;     % empirical resting membrane potentials from Maya's data
% Vrest_l3_pmc_er = 0.9*2;      % error in the resting membrane potentials
% v_e=Vrest_l3_pmc+rand(num_e,1)*Vrest_l3_pmc_er;  
% v_e_a=Vrest_l3_pmc+rand(num_e_a,1)*Vrest_l3_pmc_er;

%%% L5 ACC --> PMC
Vrest_l5_pmc = -59.8-1.2;
Vrest_l5_pmc_er = 1.2*2;
v_e=Vrest_l5_pmc+rand(num_e,1)*Vrest_l5_pmc_er; 
v_e_a=Vrest_l5_pmc+rand(num_e_a,1)*Vrest_l5_pmc_er;

%%% L3 ACC --> AMY
% Vrest_l3_amy = -64.1-1.1;
% Vrest_l3_amy_er = 1.1*2;
% v_e=Vrest_l3_amy+rand(num_e,1)*Vrest_l3_amy_er; 
% v_e_a=Vrest_l3_amy+rand(num_e_a,1)*Vrest_l3_amy_er;

%%% L5 ACC --> AMY
% Vrest_l5_amy = -60.8-1.0;
% Vrest_l5_amy_er = 1.0*2;
% v_e=Vrest_l5_amy+rand(num_e,1)*Vrest_l5_amy_er;
% v_e_a=Vrest_l5_amy+rand(num_e_a,1)*Vrest_l5_amy_er;

% slow adapters 
n_e=0*ones(num_e,1);
m_e=m_e_inf(v_e);
h_e=0*ones(num_e,1);
c_e=0*ones(num_e,1);
s_e=zeros(num_e,1);
% fast adapters
n_e_a=0*ones(num_e_a,1);
m_e_a=m_e_inf(v_e_a);
h_e_a=0*ones(num_e_a,1);
c_e_a=0*ones(num_e_a,1);
s_e_a=zeros(num_e_a,1);


%%%% I (PV), initial values
v_i=-63+rand(num_i,1)*14;     % from Zaitsev et al. 2009 (table3, LACs)
n_i=0*ones(num_i,1);
m_i=m_i_inf(v_i); 
h_i=0*ones(num_i,1); 
s_i=zeros(num_i,1);


%%%% I2 (CCK), initial values
v_i2=-63.2+rand(num_i2,1)*4;  % from Bezaire et al. 2016 (appendix1 - table1, Sara did: +/- 2)
n_i2=0*ones(num_i2,1);
m_i2=m_i2_inf(v_i2);         
h_i2=0*ones(num_i2,1); 
s_i2=zeros(num_i2,1);


%%%% initial values
s_stoch_e=zeros(num_e,1); 
s_stoch_e_a=zeros(num_e_a,1); 
s_stoch_i=zeros(num_i,1);
s_stoch_i2=zeros(num_i2,1);

num_spikes_e=0;
num_spikes_e_a=0;
num_spikes_i=0;
num_spikes_i2=0;
