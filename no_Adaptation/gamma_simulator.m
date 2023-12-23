
%%%%% Modified by Sara Ibañez from:
%%% Kopell N, Borgers C, Pervouchine D, Malerba P, Tort AB. 2010. Gamma and theta rhythms in biophysical models 
%%% of hippocampal circuits. Hippocampal Microcircuits: A Computational Modeller`s Resource Book New York, NY: Springer.


%%%%% Model for the 4 Pyramidal-PV-CCK networks without including adaptation.

% E (Pyramidal), I (PV) & I2 (CCK)

clear;
%close all
rng('default'); 

params; 

dt05=dt/2; m_steps=floor(t_final/dt+0.1); r_e=exp(-dt05/tau_d_stoch_e); r_i=exp(-dt05/tau_d_stoch_i);
r_i2=exp(-dt05/tau_d_stoch_i2);

% set synaptic conductances (by Sara: maximal conductance associated with the synapse)
g_ee=g_hat_ee*ones(num_e,num_e).*(sign(p_ee-rand(num_e,num_e))+1)/2/(p_ee*num_e);
g_ei=g_hat_ei*ones(num_e,num_i).*(sign(p_ei-rand(num_e,num_i))+1)/2/(p_ei*num_e);
g_ie=g_hat_ie*ones(num_i,num_e).*(sign(p_ie-rand(num_i,num_e))+1)/2/(p_ie*num_i);
g_ii=g_hat_ii*ones(num_i,num_i).*(sign(p_ii-rand(num_i,num_i))+1)/2/(p_ii*num_i);
g_ei2=g_hat_ei2*ones(num_e,num_i2).*(sign(p_ei2-rand(num_e,num_i2))+1)/2/(p_ei2*num_e);
g_ii2=g_hat_ii2*ones(num_i,num_i2).*(sign(p_ii2-rand(num_i,num_i2))+1)/2/(p_ii2*num_i);
g_i2i2=g_hat_i2i2*ones(num_i2,num_i2).*(sign(p_i2i2-rand(num_i2,num_i2))+1)/2/(p_i2i2*num_i2);
g_i2i=g_hat_i2i*ones(num_i2,num_i).*(sign(p_i2i-rand(num_i2,num_i))+1)/2/(p_i2i*num_i2);
g_i2e=g_hat_i2e*ones(num_i2,num_e).*(sign(p_i2e-rand(num_i2,num_e))+1)/2/(p_i2e*num_i2);

% initialize dynamic variables
initialize; 

% solve the system of Hodgkin-Huxley-like equations using the midpoint method
VE=zeros(num_e,m_steps);
%IstochE=zeros(num_e,m_steps);

for k=1:m_steps
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
    %%% E, inc
	v_e_inc=(0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +(g_ee'*s_e).*(v_rev_e-v_e)+(g_ie'*s_i).*(v_rev_i-v_e) ...
               +I_e(t_old)+g_stoch_e*s_stoch_e.*(v_rev_e-v_e) ...
               +(g_i2e'*s_i2).*(v_rev_i-v_e)) / Ce; 
           
	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    s_e_inc=0.5*(1+tanh(v_e/4)).*(1-s_e)./tau_r_e-s_e./tau_d_e; 
	%%% I, inc
    v_i_inc=0.1*(-65-v_i)+9*n_i.^4.*(-90-v_i)+35*m_i.^3.*h_i.*(55-v_i) ...
               +(g_ei'*s_e).*(v_rev_e-v_i)+(g_ii'*s_i).*(v_rev_i-v_i) ...
               +I_i(t_old)+g_stoch_i*s_stoch_i.*(v_rev_e-v_i) ...
               +(g_i2i'*s_i2).*(v_rev_i-v_i);
           
	n_i_inc=(n_i_inf(v_i)-n_i)./tau_n_i(v_i);
    h_i_inc=(h_i_inf(v_i)-h_i)./tau_h_i(v_i);
    s_i_inc=0.5*(1+tanh(v_i/4)).*(1-s_i)./tau_r_i-s_i./tau_d_i;
    %%% I2, inc
    v_i2_inc=0.1*(-65-v_i2)+9*n_i2.^4.*(-90-v_i2)+35*m_i2.^3.*h_i2.*(55-v_i2) ...        
               +(g_ei2'*s_e).*(v_rev_e-v_i2)+(g_ii2'*s_i).*(v_rev_i-v_i2) ...
               +(g_i2i2'*s_i2).*(v_rev_i-v_i2) ...
               +I_i2(t_old)+g_stoch_i2*s_stoch_i2.*(v_rev_e-v_i2); 
           
    n_i2_inc=(n_i2_inf(v_i2)-n_i2)./tau_n_i2(v_i2);    
    h_i2_inc=(h_i2_inf(v_i2)-h_i2)./tau_h_i2(v_i2);
    s_i2_inc=0.5*(1+tanh(v_i2/4)).*(1-s_i2)./tau_r_i2-s_i2./tau_d_i2;
    
    %%% E, tmp
	v_e_tmp=v_e+dt05*v_e_inc;
	n_e_tmp=n_e+dt05*n_e_inc;
	m_e_tmp=m_e_inf(v_e_tmp);
    h_e_tmp=h_e+dt05*h_e_inc;
    s_e_tmp=s_e+dt05*s_e_inc;    
	%%% I, tmp
    v_i_tmp=v_i+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
	m_i_tmp=m_i_inf(v_i_tmp);
    h_i_tmp=h_i+dt05*h_i_inc;
    s_i_tmp=s_i+dt05*s_i_inc;     
    %%% I2, tmp
    v_i2_tmp=v_i2+dt05*v_i2_inc;
	n_i2_tmp=n_i2+dt05*n_i2_inc;
	m_i2_tmp=m_i2_inf(v_i2_tmp);
    h_i2_tmp=h_i2+dt05*h_i2_inc;
    s_i2_tmp=s_i2+dt05*s_i2_inc;
    
    s_stoch_e=s_stoch_e*r_e;
    s_stoch_i=s_stoch_i*r_i;
    s_stoch_i2=s_stoch_i2*r_i2;
    
    %%% E, inc
	v_e_inc=(0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +(g_ee'*s_e_tmp).*(v_rev_e-v_e_tmp)+(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
               +I_e(t_mid)+g_stoch_e*s_stoch_e.*(v_rev_e-v_e_tmp) ...
               +(g_i2e'*s_i2_tmp).*(v_rev_i-v_e_tmp)) / Ce; 
           
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
    s_e_inc=0.5*(1+tanh(v_e_tmp/4)).*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e; 
	%%% I, inc
    v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp)+35*m_i_tmp.^3.*h_i_tmp.*(55-v_i_tmp) ...
               +(g_ei'*s_e_tmp).*(v_rev_e-v_i_tmp)+(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
               +I_i(t_mid)+g_stoch_i*s_stoch_i.*(v_rev_e-v_i_tmp) ...
               +(g_i2i'*s_i2_tmp).*(v_rev_i-v_i_tmp); 
           
 	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    s_i_inc=0.5*(1+tanh(v_i_tmp/4)).*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;   
    %%% I2, inc
    v_i2_inc=0.1*(-65-v_i2_tmp)+9*n_i2_tmp.^4.*(-90-v_i2_tmp)+35*m_i2_tmp.^3.*h_i2_tmp.*(55-v_i2_tmp) ...
               +(g_ei2'*s_e_tmp).*(v_rev_e-v_i2_tmp)+(g_ii2'*s_i_tmp).*(v_rev_i-v_i2_tmp) ...
               +(g_i2i2'*s_i2_tmp).*(v_rev_i-v_i2_tmp) ...
               +I_i2(t_mid)+g_stoch_i2*s_stoch_i2.*(v_rev_e-v_i2_tmp); 	
           
    n_i2_inc=(n_i2_inf(v_i2_tmp)-n_i2_tmp)./tau_n_i2(v_i2_tmp);
    h_i2_inc=(h_i2_inf(v_i2_tmp)-h_i2_tmp)./tau_h_i2(v_i2_tmp);
    s_i2_inc=0.5*(1+tanh(v_i2_tmp/4)).*(1-s_i2_tmp)./tau_r_i2-s_i2_tmp./tau_d_i2;
    
    %%% ve, vi
    v_e_old=v_e;
    v_i_old=v_i;
    v_i2_old=v_i2;
    
    %%% e
	v_e=v_e+dt*v_e_inc;
	n_e=n_e+dt*n_e_inc;
	m_e=m_e_inf(v_e);
    h_e=h_e+dt*h_e_inc;
    s_e=s_e+dt*s_e_inc;
    %%%% i
	v_i=v_i+dt*v_i_inc;
	n_i=n_i+dt*n_i_inc;
	m_i=m_i_inf(v_i);
    h_i=h_i+dt*h_i_inc;
    s_i=s_i+dt*s_i_inc;
    %%%% i2
	v_i2=v_i2+dt*v_i2_inc;
	n_i2=n_i2+dt*n_i2_inc;
	m_i2=m_i2_inf(v_i2);
    h_i2=h_i2+dt*h_i2_inc;
    s_i2=s_i2+dt*s_i2_inc;
    
    s_stoch_e=s_stoch_e*r_e; 
    s_stoch_i=s_stoch_i*r_i;
    s_stoch_i2=s_stoch_i2*r_i2;
    
    u_e=rand(num_e,1); 
    u_i=rand(num_i,1);
    u_i2=rand(num_i2,1);
    %a_rd=rand(20,1);

	s_stoch_e=s_stoch_e+max(sign(f_stoch_e*dt/1000-u_e),0).*(1-s_stoch_e);
	s_stoch_i=s_stoch_i+max(sign(f_stoch_i*dt/1000-u_i),0).*(1-s_stoch_i);
 	s_stoch_i2=s_stoch_i2+max(sign(f_stoch_i2*dt/1000-u_i2),0).*(1-s_stoch_i2);

	% determine which and how many e- and i-cells spiked in the current time step

    which_e=find(v_e_old<-5 & v_e >=-5); 
    which_i=find(v_i_old<-5 & v_i >=-5);
    which_i2=find(v_i2_old<-5 & v_i2 >=-5);
    
    l_e=length(which_e); 
    l_i=length(which_i);
    l_i2=length(which_i2);
    
    if l_e>0
        range=num_spikes_e+1:num_spikes_e+l_e;
		i_e_spikes(range)=which_e; 
        t_e_spikes(range)= ...
            (v_e(which_e)*(k-1)*dt-v_e_old(which_e)*k*dt)./(v_e(which_e)-v_e_old(which_e));
        num_spikes_e=num_spikes_e+l_e;
    end 

    if l_i>0
        range=num_spikes_i+1:num_spikes_i+l_i;
 		i_i_spikes(range)=which_i; 
        t_i_spikes(range)= ...
            (v_i(which_i)*(k-1)*dt-v_i_old(which_i)*k*dt)./(v_i(which_i)-v_i_old(which_i));
        num_spikes_i=num_spikes_i+l_i;
    end 
    
    if l_i2>0
        range=num_spikes_i2+1:num_spikes_i2+l_i2;
 		i_i2_spikes(range)=which_i2; 
        t_i2_spikes(range)= ...
            (v_i2(which_i2)*(k-1)*dt-v_i2_old(which_i2)*k*dt)./(v_i2(which_i2)-v_i2_old(which_i2));
        num_spikes_i2=num_spikes_i2+l_i2;
    end 

    VE(:,k)=v_e;
%    IstochE(:,k)=g_stoch_e*s_stoch_e.*(v_rev_e-v_e_tmp);

end


% plot the spike rastergram
figure;
rastergram3;

% plot action potentials
tt = 0:(t_final/m_steps):t_final-(t_final/m_steps);

figure;
plot(tt,VE(num_e/2,:))
title('Ve')
xlabel('time (ms)')
ylabel('membrane potential  (mV)')






