function m_i2_inf=m_i2_inf(v);
alpha_m=0.32*(v+54)./(1-exp(-(v+54)/4));
beta_m=0.28*(v+27)./(exp((v+27)/5)-1);
m_i2_inf=alpha_m./(alpha_m+beta_m);
