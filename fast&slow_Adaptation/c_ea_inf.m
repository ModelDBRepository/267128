function c_ea_inf=c_ea_inf(v);
% alpha_c=0.02*( v-20)./(1-exp(-(v-20)/9));
% beta_c=0.002*(v-20)./(exp((v-20)/9)-1);
% alpha_c=0.02*(v-10)./(1-exp(-(v-10)/9));
% beta_c=0.002*(v-10)./(exp((v-10)/9)-1);
alpha_c=0.0001*(v+30)./(1-exp(-(v+30)/9));
beta_c=0.0001*(v+30)./(exp((v+30)/9)-1);
c_ea_inf=alpha_c./(alpha_c+beta_c);
