function n_i2_inf=n_i2_inf(v);
alpha_n=0.032*(v+52)./(1-exp(-(v+52)/5));
beta_n=0.5*exp(-(v+57)/40);
n_i2_inf=alpha_n./(alpha_n+beta_n);
