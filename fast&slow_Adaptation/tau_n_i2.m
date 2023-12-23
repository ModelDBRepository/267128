function tau_n_i2=tau_n_i2(v);
alpha_n=0.032*(v+52)./(1-exp(-(v+52)/5));
beta_n=0.5*exp(-(v+57)/40);
tau_n_i2=1./(alpha_n+beta_n);
