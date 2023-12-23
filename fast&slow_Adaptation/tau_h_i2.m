function tau_h_i2=tau_h_i2(v);
alpha_h=0.128*exp(-(v+50)/18);
beta_h=4./(1+exp(-(v+27)/5));
tau_h_i2=1./(alpha_h+beta_h);
