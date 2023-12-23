if num_spikes_i2>0, plot(t_i2_spikes,i_i2_spikes,'.','color',[0,0.8,0.8]); hold on; end
if num_spikes_i>0, plot(t_i_spikes,i_i_spikes+num_i2,'.b'); hold on; end
if num_spikes_e_a>0, plot(t_ea_spikes,i_ea_spikes+num_i2+num_i,'.m');  hold on; end
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i2+num_i+num_e_a,'.r');  hold on; end

plot([0,t_final],[num_i2,num_i2],'--k','Linewidth',1); 
plot([0,t_final],[num_i2+num_i,num_i2+num_i],'--k','Linewidth',1);
plot([0,t_final],[num_i2+num_i+num_e_a,num_i2+num_i+num_e_a],'--k','Linewidth',1);

hold off;

set(gca,'Fontsize',12); 
axis([0,t_final,0,num_e+num_e_a+num_i+num_i2+1]); 
title('Pyramidal and inhibitory cells')
xlabel('time (ms)')

% if num_e<4 || num_i2<4 
%     set(gca,'Ytick',num_e+num_i2+num_i);
% else
%  	set(gca,'Ytick',[num_i2,num_e+num_i2+num_i]);
% end

if num_e_a > 0
    set(gca,'Ytick',[num_i2,num_i2+num_i,num_e_a+num_i2+num_i,num_e+num_e_a+num_i2+num_i]);
else
    set(gca,'Ytick',[num_i2,num_i2+num_i,num_e+num_i2+num_i]);
end

shg;  % show most recent graph window




