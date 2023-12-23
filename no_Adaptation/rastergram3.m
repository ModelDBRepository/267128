if num_spikes_i2>0, plot(t_i2_spikes,i_i2_spikes,'.','color',[0,0.8,0.8]); hold on; end
if num_spikes_i>0, plot(t_i_spikes,i_i_spikes+num_i2,'.b'); hold on; end
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i2+num_i,'.r');  hold on; end
if num_e>0 && num_i>0, plot([0,t_final],[num_i2,num_i2],'--k','Linewidth',1); end
if num_e>0 && num_i2>0, plot([0,t_final],[num_i2+num_i,num_i2+num_i],'--k','Linewidth',1); end

hold off;

set(gca,'Fontsize',12); 
axis([0,t_final,0,num_e+num_i+num_i2+1]); 
title('Pyramidal and CCK')
xlabel('time (ms)')

% if num_e<4 || num_i2<4 
%     set(gca,'Ytick',num_e+num_i2+num_i);
% else
%  	set(gca,'Ytick',[num_i2,num_e+num_i2+num_i]);
% end

set(gca,'Ytick',[num_i2,num_i2+num_i,num_e+num_i2+num_i]);

shg;  % show most recent graph window




