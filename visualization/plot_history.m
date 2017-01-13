% plot global energy (non)conservation
% plots history.dat output files
% track.dat files contain global particle and field energy data 
% assumes history.dat file is in current working directory
%
% input:
% do_plot = 1 to make plots, do_plot = 0 to omit plotting
% (do_plot = 0 simply returns values read from the file)
% default: do_plot = 0
%
% do_save = 1 to save plots as .png, do_save = 0 to omit saving
% note: do_save takes precedence over do_plot, such that if do_save = 1
% and do_plot = 0, figures will be created and then immediately closed
% after saving.
% default: do_save = 1
%
% returns:
% t is a vector of times of length nt when particle trajectory was recorded
% px, py, pz are total particle momenta
% K is total particle energy 
% B2,E2 are magnetic and electric field energy 
% En is total energy (kinetic + EM)
%
% creates files (if do_save == 1): 
% history.png plot of momentum and energy vs time 

function [t,px,py,pz,K,B2,E2,En] = plot_history(do_plot,do_save)

if nargin < 2
    do_save = 1;
end
if nargin < 1
    do_plot = 0;
end

dname = './';
fname = [dname 'history.dat'];

% second line to deal with new header in track.dat files
% dat=dlmread(fname);
dat = dlmread(fname,'',1,0); 

t=dat(:,1);
px = dat(:,2);
py = dat(:,3);
pz = dat(:,4);
K = dat(:,5);
B2 = dat(:,6);
E2 = dat(:,7);
En = dat(:,8);

if (do_save || do_plot)
    do_close = ~do_plot;
     
    % parameters for plots of components
    f = figure;
    FS = 10;
    LW = 2; 
    legloc = 'best'; 
    nrows=2; 
    ncols=1;
    
    % plot 3 components of momentum
    subplot(nrows,ncols,1); 
    plot_components(t,[px py pz],LW);
    set(gca,'fontsize',FS); box on;
%     xlabel('Time');
    ylabel('Momentum');
    legend('p_x','p_y','p_z','Location',legloc);
    axis tight; box on ; 
    
    % plot E,B,K energies 
    subplot(nrows,ncols,2); 
    plot_components(t,[K E2 B2],LW); hold on ; 
        set(gca,'fontsize',FS); box on;
    plot(t,En,'k','linewidth',LW); hold off; 
    % use log scale while energy is blowing up
    set(gca,'yscale','log'); 
    xlabel('Time'); 
    ylabel('Energy'); 
    legend('Kinetic','Electric','Magnetic','Total','Location',legloc); 
    axis tight; box on; 
    save_and_close('history.png',f,do_save,do_close); 

end