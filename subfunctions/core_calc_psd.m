
% calc_psd
% calculates PSD for each ROI
% 
% Written by Daniel Buscombe, various times in 2012-2014
% while at
% School of Marine Science and Engineering, University of Plymouth, UK
% then
% Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 
% please contact:
% dbuscombe@usgs.gov
% for lastest code version please visit:
% https://github.com/dbuscombe-usgs
% see also (project blog):
% http://dbuscombe-usgs.github.com/
%====================================
%   This function is part of 'dgs-core-gui' software
%   This software is in the public domain because it contains materials that originally came 
%   from the United States Geological Survey, an agency of the United States Department of Interior. 
%   For more information, see the official USGS copyright policy at 
%   http://www.usgs.gov/visual-id/credit_usgs.html#copyright
%====================================

dofilt=0;
density=10;
start_size=3;

winsize = 50; % size of window in vertical direction

MotherWav='Morlet';
Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
    'Dj',1/8,... %8, ...    % this will do dj sub-octaves per octave
    'S0',start_size,...    % this says start at a scale of X pixels
    'J1',[],...
    'Mother',MotherWav);

if sample(ix).num_roi>0
    
    [vr,P,scale]=core_get_psd(sample(ix).roi{1},density,Args,ix,winsize);
    
    sample(ix).dist=P;
    sample(ix).scale=scale.*sample(ix).resolution;
    
    sample(ix).percentiles=zeros(size(P,2),9);
    sample(ix).geom_moments=zeros(size(P,2),4);
    sample(ix).arith_moments=zeros(size(P,2),4);
    
    for l=1:size(P,1)
        [sample(ix).percentiles(l,:),sample(ix).geom_moments(l,:),...
            sample(ix).arith_moments(l,:)]=gsdparams(P(l,:),sample(ix).scale);
        sample(ix).geom_moments(l,2) = 1000*2^-sample(ix).geom_moments(l,2);
    end
    
    sample(ix).locations=linspace(1,size(sample(ix).data,1),size(P,1));
    %[1:density:size(sample(ix).data,1)];
    
    h=findobj('Tag','plot_axes');
    axes(h)
    cla(ax2)
    
    pcolor(sample(ix).scale,sample(ix).locations,sample(ix).dist')
    shading flat
    
    if sample(ix).resolution==1
        xlabel('Size (Pixels)')
    else
        xlabel('Size (mm)')
    end
    ylabel('Row')
    axis tight
    set(gca,'ydir','normal')
    hold on
    plot(sample(ix).arith_moments(:,1),sample(ix).locations,...
        'r','linewidth',2)
    set(gca,'xscale','log'), grid off
    
    
    axes(ax)
    set(findobj('tag','current_image'),'userdata',sample);
    
    
else
    
    uiwait(msgbox('Create ROI first!','Warning','modal'));
    
end


