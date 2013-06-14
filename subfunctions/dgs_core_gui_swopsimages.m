
% dgs_gui_swopsimages
% callback for main program, swops images
% 
% Written by Daniel Buscombe, various times in 2012 and 2013
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

% get current image
ix=get(findobj('tag','PickImage'),'value');
% and current image data
sample=get(findobj('tag','current_image'),'userdata');

% % get next image
if isempty(sample(ix).data)
    sample(ix).data=imread([image_path char(image_name(ix))]);
    
    if numel(size(sample(ix).data))==3
        sample(ix).data=double(0.299 * sample(ix).data(:,:,1) + 0.5870 * ...
            sample(ix).data(:,:,2) + 0.114 * sample(ix).data(:,:,3));
    else
        sample(ix).data=double(sample(ix).data);
    end
    
end

set(findobj('tag','current_image'),'cdata',sample(ix).data);

% set resolution bar to whatever the resolution actually is
set(findobj('tag','res'),'String',num2str(sample(ix).resolution));

% remove plotted roi lines from previous image
chx = get(ax,'Children');
if length(chx)>=2
    chx(end)=[];
    delete(chx)
end

% update title
% a=get(findobj('tag','im_axes1'),'title');
% set(get(findobj('tag','im_axes1'),'title'),'string',char(sample(ix).name));

h=findobj('tag','current_image');
set(h,'cdata',sample(ix).data); % make first image appear

[Nv,Nu,blank] = size(sample(ix).data);
set(h,'xdata',1:Nu); % scales and labels
set(h,'ydata',1:Nv);

set(ax,'ylim',[1,size(sample(ix).data,1)])
set(ax,'xlim',[1,size(sample(ix).data,2)])

% if navigating back, draw roi lines back on
if sample(ix).num_roi>0
    for k=1:sample(ix).num_roi
        sample(ix).roi_line{k} = line(sample(ix).roi_x{k},sample(ix).roi_y{k},'color','red','linewidth',2);
    end
end

% first set axes ticks to be increments of 500
set(ax,'ytick',linspace(1,size(sample(ix).data,1),2))
set(ax,'xtick',linspace(1,size(sample(ix).data,2),2))
% scale current x and y labels
set(ax,'xticklabels',num2str(get(ax,'xtick')'.*sample(ix).resolution))
set(ax,'yticklabels',num2str(get(ax,'ytick')'.*sample(ix).resolution))
% axis tight

if isfield(sample(ix),'roi_line')
    if iscell(sample(ix).roi_line)
    sample(ix).roi_line{1} = line(sample(ix).roi_x{1},...
        sample(ix).roi_y{1},'color','red','linewidth',5);
    end
end

axes(ax)


if ~isempty(sample(ix).dist)
    
    
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
    
else
    cla(ax2)
    title('')
end


set(findobj('tag','current_image'),'userdata',sample);

% clear chx k n Nu Nv mag im auto nlags l centx centy tmpimage h




