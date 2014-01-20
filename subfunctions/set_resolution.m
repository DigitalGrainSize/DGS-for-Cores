
% set_resolution
% set resolution to whatever is in the text box and update axes data and
% plots accordingly
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
ButtonName = questdlg('Do this for all images?','Same resolution?', ...
    'Yes','No,', 'Yes');

if strcmp(ButtonName,'Yes')
    
    for ii=1:length(sample)
        
        sample(ii).resolution=str2double(get(findobj('tag','res'),'String'));

        if ~isempty(sample(ii).dist) && sample(ii).resolution==1
            orig_scale{ii}=sample(ii).scale;
        elseif ~isempty(sample(ii).dist) && sample(ii).resolution~=1
            orig_scale{ii}=sample(ii).scale.*1/sample(ii).resolution;
        end
                
        
        if sample(ii).resolution~=1
            % if there is a distribution but resolution is not 1
            % we need to replot ax2 (plot), replot ax3 (auto image)
            
            if ~isempty(sample(ii).dist)
                sample(ii).scale=sample(ii).scale.*sample(ix).resolution;
                
                sample(ii).percentiles=zeros(size(sample(ii).dist,2),9);
                sample(ii).geom_moments=zeros(size(sample(ii).dist,2),4);
                sample(ii).arith_moments=zeros(size(sample(ii).dist,2),4);
                
                for l=1:size(sample(ii).dist,2)
                    [sample(ii).percentiles(l,:),sample(ii).geom_moments(l,:),...
                        sample(ii).arith_moments(l,:)]=gsdparams(sample(ii).dist(:,l),sample(ii).scale);
                    sample(ii).geom_moments(l,2) = 1000*2^-sample(ii).geom_moments(l,2);
                end
                
                axes(ax2)
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
            end
            
            
        else %if sample(ii).resolution==1
            
            % if there is a distribution and resolution is 1
            % we need to replot ax2 (plot), replot ax3 (auto image)
            
            if ~isempty(sample(ii).dist)
                
                sample(ii).scale=sample(ii).scale.*sample(ix).resolution;
                
                sample(ii).percentiles=zeros(size(sample(ii).dist,2),9);
                sample(ii).geom_moments=zeros(size(sample(ii).dist,2),4);
                sample(ii).arith_moments=zeros(size(sample(ii).dist,2),4);
                
                for l=1:size(sample(ii).dist,2)
                    [sample(ii).percentiles(l,:),sample(ii).geom_moments(l,:),...
                        sample(ii).arith_moments(l,:)]=gsdparams(sample(ii).dist(:,l),sample(ii).scale);
                    sample(ii).geom_moments(l,2) = 1000*2^-sample(ii).geom_moments(l,2);
                end
                
                axes(ax2)
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
                
            end % empty
            
            
        end % res ==1
        
    end % ii
    
    if sample(ix).resolution~=1
        % turn labels to 'mm' rather than 'pixels'
        set(get(ax,'xlabel'),'string','mm')
        set(get(ax,'ylabel'),'string','mm')
        
        % first set axes ticks to be increments of 500
        set(ax,'ytick',(500:500:size(sample(ix).data,1)))
        set(ax,'xtick',(500:500:size(sample(ix).data,2)))
        % scale current x and y labels
        set(ax,'xticklabels',num2str(get(ax,'xtick')'.*sample(ix).resolution))
        set(ax,'yticklabels',num2str(get(ax,'ytick')'.*sample(ix).resolution))
        
    else
        % turn labels to 'Pixels' rather than 'mm'
        set(get(ax,'xlabel'),'string','Pixels')
        set(get(ax,'ylabel'),'string','Pixels')
        
        % first set axes ticks to be increments of 500
        set(ax,'ytick',(500:500:size(sample(ix).data,1)))
        set(ax,'xtick',(500:500:size(sample(ix).data,2)))
        % scale current x and y labels
        set(ax,'xticklabels',num2str(get(ax,'xtick')'.*sample(ix).resolution))
        set(ax,'yticklabels',num2str(get(ax,'ytick')'.*sample(ix).resolution))
       
    end
    
    
else % no just this one
    
    if ~isempty(sample(ix).dist) && sample(ix).resolution==1
        orig_scale{ix}=sample(ix).scale;
    elseif ~isempty(sample(ix).dist) && sample(ix).resolution~=1
        orig_scale{ix}=sample(ix).scale.*1/sample(ix).resolution;
    end
    
    sample(ix).resolution=str2double(get(findobj('tag','res'),'String'));
    
    
    if sample(ix).resolution~=1
        % if there is a distribution but resolution is not 1
        % we need to replot ax2 (plot), replot ax3 (auto image)
        
        if ~isempty(sample(ix).dist)
            
            sample(ix).scale=sample(ix).scale.*sample(ix).resolution;
            
            sample(ix).percentiles=zeros(size(sample(ix).dist,2),9);
            sample(ix).geom_moments=zeros(size(sample(ix).dist,2),4);
            sample(ix).arith_moments=zeros(size(sample(ix).dist,2),4);
            
            for l=1:size(sample(ix).dist,2)
                [sample(ix).percentiles(l,:),sample(ix).geom_moments(l,:),...
                    sample(ix).arith_moments(l,:)]=gsdparams(sample(ix).dist(:,l),sample(ix).scale);
                sample(ix).geom_moments(l,2) = 1000*2^-sample(ix).geom_moments(l,2);
            end
            
            axes(ax2)
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
            
        end
        
        % turn labels to 'mm' rather than 'pixels'
        set(get(ax,'xlabel'),'string','mm')
        set(get(ax,'ylabel'),'string','mm')
        
        % first set axes ticks to be increments of 500
        set(ax,'ytick',(500:500:size(sample(ix).data,1)))
        set(ax,'xtick',(500:500:size(sample(ix).data,2)))
        % scale current x and y labels
        set(ax,'xticklabels',num2str(get(ax,'xtick')'.*sample(ix).resolution))
        set(ax,'yticklabels',num2str(get(ax,'ytick')'.*sample(ix).resolution))
        
        
    else %if sample(ix).resolution==1
        
        % if there is a distribution and resolution is 1
        % we need to replot ax2 (plot), replot ax3 (auto image)
        
        if ~isempty(sample(ix).dist)
            sample(ix).scale=sample(ix).scale.*sample(ix).resolution;
            
            sample(ix).percentiles=zeros(size(sample(ix).dist,2),9);
            sample(ix).geom_moments=zeros(size(sample(ix).dist,2),4);
            sample(ix).arith_moments=zeros(size(sample(ix).dist,2),4);
            
            for l=1:size(sample(ix).dist,2)
                [sample(ix).percentiles(l,:),sample(ix).geom_moments(l,:),...
                    sample(ix).arith_moments(l,:)]=gsdparams(sample(ix).dist(:,l),sample(ix).scale);
                sample(ix).geom_moments(l,2) = 1000*2^-sample(ix).geom_moments(l,2);
            end
            
            axes(ax2)
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
            
        end
        
        % turn labels to 'Pixels' rather than 'mm'
        set(get(ax,'xlabel'),'string','Pixels')
        set(get(ax,'ylabel'),'string','Pixels')
        
        % first set axes ticks to be increments of 500
        set(ax,'ytick',(500:500:size(sample(ix).data,1)))
        set(ax,'xtick',(500:500:size(sample(ix).data,2)))
        % scale current x and y labels
        set(ax,'xticklabels',num2str(get(ax,'xtick')'.*sample(ix).resolution))
        set(ax,'yticklabels',num2str(get(ax,'ytick')'.*sample(ix).resolution))
        
        
    end
    
end

% finally make ax the current axes, and submit sample to userdata
axes(ax)
set(findobj('tag','current_image'),'userdata',sample);

