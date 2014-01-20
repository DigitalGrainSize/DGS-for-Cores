
% un_flatten
% reverses any previous flattening on this image
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
%   This function is part of 'dgs-gui' software
%   This software is in the public domain because it contains materials that originally came 
%   from the United States Geological Survey, an agency of the United States Department of Interior. 
%   For more information, see the official USGS copyright policy at 
%   http://www.usgs.gov/visual-id/credit_usgs.html#copyright
%====================================

% get current image
% and current image data
sample=get(findobj('tag','current_image'),'userdata');

if sample(ix).flattened %|| ~sample(ix).filtered
    
    
    if isfield(sample(ix),'orig_data')
        sample(ix).data=sample(ix).orig_data;
        sample(ix).orig_data=[];
    end
    sample(ix).flattened = 0;
    
    set(findobj('tag','current_image'),'cdata',sample(ix).data);
    set(findobj('tag','current_image'),'userdata',sample);

    if isempty(sample(ix).data)
        dgs_core_gui_swopsimages
    end
    
    sample=get(findobj('tag','current_image'),'userdata');
    
    for k=1:sample(ix).num_roi
        sample(ix).roi{k}=sample(ix).data(min(sample(ix).roi_y{k}):...
            max(sample(ix).roi_y{k}),...
            min(sample(ix).roi_x{k}):...
            max(sample(ix).roi_x{k}));
    end
    
    set(findobj('tag','current_image'),'cdata',sample(ix).data);
    set(findobj('tag','current_image'),'userdata',sample);
    
    clear Nu Nv mag im auto nlags l centx centy
    
    uiwait(msgbox('filter removed',' '));
    
else
    uiwait(msgbox('image not been flattened',' '));
    
end
