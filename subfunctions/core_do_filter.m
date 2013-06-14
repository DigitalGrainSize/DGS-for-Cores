
% do_filter
% flattens current image
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

if length(sample)>1
    
    ButtonName = questdlg('Filter all images?','Filter all?', ...
        'Yes','No, just this one', 'Yes');
    
    if strcmp(ButtonName,'Yes')
        wh = waitbar(0,'Please wait, filtering images ...');
        
        for ii=1:length(sample)
            
            if ~sample(ii).filtered
                
                disp('Filtering image ...')
                
                % read data in if not already done so
                if isempty(sample(ii).data)
                    sample(ii).data=imread([image_path char(image_name(ii))]);
                    
                if numel(size(sample(ii).data))==3
                    sample(ii).data=double(0.299 * sample(ii).data(:,:,1) + 0.5870 * ...
                        sample(ii).data(:,:,2) + 0.114 * sample(ii).data(:,:,3));
                else
                    sample(ii).data=double(sample(ii).data);
                end

                    
                end
                
                [rows,cols] = size(sample(ii).data);
                sample(ix).orig_data=sample(ii).data;
                
                % filter parameters
                boost=sample(ii).filt1;
                CutOff=sample(ii).filt2;
                order=sample(ii).filt3;
                
                try
                    sample(ii).data= normalise(sample(ii).data);   % Rescale values 0-1 (and cast  to `double' if needed).
                    FFTlogIm = fft2(log(sample(ii).data+.01)); % Take FFT of log (with offset
                    % to avoid log of 0).
                    hb = highboostfilter([rows cols], CutOff, order, boost);
                    sample(ii).data = exp(real(ifft2(FFTlogIm.*hb)));  % Apply the filter, invert
                    % fft, and invert the log.
                    sample(ii).data=rescale(sample(ii).data,0,255);
                catch
                    disp('Error in filtering image. Continuing without filtering')
                end
                
                sample(ii).filtered = 1;
                
                for k=1:sample(ii).num_roi
                    sample(ii).roi{k}=sample(ii).data(min(sample(ii).roi_y{k}):...
                        max(sample(ii).roi_y{k}),...
                        min(sample(ii).roi_x{k}):...
                        max(sample(ii).roi_x{k}));
                end
                
                
            end
            waitbar(ii/length(sample),wh)
            
        end
        close(wh)
        
        for ii=1:length(sample)
            sample(ii).dist=[];
            sample(ii).percentiles=[];
            sample(ii).arith_moments=[];
            sample(ii).geom_moments=[];
        end
        %     set(findobj('tag','current_image'),'userdata',sample);
        
        
    else % no, just this image
        
        if ~sample(ix).filtered
            
            h = waitbar(0,'Please wait...');
            
            disp('Filtering image ...')
            [rows,cols] = size(sample(ix).data);
            sample(ix).orig_data=sample(ix).data;
            
            % filter parameters
            boost=sample(ix).filt1;
            CutOff=sample(ix).filt2;
            order=sample(ix).filt3;
            
            try
                sample(ix).data= normalise(sample(ix).data);   % Rescale values 0-1 (and cast  to `double' if needed).
                waitbar(.3,h)
                FFTlogIm = fft2(log(sample(ix).data+.01)); % Take FFT of log (with offset
                % to avoid log of 0).
                waitbar(.6,h)
                hb = highboostfilter([rows cols], CutOff, order, boost);
                waitbar(.8,h)
                sample(ix).data = exp(real(ifft2(FFTlogIm.*hb)));  % Apply the filter, invert
                % fft, and invert the log.
                waitbar(.9,h)
                sample(ix).data=rescale(sample(ix).data,0,255);
            catch
                disp('Error in filtering image. Continuing without filtering')
            end
            close(h)
            
            sample(ix).filtered = 1;
            
            for k=1:sample(ix).num_roi
                sample(ix).roi{k}=sample(ix).data(min(sample(ix).roi_y{k}):...
                    max(sample(ix).roi_y{k}),...
                    min(sample(ix).roi_x{k}):...
                    max(sample(ix).roi_x{k}));
            end
            
            if ~isempty(sample(ix).dist)
                sample(ix).dist=[];
                sample(ix).percentiles=[];
                sample(ix).arith_moments=[];
                sample(ix).geom_moments=[];
            end
            
        end
        
        
    end
    
else
    
    if ~sample(ix).filtered
        
        h = waitbar(0,'Please wait...');
        
        disp('Filtering image ...')
        [rows,cols] = size(sample(ix).data);
        sample(ix).orig_data=sample(ix).data;
        
        % filter parameters
        boost=sample(ix).filt1;
        CutOff=sample(ix).filt2;
        order=sample(ix).filt3;
        
        try
            sample(ix).data= normalise(sample(ix).data);   % Rescale values 0-1 (and cast  to `double' if needed).
            waitbar(.3,h)
            FFTlogIm = fft2(log(sample(ix).data+.01)); % Take FFT of log (with offset
            % to avoid log of 0).
            waitbar(.6,h)
            hb = highboostfilter([rows cols], CutOff, order, boost);
            waitbar(.8,h)
            sample(ix).data = exp(real(ifft2(FFTlogIm.*hb)));  % Apply the filter, invert
            % fft, and invert the log.
            waitbar(.9,h)
            sample(ix).data=rescale(sample(ix).data,0,255);
        catch
            disp('Error in filtering image. Continuing without filtering')
        end
        close(h)
        
        sample(ix).filtered = 1;
        
        for k=1:sample(ix).num_roi
            sample(ix).roi{k}=sample(ix).data(min(sample(ix).roi_y{k}):...
                max(sample(ix).roi_y{k}),...
                min(sample(ix).roi_x{k}):...
                max(sample(ix).roi_x{k}));
        end
        
        if ~isempty(sample(ix).dist)
            sample(ix).dist=[];
            sample(ix).percentiles=[];
            sample(ix).arith_moments=[];
            sample(ix).geom_moments=[];
        end
        
    end
    
    
    
end


set(findobj('tag','current_image'),'cdata',sample(ix).data);
set(findobj('tag','current_image'),'userdata',sample);
disp('... done!')

clear h k Nu Nv mag im auto nlags l centx centy


