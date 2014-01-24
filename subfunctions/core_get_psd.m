
function [vr,P,scale]=core_get_psd(himt,density,Args,ii, ywin)
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
% 
% for testing/debugging
% addpath(genpath(pwd))
% clear all;clc
% 
% dofilt=0;
% density=10;
% start_size=3;
% 
% MotherWav='Morlet';
% Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
%     'Dj',1/8,... %8, ...    % this will do dj sub-octaves per octave
%     'S0',start_size,...    % this says start at a scale of X pixels
%     'J1',[],...
%     'Mother',MotherWav);
% 
% ii=1;
% himt=double(imread('./images/313-M0027A-014H-01_scan.tiff_crop.tif'));
% 
% ywin = 50; 

xwin = size(himt,2);

% prc_overlapy = 50; % overlap in percent
% prc_overlapx = 99; % overlap in percent

% overlap in number of rows/columns
yshift = density; %round((xwin/100)*(100-prc_overlapy))
xshift = density; %round((xwin/100)*(100-prc_overlapx))

[Ny, Nx] = size(himt);

% Round xwin and ywin to nearest odd integer
xwin = fix(xwin); ywin = fix(ywin);
if rem(xwin,2)==0, xwin=xwin+1; end
if rem(ywin,2)==0, ywin=ywin+1; end

% Make x and y index vectors that define the window
xwin = (1:xwin)-(xwin+1)/2;
ywin = (1:ywin)-(ywin+1)/2;

% Identify window center locations
xshift = round(xshift);
yshift = round(yshift);

ic = 1:yshift:Ny;
jc = 1:xshift:Nx;

P1 = cell(length(ic),length(jc));

v = zeros(length(ic),length(jc));

h = waitbar(0,['Please wait... processing image ',num2str(ii)]);

for i = 1:length(ic)
    for j = 1:length(jc)
        
        x = himt(ic(i)+ywin((ic(i)+ywin)>0 & (ic(i)+ywin)<=Ny),...
            jc(j)+xwin((jc(j)+xwin)>0 & (jc(j)+xwin)<=Nx));
        
        x=x(1,:); x=x(:);
        tr=polyval(polyfit([1:length(x)]',x,1),[1:length(x)]');
        x=x-tr;
        
        cols = length(x);
        
        J1=fix((log(cols*1/Args.S0)/log(2))/Args.Dj);
        
        if Args.Pad == 1
            base2 = fix(log(cols)/log(2) + 0.4999);   % power of 2 nearest to N
            Y=zeros(1,cols+(2^(base2+1)-cols));
            Y(1:cols) = x - mean(x);
        else
            Y=x;
        end
        n = length(Y);
        
        k = [1:fix(n/2)];
        k = k.*((2.*pi)/(n*1));
        k = [0., k, -k(fix((n-1)/2):-1:1)];
        
        %....compute FFT of the (padded) time series
        f = fft(Y);    % [Eqn(3)]
        
        %....construct SCALE array & empty PERIOD & WAVE arrays
        scale = Args.S0*2.^((0:J1)*Args.Dj);
        
        wave = zeros(J1+1,n);  % define the wavelet array
        wave = wave + 1i*wave;  % make it complex
        % loop through all scales and compute transform
        for a1 = 1:J1+1
            [daughter,fourier_factor,coi,dofmin]=wave_bases(Args.Mother,k,scale(a1),-1);
            wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
        end
        
        wave = wave(:,1:cols);  % get rid of padding before returning
        
        sinv=1./(scale');
        wave=sinv(:,ones(1,cols)).*(abs(wave).^2);
        
        twave=zeros(size(wave));
        npad=2.^ceil(log2(cols));
        k = 1:fix(npad/2);
        k = k.*((2.*pi)/npad);
        k2=[0., k, -k(fix((npad-1)/2):-1:1)].^2;
        
        snorm=scale./1;
        for ii=1:size(wave,1)
            F=exp(-.5*(snorm(ii)^2)*k2);
            smooth=ifft(F.*fft(wave(ii,:),npad));
            twave(ii,:)=smooth(1:cols);
        end
        
        if isreal(wave)
            twave=real(twave);
        end
        
        P1{i,j}=var(twave,[],2);
        
        v(i,j) = sum(P1{i,j}./sum(P1{i,j}) .* scale');
        
        %keep P1 v himt Args scale h i j ic jc ii xwin ywin Ny Nx
        
    end
    
    waitbar(i/length(ic),h)
    
end

close(h)

vr = imresize(v,[Ny Nx]);
clear v

%imagesc(vr); axis image; colorbar

mindim = min(min(cell2mat(cellfun(@length,P1, 'UniformOutput',0))));

for i = 1:length(ic)
    for j = 1:length(jc)
        P1{i,j}=P1{i,j}(1:mindim);
        P1{i,j}=P1{i,j}./sum(P1{i,j});
    end
end

P=zeros(length(ic),length(scale));
for i = 1:length(ic)
    P(i,:) = mean(cell2mat({P1{i,:}}),2);
end

% P1 = cell2mat(P1);

