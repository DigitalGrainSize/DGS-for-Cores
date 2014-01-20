
function [P1,scale]=core_get_psd(himt,density,Args,ii)
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

% ii=1;
% density =10;
% himt=double(imread('./images/313-M0027A-014H-01_scan.tiff_crop.tif'));

blockSize = 500; % Rows in block.
prc_overlap = 50; % overlap in percent

% cut image into blocks with specified overlap

[rows blockSizeC] = size(himt);

% overlap in number of rows
overlap = round((blockSize/100)*(100-prc_overlap));

wholeBlockRows = floor(rows / blockSize);
blockVectorR = [blockSize * ones(1, wholeBlockRows), rem(rows, blockSize)];
blockVectorC = [blockSizeC * ones(1, 1), rem(blockSizeC, blockSizeC)];
ca1 = mat2cell(himt, blockVectorR, blockVectorC);

ca1(cellfun(@isempty,ca1))=[];
ca1 = ca1';

% shorten image by overlap pixels
himt = himt(overlap:end,:);
[rows blockSizeC] = size(himt);

wholeBlockRows = floor(rows / blockSize);
blockVectorR = [blockSize * ones(1, wholeBlockRows), rem(rows, blockSize)];
blockVectorC = [blockSizeC * ones(1, 1), rem(blockSizeC, blockSizeC)];
ca2 = mat2cell(himt, blockVectorR, blockVectorC);

ca2(cellfun(@isempty,ca2))=[];
ca2 = ca2';

% interleave the 2 cells
ca = cell(size(ca1,1)+size(ca2,1),1);
counter=1;
for k=1:2:size(ca,1)
    ca{k,1} = ca1{counter,1};
    counter = counter+1;
end

counter=1;
for k=2:2:size(ca,1)
    ca{k,1} = ca2{counter,1};
    counter = counter+1;
end

% process each block by row density 'density'
W1=cell(1,size(ca,1));

h = waitbar(0,['Please wait... processing image ',num2str(ii)]);

for kk=1:size(ca,1)
    
    himt = ca{kk,1};
    [rows,cols] = size(himt);
    
    W1{kk}=cell(1,size([1:density:rows],2));
    
    for j=1:density:rows
        
        x=himt(j,:); x=x(:);
        tr=polyval(polyfit([1:length(x)]',x,1),[1:length(x)]');
        x=x-tr;
        
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
        scale=scale./2;
        
        if isreal(wave)
            twave=real(twave);
        end
        
        W1{kk}{j}=var(twave,[],2);
        
        keep W1 j himt rows cols Args scale h density ca kk
        
        
    end
    waitbar(kk/size(ca,1),h)
    
    
end
close(h)


P1 = cell(size(ca,1),1);
for kk = 1:size(ca,1)
    W1{kk}=cell2mat(W1{kk});
    P1{kk} = ndnanfilter(W1{kk},@rectwin,[1 10]);
    for k=1:size(P1{kk},2)
        P1{kk}(:,k)=P1{kk}(:,k)./sum(P1{kk}(:,k));
    end
    %P1{kk} = mean(P1{kk},2);
end

P2 = [cell2mat(P1(1:end-2)'),cell2mat(P1(end-1)),cell2mat(P1(end))]';

% average by density rows
[rows blockSizeC] = size(P2);
blockSize = density;

wholeBlockRows = floor(rows / blockSize);
blockVectorR = [blockSize * ones(1, wholeBlockRows), rem(rows, blockSize)];
blockVectorC = [blockSizeC * ones(1, 1), rem(blockSizeC, blockSizeC)];
P1 = mat2cell(P2, blockVectorR, blockVectorC);

P1(cellfun(@isempty,P1))=[];
P1 = P1';
P1 = cellfun(@mean,P1, 'UniformOutput',0);
P1 = cell2mat(P1)';

