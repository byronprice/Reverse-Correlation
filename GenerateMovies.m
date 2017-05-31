% GenerateMovies.m

for jj=12:100
    beta = -3;spaceExp = 2;timeExp = 2;
    movie_FrameRate = 60; % hz
    movieTime_Seconds = 5*60;
    numStimuli = movieTime_Seconds*movie_FrameRate;
    maxPix = 2350;
    minPix = 1425;
    mmPerPixel = 0.2363;
    
    screenPix_to_effPix = 25;
    effectivePixels = [minPix/screenPix_to_effPix,maxPix/screenPix_to_effPix];
    
    % samplingSpace = screenPix_to_effPix*mmPerPixel/10;
    % nyquistSpace = (1/atand(samplingSpace/25))/2;
    % nyquistTime = movie_FrameRate/2;
    % u = (u./max(u)).*nyquistSpace;
    % v = (v./max(v)).*nyquistSpace;
    % t = (t./max(t)).*nyquistTime;
    
    DIM = [effectivePixels(1),effectivePixels(2),numStimuli];
    
    u = [(0:floor(DIM(1)*2/2)) -(ceil(DIM(1)*2/2)-1:-1:1)]'/(DIM(1)*2);
    v = [(0:floor(DIM(2)*2/2)) -(ceil(DIM(2)*2/2)-1:-1:1)]'/(2*DIM(2));
    t = [(0:floor(DIM(3)*2/2)) -(ceil(DIM(3)*2/2)-1:-1:1)]'/(2*DIM(3));
    [V,U,T] = meshgrid(v,u,t);
    S_f = (U.^spaceExp+V.^spaceExp+T.^timeExp).^(beta/2);
    clear U V T u v t;
    S_f(S_f==inf) = 0;
    % S_f = S_f.^0.5;
    phi = rand([DIM(1)*2,DIM(2)*2,DIM(3)*2]);
    tempFFT = S_f.^0.5.*(cos(2*pi*phi)+1i*sin(2*pi*phi));
    X = ifftn(tempFFT);
    X = real(X);
    
    % get unbiased movie
    S_f(S_f==0) = 1;
    
    desiredMin = 0;
    desiredMax = 255;
    Grey = 127;%desiredStd = 38;
    for ii=1:numStimuli*2
        temp = X(:,:,ii);
        currentMax = max(temp(:));
        currentMin = min(temp(:));
        temp = (desiredMax-desiredMin)./(currentMax-currentMin).*(temp-currentMax)+desiredMax;
        difference = mean(temp(:))-Grey;
        X(:,:,ii) = temp-difference;
    end
    
    S = uint8(X);
    unbiasedS = real(ifftn(fftn(double(S))./S_f));
    clear X Y S_f;
    S = S(1:DIM(1),1:DIM(2),1:numStimuli);
    unbiasedS = unbiasedS(1:DIM(1),1:DIM(2),1:numStimuli);
    
    fileName = sprintf('5Min_PinkNoiseMovie%d.mat',jj);
    save(fileName,'S','screenPix_to_effPix','maxPix','minPix','beta','movie_FrameRate',...
        'mmPerPixel','numStimuli','movieTime_Seconds','DIM');
    
    fileName = sprintf('5Min_UnbiasedPinkNoiseMovie%d.mat',jj);
    save(fileName,'unbiasedS','screenPix_to_effPix','maxPix','minPix','beta','movie_FrameRate',...
        'movieTime_Seconds','mmPerPixel','DIM');
    
    pause(2);clearvars -except jj;pause(2);
end
% from Dong 2001 ... Spatiotemporal Inseparability of Natural Images and
%   Visual Sensitivities

% power spectrum R is
%  f is spatial frequency, w is temporal frequency
%R(f,w) = (1/f^3.3)*P(w/f)  ... P(v) ~ 1/v^2

% vhat = 0.6 , r1 = 2 meters r2 = 40 meters n = 3.7, m = 2.3
%  v0 = 1.02
% R(f,w) = K*vhat/(2*f^(m-1)*w^2)*...
%   [(n-2)/(x+1)^(n-1)-(n-1)/(x+1)^(n-2)]
%  last part would be a subtraction with
%   x evaluated at w*r2/(f*v0) and
%   w*r1/(f*v0)
% u = (u./max(u)).*nyquistSpace;u(49:end) = u(49:end)-min(u);
% v = (v./max(v)).*nyquistSpace;v(30:end) = v(30:end)-min(v);
% t = (t./max(t)).*nyquistTime;t(1802:end) = t(1802:end)-min(t);
% vhat = 0.6;v0 = 1.02;r1 = 2;r2 = 40;n=3.7;m=2.3;
% % r1 and r2 in units of meters ... must convert to relevant units
% f = sqrt(U.^2+V.^2);w = T;f(f==0) = 1;w(w==0) = 1;
% x1 = w.*r2./(f.*v0);
% x2 = w.*r1./(f.*v0);
% S_f = (vhat./(2.*(f.^(m-1)).*(w.^2))).*...
%     (((n-2)./(x1+1).^(n-1)-(n-1)./(x2+1).^(n-2))-((n-2)./(x2+1).^(n-1)-(n-1)./(x2+1).^(n-2)));
