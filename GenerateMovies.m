% GenerateMovies.m
beta = -3;spaceExp = 2;timeExp = 2;
for jj=1:100
    numStimuli = 15*60*60;
    maxPix = 2350;
    minPix = 1425;
    
    screenPix_to_effPix = 25;
    effectivePixels = [maxPix/screenPix_to_effPix,minPix/screenPix_to_effPix];
    
    DIM = [effectivePixels(1),effectivePixels(2),numStimuli];
    
    u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
    v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
    t = [(0:floor(DIM(3)*2/2)) -(ceil(DIM(3)*2/2)-1:-1:1)]'/(2*DIM(3));
    [U,V,T] = meshgrid(u,v,t);
    S_f = single((U.^spaceExp+V.^spaceExp+T.^timeExp).^(beta/2));
    clear U V T u v t;
    S_f(S_f==inf) = 0;
    % S_f = S_f.^0.5;
    phi = rand([DIM(2),DIM(1),DIM(3)*2],'single');
    X = ifftn(S_f.^0.5.*(cos(2*pi*phi)+1i*sin(2*pi*phi)));
    X = real(X);
    
    % get unbiased movie
    tempFFT = fftn(X);S_f(S_f==0) = 1;
    unbiasedS = ifftn(tempFFT./S_f);
    unbiasedS = unbiasedS(:,:,1:numStimuli);
    
    X = X(:,:,1:numStimuli);
    desiredMin = 0;
    desiredMax = 255;
    Grey = 127;%desiredStd = 38;
    for ii=1:numStimuli
        temp = X(:,:,ii);
        currentMax = max(temp(:));
        currentMin = min(temp(:));
        temp = (desiredMax-desiredMin)./(currentMax-currentMin).*(temp-currentMax)+desiredMax;
        difference = mean(temp(:))-Grey;
        X(:,:,ii) = temp-difference;
        
        temp = unbiasedS(:,:,ii);
        currentMax = max(temp(:));
        currentMin = min(temp(:));
        temp = (desiredMax-desiredMin)./(currentMax-currentMin).*(temp-currentMax)+desiredMax;
        difference = mean(temp(:))-Grey;
        unbiasedS(:,:,ii) = temp-difference;
    end
    % currentMax = max(max(max(X)));
    % currentMin = min(min(min(X)));
    % X = (desiredMax-desiredMin)./(currentMax-currentMin).*(X-currentMax)+desiredMax;
    downSampleFactor = 1;
    S = uint8(X(:,:,1:downSampleFactor:end));
    numStimuli = numStimuli/downSampleFactor;
    clear X Y S_f;
    unbiasedS = uint8(unbiasedS);
    movie_FrameRate = 60;
    fileName = sprintf('15Min_PinkNoiseMovie%d.mat',jj);
    save(fileName,'S','screenPix_to_effPix','maxPix','minPix','beta','movie_FrameRate');
    
    fileName = sprintf('15Min_UnbiasedPinkNoiseMovie%d.mat',jj);
    save(fileName,'unbiasedS','screenPix_to_effPix','maxPix','minPix','beta','movie_FrameRate');
end