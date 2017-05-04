% SIMULATE DATA and test algorithm
numStimuli = 10000;train = round(numStimuli*0.7);
N = 30;
DIM = [N,N,numStimuli];beta = -3;

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
t = [(0:floor(DIM(3)*2/2)) -(ceil(DIM(3)*2/2)-1:-1:1)]'/(2*DIM(3));
[U,V,T] = meshgrid(u,v,t);
S_f = single((U.^2+V.^2+T.^2).^(beta/2));
clear U V T u v t;
S_f(S_f==inf) = 0;
% S_f = S_f.^0.5;
phi = rand([DIM(2),DIM(1),DIM(3)*2],'single');
tempFFT = S_f.^0.5.*(cos(2*pi*phi)+1i*sin(2*pi*phi));
X = real(ifftn(tempFFT));

% get unbiased movie
S_f(S_f==0) = 1;

desiredMin = 0;
desiredMax = 255;
Grey = 127;%desiredStd = 38;
parfor ii=1:numStimuli*2
    temp = X(:,:,ii);
    currentMax = max(temp(:));
    currentMin = min(temp(:));
    temp = (desiredMax-desiredMin)./(currentMax-currentMin).*(temp-currentMax)+desiredMax;
    difference = mean(temp(:))-Grey;
    X(:,:,ii) = temp-difference;
end


S = uint8(X);
unbiasedS = real(ifftn(fftn(single(S))./S_f));
clear X Y S_f phi tempFFT;
unbiasedS = unbiasedS(:,:,1:numStimuli);
S = S(:,:,1:numStimuli);

gaborFun = @(x,y,t,xc,yc,sigma,k,n,gama,spatFreq,v,theta,phi) ...
    exp(-((x-xc).^2+gama.*gama.*(y-yc).^2)/(2*sigma*sigma))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)-phi)...
    .*(k.*t).^n.*exp(-k.*t).*(1/gamma(n+1)-(k.*t).^2./(gamma(n+3)));

x = linspace(-15,15,N);y = linspace(-15,15,N);
t = 150:-1000/50:0;
[x,y,t] = meshgrid(x,y,t);
gabor = gaborFun(x,y,t,0,0,3,0.4,20,1,0.05,0,60*pi/180,0);
gaborEnergy = sum(sum(sum(gabor.*gabor)));

timePoints = size(t,3);

r = zeros(numStimuli-(timePoints-1),1);
filterOutput = zeros(numStimuli-(timePoints-1),1);
newS = zeros(numStimuli-(timePoints-1),N*N*timePoints,'single');
for ii=timePoints:numStimuli
    Vid = single(S(:,:,ii-(timePoints-1):ii));
    filteredVid = Vid.*gabor;
    gaborOutput = sum(filteredVid(:)).*gaborEnergy;
    lambda = exp(25*gaborOutput);
    r(ii-(timePoints-1)) = poissrnd(lambda);
    filterOutput(ii-(timePoints-1)) = gaborOutput;
    temp = single(unbiasedS(:,:,ii-(timePoints-1):ii));
    newS(ii-(timePoints-1),:) = temp(:);
end
clear unbiasedS S;
L = zeros(N*N*timePoints,N*N*timePoints,'single');

%operator = [0,-1,0;-1,4,-1;0,-1,0];
bigCount = 1;
for kk=1:timePoints
    for jj=1:N
        for ii=1:N
            tempMat = zeros(N,N,timePoints);
            tempMat(ii,jj,kk) = 6;
            if ii > 1
                tempMat(ii-1,jj,kk) = -1;
            end
            if ii < N
                tempMat(ii+1,jj,kk) = -1;
            end
            if jj > 1
                tempMat(ii,jj-1,kk) = -1;
            end
            if jj < N
                tempMat(ii,jj+1,kk) = -1;
            end
            if kk > 1
                tempMat(ii,jj,kk-1) = -1;
            end
            if kk < timePoints
                tempMat(ii,jj,kk+1) = -1;
            end
            L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
        end
    end
end

bigLambda = logspace(1,4,10);
RMS = zeros(length(bigLambda),1);
tempF = zeros(length(bigLambda),N*N*timePoints);
for ii=1:length(bigLambda)
    A = [newS(1:train,:);bigLambda(ii).*L];
    constraints = [r(1:train);zeros(N*N*timePoints,1)];
    fhat = pinv(double(A))*double(constraints);
    tempF(ii,:) = fhat;
    RMS(ii) = norm(r(train+1:end)-newS(train+1:end,:)*fhat)./sqrt(N*N*timePoints);
    display(ii);
end
[~,bestMap] = min(RMS);

constraints = [double(r);zeros(N*N*timePoints,1)];
A = [double(newS);double(bigLambda(bestMap).*L)];
fhat = pinv(A)*constraints;

%fhat = tempF(bestMap,:);

fhat = reshape(fhat,[N,N,timePoints]);
minVal = min(fhat(:));maxVal = max(fhat(:));
min2 = min(gabor(:));max2 = max(gabor(:));
for ii=1:timePoints
   subplot(1,2,1);imagesc(fhat(:,:,ii));caxis([minVal maxVal]);
   subplot(1,2,2);imagesc(gabor(:,:,ii));caxis([min2 max2]);
   pause(0.5);
end