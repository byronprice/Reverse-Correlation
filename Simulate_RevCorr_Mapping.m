% SIMULATE DATA and test algorithm
numStimuli = 5000;
N = 40;


% uniform random numbers seem to be preferrable to standard normals
newS = zeros(numStimuli,N*N);
uncorrectedS = zeros(numStimuli,N*N);
DIM = [N,N];beta = -2;

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
% Reproduce these frequencies along ever row
U = repmat(u,1,DIM(2)); 
% v is the set of frequencies along the second dimension.  For a square
% region it will be the transpose of u
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]/DIM(2);
% Reproduce these frequencies along ever column
V = repmat(v,DIM(1),1);

% [U,V] = meshgrid(u,v); U = U'; V = V';
% Generate the power spectrum
S_f = (U.^2 + V.^2).^(beta/2);

% Set any infinities to zero
S_f(S_f==inf) = 0;

tempSf = S_f; tempSf(tempSf==0) = 1;

gaborFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3)).*sin(2*pi*0.05.*(x.*cos(pi/4)+y.*sin(pi/4)));
gaussFilter = @(x,y,xCen,yCen,std) exp(-(x-xCen).^2./(2*std*std)-(y-yCen).^2./(2*std*std));

x = linspace(-20,20,N);y = linspace(-20,20,N);
[X,Y] = meshgrid(x,y);
gabor = gaborFilter(X,Y);
gauss = gaussFilter(X,Y,0,0,3);
r = zeros(numStimuli,1);
filterOutput = zeros(numStimuli,1);

desiredMax = 255;
desiredMin = 0;
parfor ii=1:numStimuli
    X = spatialPattern(DIM,beta);
    currentMax = max(X(:));
    currentMin = min(X(:));
    Y = (desiredMax)./(currentMax-currentMin).*(X-currentMax)+desiredMax;
    meanVal = mean(Y(:));difference = meanVal-127;
    Y = Y-difference;%Y = Y-127;
    
    gaussIm = 1;%conv2(Y,gauss,'same');
    contrastIm = Y./gaussIm;
    filteredIm = contrastIm.*gabor;
    gaborOutput = sum(filteredIm(:));
    lambda = exp(gaborOutput./(N*N));
    r(ii) = poissrnd(lambda);
    filterOutput(ii) = gaborOutput;
    
    % correct for biased power spectrum
    tempFFT = fft2(Y);
    corrected = ifft2(tempFFT./tempSf);
    newS(ii,:) = corrected(:);
    uncorrectedS(ii,:) = Y(:);
end

L = zeros(N*N,N*N,'single');

%operator = [0,-1,0;-1,4,-1;0,-1,0];
bigCount = 1;
for jj=1:N
    for ii=1:N
        tempMat = zeros(N,N);
        tempMat(ii,jj) = 4;
        if ii > 1
            tempMat(ii-1,jj) = -1;
        end
        if ii < N
            tempMat(ii+1,jj) = -1;
        end
        if jj > 1
            tempMat(ii,jj-1) = -1;
        end
        if jj < N
            tempMat(ii,jj+1) = -1;
        end
        L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
    end
end

numBasis = N*N/4;
xCen = linspace(-20,20,N/2);
yCen = linspace(-20,20,N/2);

basisFuns = zeros(N*N,numBasis);
basisStd = 0.85;
count = 1;
for ii=1:length(xCen)
    for jj=1:length(yCen)
        gauss = gaussFilter(X,Y,xCen(ii),yCen(jj),basisStd);
        basisFuns(:,count) = gauss(:)./sum(gauss(:));
        count = count+1;
    end
end

Design = newS*basisFuns;
[b,dev,stats] = glmfit(Design,r,'poisson');
yhat = glmval(b,basisFuns,'log');
figure();imagesc(reshape(yhat,[N,N]));title('Gaussian Basis Functions: Poisson GLM');
[b2,dev2,stats2] = glmfit(Design,r,'normal');
yhat = glmval(b2,basisFuns,'identity');
figure();imagesc(reshape(yhat,[N,N]));
title('Gaussian Basis Functions: Normal GLM');

bigLambda = [0,5e1,1e2,1e3,1e4,5e4,1e5,5e5];
RMS = zeros(length(bigLambda),1);
parfor ii=1:length(bigLambda)
    A = [newS;bigLambda(ii).*L];
    constraints = [r;zeros(N*N,1)];
    fhat = pinv(A)*constraints;
    RMS(ii) = (N*N)^(-0.5)*norm(r-newS*fhat);
%     [b,dev,~] = glmfit(A,constraints,'poisson');
%     figure();imagesc(reshape(b(2:end),[N,N]));
%     RMS(ii) = dev;
end
rmsDiff = diff(RMS);lambdaDiff = diff(bigLambda)';
deltaRMSdeltaLambda = rmsDiff./lambdaDiff;
[maxVal,index] = max(abs(deltaRMSdeltaLambda));
onepercent = 0.01*maxVal;
firstBelow = find(abs(deltaRMSdeltaLambda(index:end))<onepercent,1);
bestMap = index+firstBelow-1;

constraints = [r;zeros(N*N,1)];
A = [newS;bigLambda(bestMap).*L];
fhat = pinv(A)*constraints;

X = spatialPattern(DIM,beta);
currentMax = max(X(:));
currentMin = min(X(:));
Y = (desiredMax)./(currentMax-currentMin).*(X-currentMax)+desiredMax;
meanVal = mean(Y(:));difference = meanVal-127;
Y = Y-difference;
figure();subplot(3,1,1);imagesc(Y);
title('Example Image');
subplot(3,1,2);imagesc(gabor);title('Gabor Filter');
subplot(3,1,3);imagesc(reshape(fhat,[N,N]));
title(sprintf('Unbiased RF Estimate, Lambda: %3.0e',bigLambda(bestMap)));