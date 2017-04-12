% SIMULATE DATA and test algorithm
numStimuli = 5000;
N = 40;


% uniform random numbers seem to be preferrable to standard normals
newS = zeros(numStimuli,N*N);
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


gaborFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3)).*sin(2*pi*0.05.*(x.*cos(pi/4)+y.*sin(pi/4)));
gaussFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3));

x = linspace(-20,20,N);y = linspace(-20,20,N);
[X,Y] = meshgrid(x,y);
gabor = gaborFilter(X,Y);
gauss = gaussFilter(X,Y);
r = zeros(numStimuli,1);
filterOutput = zeros(numStimuli,1);
parfor ii=1:numStimuli
    X = spatialPattern(DIM,beta);
    X = X-min(min(X));
    X = (X./max(max(X))).*255;
    Y = X(:);
    meanVal = mean(Y);difference = meanVal-127;
    newS(ii,:) = Y-difference;
    gaussOutput = 1;%sum(sum(conv2(tempIm,gauss)));
    gaborOutput = sum(newS(ii,:)'.*gabor(:));
    lambda = exp((gaborOutput/gaussOutput)./(N*N));
    r(ii) = poissrnd(lambda);
    filterOutput(ii) = gaborOutput;
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
bigLambda = [0,5e1,1e2,1e3,1e4,5e4,1e5,1e6];
RMS = zeros(length(bigLambda),1);
parfor ii=1:length(bigLambda)
    A = [newS;bigLambda(ii).*L];
    constraints = [r;zeros(N*N,1)];
    fhat = pinv(A)*constraints;
    RMS(ii) = (N*N)^(-0.5)*norm(r-newS*fhat);
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

figure();subplot(4,1,1);imagesc(reshape(newS(10,:),[N,N]));
title('Example Image');
subplot(4,1,2);imagesc(gabor);title('Gabor Filter');
subplot(4,1,3);imagesc(reshape(fhat,[N,N]));
title(sprintf('Biased RF, Lambda: %3.0e',bigLambda(bestMap)));

S_f(S_f==0) = 1;
tempFFT = fft2(reshape(fhat,[N,N]));
correctedFhat = ifft2(numStimuli.*tempFFT./S_f);
correctedFhat = correctedFhat.*conj(correctedFhat);
subplot(4,1,4);imagesc(correctedFhat);
title(sprintf('Unbiased RF'));