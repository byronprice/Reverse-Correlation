% SIMULATE DATA and test algorithm
numStimuli = 5000;train = round(numStimuli*0.8);
N = 40;

DIM = [N,N,numStimuli];beta = -3.3;

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
X = ifftn(tempFFT);
X = real(X);

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
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)+phi)...
    .*(k.*t).^n.*exp(-k.*t).*(1/gamma(n+1)-(k.*t).^2./(gamma(n+3)));

x = linspace(-20,20,N);y = linspace(-20,20,N);
t = linspace(200,0,1000/60);
[x,y,t] = meshgrid(x,y,t);
gabor = gaborFun(x,y,t,0,0,3,0.4,20,1,0.05,0,60*pi/180,0);
gaborEnergy = sum(sum(sum(gabor.*gabor)));

r = zeros(numStimuli-15,1);
filterOutput = zeros(numStimuli-15,1);
newS = zeros(numStimuli-15,N*N*16);
parfor ii=16:numStimuli
    Vid = single(S(:,:,ii-15:ii));
    filteredVid = Vid.*gabor;
    gaborOutput = sum(filteredVid(:));
    lambda = exp(gaborOutput.*gaborEnergy);
    r(ii-15) = poissrnd(lambda);
    filterOutput(ii-15) = gaborOutput;
    temp = single(unbiasedS(:,:,ii-15:ii));
    newS(ii-15,:) = temp(:);
end
clear unbiasedS S;
L = zeros(N*N*16,N*N*16,'single');

%operator = [0,-1,0;-1,4,-1;0,-1,0];
bigCount = 1;
for kk=1:16
    for jj=1:N
        for ii=1:N
            tempMat = zeros(N,N,16);
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
            if kk < 16
                tempMat(ii,jj,kk+1) = -1;
            end
            L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
        end
    end
end

% numBasis = N*N/4;
% xCen = linspace(-20,20,N/2);
% yCen = linspace(-20,20,N/2);
% 
% normalDeviance = zeros(8,1);
% poissonDeviance = zeros(8,1);
% bigCount = 1;
% for basisStd = [0.5,1,2,5,10,50]
%     basisFuns = zeros(N*N,numBasis);
%     count = 1;
%     for ii=1:length(xCen)
%         for jj=1:length(yCen)
%             gauss = gaussFilter(x,y,xCen(ii),yCen(jj),basisStd);
%             basisFuns(:,count) = gauss(:)./sum(gauss(:));
%             count = count+1;
%         end
%     end
%     
%     Design = newS*basisFuns;
%     [b,dev,stats] = glmfit(Design,r,'poisson','link','identity');
%     poissonDeviance(bigCount) = dev;
%     yhat = glmval(b,basisFuns,'identity');
%     figure();imagesc(reshape(yhat,[N,N]));title(sprintf('Poisson GLM: STD %3.1f',basisStd));
%     [b2,dev2,stats2] = glmfit(Design,r,'normal');
%     normalDeviance(bigCount) = dev2;bigCount = bigCount+1;
%     yhat = glmval(b2,basisFuns,'identity');
%     figure();imagesc(reshape(yhat,[N,N]));
%     title(sprintf('Normal GLM: STD %3.1f',basisStd));
% end
% figure();subplot(2,1,1);plot(normalDeviance);subplot(2,1,2);plot(poissonDeviance);

bigLambda = logspace(0,5,10);
RMS = zeros(length(bigLambda),1);
for ii=1:length(bigLambda)
    A = [newS(1:train,:);bigLambda(ii).*L];
    constraints = [r(1:train);zeros(N*N*16,1)];
    fhat = pinv(A)*constraints;
    RMS(ii) = norm(r(train+1:end)-newS(train+1:end,:)*fhat)./sqrt(N*N*16);
end
[~,bestMap] = min(RMS);

constraints = [r;zeros(N*N*16,1)];
A = [newS;bigLambda(bestMap).*L];
fhat = pinv(A)*constraints;

fhat = reshape(fhat,[N,N,16]);
minVal = min(fhat(:));maxVal = max(fhat(:));
for ii=1:16
   imagesc(fhat(:,:,ii));caxis([minVal maxVal]);pause(0.5);
end