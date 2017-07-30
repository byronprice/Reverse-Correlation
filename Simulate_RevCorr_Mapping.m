% SIMULATE DATA and test algorithm
warning('off','all');
numStimuli = 5000;train = round(numStimuli*0.7);
N = 30;
DIM = [N,N,numStimuli];beta = -2;

u = [(0:floor(DIM(1)*2/2)) -(ceil(DIM(1)*2/2)-1:-1:1)]'/(2*DIM(1));
v = [(0:floor(DIM(2)*2/2)) -(ceil(DIM(2)*2/2)-1:-1:1)]'/(2*DIM(2));
t = [(0:floor(DIM(3)*2/2)) -(ceil(DIM(3)*2/2)-1:-1:1)]'/(2*DIM(3));
[V,U,T] = meshgrid(v,u,t);
S_f = (U.^2+V.^2+T.^2).^(beta/2);
clear U V T u v t;
S_f(S_f==inf) = 0;
% S_f = S_f.^0.5;
phi = rand([DIM(1)*2,DIM(2)*2,DIM(3)*2],'single');
tempFFT = S_f.*(cos(2*pi*phi)+1i*sin(2*pi*phi));
X = real(ifftn(tempFFT));

% get unbiased movie
S_f = 1./S_f;
S_f(S_f==inf) = 0;

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
unbiasedS = real(ifftn(fftn(double(S)).*S_f));
clear X Y S_f phi tempFFT;
unbiasedS = unbiasedS(1:DIM(1),1:DIM(2),1:numStimuli);
S = S(1:DIM(1),1:DIM(2),1:numStimuli);

a = 0;b = 255;
% for ii=1:numStimuli
%    temp = unbiasedS(:,:,ii); 
%    currentMin = min(temp(:));currentMax = max(temp(:));
%    temp = ((b-a).*(temp-currentMin))/(currentMax-currentMin)+a;
% %    subplot(2,1,1);imagesc(unbiasedS(:,:,ii));caxis([min(unbiasedS(:)) max(unbiasedS(:))]);
% %    subplot(2,1,2);imagesc(temp);caxis([a b]);pause(0.1);
%    unbiasedS(:,:,ii) = temp;
% end
temp = unbiasedS;
currentMin = min(temp(:));currentMax = max(temp(:));
unbiasedS = ((b-a).*(temp-currentMin))/(currentMax-currentMin)+a;
clear temp;

% S(S<60) = 0;S(S>=60 & S<196) = 127;S(S>=196) = 255;

gaborFun = @(x,y,t,k,n,v,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
    exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
    ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc)+v.*t)-phi)...
    .*(k.*t).^n.*exp(-k.*t).*(1/gamma(n+1)-(k.*t).^2./(gamma(n+3)));

nonLinFun = @(x,base,slope,rise,drop) rise./(1+exp(-(x-base)*slope))-drop;

x = linspace(-15,15,N);y = linspace(-15,15,N);
t = linspace(150,0,30);
[X,Y,T] = meshgrid(x,y,t);
gabor = gaborFun(X,Y,T,0.4,15,0,pi/6,0,0,4,3,0.1,60*pi/180,0);
gaborEnergy = sum(gabor(:).*gabor(:));
timePoints = size(T,3);

base = 50;slope = 0.1;rise = 3;drop = 0;

r = zeros(numStimuli-(timePoints-1),1);
filterOutput = zeros(numStimuli-(timePoints-1),1);
newS = zeros(numStimuli-(timePoints-1),N*N*timePoints);
for ii=timePoints:numStimuli
    Vid = double(S(:,:,ii-(timePoints-1):ii));
    filteredVid = Vid.*gabor;
    gaborOutput = sum(filteredVid(:));
    lambda = nonLinFun(gaborOutput,base,slope,rise,drop);
    r(ii-(timePoints-1)) = poissrnd(lambda);
    filterOutput(ii-(timePoints-1)) = gaborOutput;
    temp = unbiasedS(:,:,ii-(timePoints-1):ii);
%     temp = S(:,:,ii-(timePoints-1):ii);
    newS(ii-(timePoints-1),:) = temp(:);
end
clear unbiasedS S;

% GLM approach
% gaussFun = @(x,y,t,xc,yc,tc,stdxy,stdt) exp(-((x-xc).*(x-xc))./(2*stdxy*stdxy)-...
%     ((y-yc).*(y-yc))./(2*stdxy*stdxy)-((t-tc).*(t-tc))./(2*stdt*stdt));
% centerSpace = x(1:2:end);
% centerTime = t(1:2:end);
% numBasisSpace = length(centerSpace);
% numBasisTime = length(centerTime);
% basisSpaceStds = [1,5,6,10,15,20,25];
% basisTimeStds = [5,10,20,25,30,40,50,75,100];
% 
% numSpaceStds = length(basisSpaceStds);
% numTimeStds = length(basisTimeStds);
% 
% fullSize = N*N*timePoints;
% totalParams = numBasisSpace^2*numBasisTime;
% basisFuns = zeros(fullSize,totalParams);
% 
% heldOutDeviance = zeros(numSpaceStds,numTimeStds);
% for yy=1:numSpaceStds
%     for zz=1:numTimeStds
%         count = 1;
%         for ii=1:numBasisTime
%             for jj=1:numBasisSpace
%                 for kk=1:numBasisSpace
%                     temp = gaussFun(X,Y,T,centerSpace(jj),centerSpace(kk),...
%                         centerTime(ii),basisSpaceStds(yy),basisTimeStds(zz));
%                     basisFuns(:,count) = temp(:)./max(temp(:));
% %                     subplot(2,1,1);imagesc(temp(:,:,(ii-1)*2+1));
% %                     temp2 = temp(:,:,(ii-1)*2+1);
% %                     [~,idx] = max(temp2(:));
% %                     [row,col] = ind2sub(size(temp2),idx);
% %                     subplot(2,1,2);plot(squeeze(temp(row,col,:)));
% %                     pause(0.1);
%                     count = count+1;
%                 end
%             end
%         end
%         tic;
%         spikeTrain = r(1:train);
%         design = newS(1:train,:)*basisFuns;
%         [b,~,~] = glmfit(design,spikeTrain,'poisson');
%         
%         testDesign = [ones(length(r(train+1:end)),1),newS(train+1:end,:)*basisFuns];
%         testTrain = r(train+1:end);
%         mu = exp(testDesign*b);
%         temp = testTrain.*log(testTrain./mu)-(testTrain-mu);
%         temp(isnan(temp)) = mu(isnan(temp));
%         toc;
%         heldOutDeviance(yy,zz) = sum(temp);
%         display(heldOutDeviance(yy,zz));
%     end
% end
% 
% [~,idx] = min(heldOutDeviance(:));
% [row,col] = ind2sub(size(heldOutDeviance),idx);
% 
% bestSpaceStd = basisSpaceStds(row);
% bestTimeStd = basisTimeStds(col);
% 
% count = 1;
% for ii=1:numBasisTime
%     for jj=1:numBasisSpace
%         for kk=1:numBasisSpace
%             temp = gaussFun(X,Y,T,centerSpace(jj),centerSpace(kk),...
%                 centerTime(ii),bestSpaceStd,bestTimeStd);
%             basisFuns(:,count) = temp(:)./max(temp(:));
%             count = count+1;
%         end
%     end
% end
% 
% spikeTrain = r;
% design = newS*basisFuns;
% [b,~,~] = glmfit(design,spikeTrain,'poisson');
% 
% fhat = basisFuns*b(2:end);
% fhat = reshape(fhat,[N,N,timePoints]);
% minVal = min(fhat(:));maxVal = max(fhat(:));
% min2 = min(gabor(:));max2 = max(gabor(:));
% for ii=1:timePoints
%    subplot(1,2,1);imagesc(fhat(:,:,ii));caxis([minVal maxVal]);
%    subplot(1,2,2);imagesc(gabor(:,:,ii));caxis([min2 max2]);
%    pause(1);
% end

% regularized pseudoinverse approach
L = sparse(N*N*timePoints,N*N*timePoints);

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

bigLambda = logspace(0,4,10);
RMS = zeros(length(bigLambda),1);
% tempF = zeros(length(bigLambda),N*N*timePoints);
for ii=1:length(bigLambda)
    tic;
    A = [newS(1:train,:);bigLambda(ii).*L];
    constraints = [r(1:train);sparse(N*N*timePoints,1)];
    fhat = A\double(constraints);
%     tempF(ii,:) = fhat;
    RMS(ii) = norm(r(train+1:end)-newS(train+1:end,:)*fhat)./sqrt(N*N*timePoints);
    toc;
    display(RMS(ii));
end
[~,bestMap] = min(RMS);

constraints = [r;sparse(N*N*timePoints,1)];
A = [newS;bigLambda(bestMap).*L];
fhat = A\constraints;

%fhat = tempF(bestMap,:);

fhat = reshape(fhat,[N,N,timePoints]);
minVal = min(fhat(:));maxVal = max(fhat(:));
min2 = min(gabor(:));max2 = max(gabor(:));
for ii=1:timePoints
   subplot(1,2,1);imagesc(fhat(:,:,ii));caxis([minVal maxVal]);
   subplot(1,2,2);imagesc(gabor(:,:,ii));caxis([min2 max2]);
   pause(1);
end