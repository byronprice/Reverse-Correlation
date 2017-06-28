function [PosteriorSamples,PosteriorMean,PosteriorInterval,Likelihood] =...
    BayesianFitSRFQuick(AnimalName,Date,NoiseType)
%BayesianFitSRFQuick.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatial receptive field for cells in V1
%  Images are displayed in image coordinates
%  MCMC
%    reduced model from BayesianFitSRF.m, which at present would require
%    3.5 seconds / iteration to run, with an estimated run-time of 64 days
%    to get 2000 posterior samples

% declare global variables
global numStimuli totalMillisecs ...
    unbiasedS numParameters ...
    movement X Y spikeTrain;
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','adfreq','allad','svStrobed','tsevs')
load(StimulusFileName)

% gaborFun = @(x,y,B,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
%     B.*exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
%     ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
%     .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)-sin(theta-pi/2).*(y-yc))-phi);
% 
% %nonLinFun = @(x,baseX,scale) exp((x-baseX)/scale);
% nonLinFun = @(x,base,slope,rise) rise./(1+exp(-(x-base).*slope));
%nonLinFun = @(x,slope,intercept) max(0,slope.*x+intercept);

% horzDegrees = atan((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
% vertDegrees = atan((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);

DIM = zeros(2,1);
if size(effectivePixels,2) == 1
    DIM(1) = sqrt(effectivePixels);
    DIM(2) = sqrt(effectivePixels);
elseif size(effectivePixels,2) == 2
    DIM(1) = effectivePixels(2);
    DIM(2) = effectivePixels(1);
end
    
xaxis = linspace(-round(screenPix_to_effPix*DIM(2)/2)+1,...
    round(screenPix_to_effPix*DIM(2)/2),DIM(2));
yaxis = linspace(round(3*screenPix_to_effPix*DIM(1)/4),...
    -round(screenPix_to_effPix*DIM(1)/4)+1,DIM(1));

[X,Y] = meshgrid(xaxis,yaxis);

% CREATE UNBIASED VERSION OF MOVIE BY DIVIDING OUT POWER SPECTRUM
%  USED TO GENERATE THE MOVIE

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
[V,U] = meshgrid(v,u);
S_f = (U.^2+V.^2).^(beta/2);

S_f(S_f==inf) = 1;
a = 0;b = 255;
unbiasedS = zeros(size(S));
for ii=1:numStimuli
    temp = reshape(real(ifftn(fftn(double(reshape(S(ii,:),[DIM(1),DIM(2)])))./S_f)),[DIM(1)*DIM(2),1])';
    currentMin = min(temp);currentMax = max(temp);
    temp = ((b-a).*(temp-currentMin))/(currentMax-currentMin)+a;
    unbiasedS(ii,:) = temp;
end

clear S S_f U V u v;

% REORGANIZE SPIKING DATA
temp = ~cellfun(@isempty,allts);
Chans = find(sum(temp,1));numChans = length(Chans);
totalUnits = sum(sum(temp));

temp = cell(totalUnits,1);
count = 1;
for ii=1:numChans
   for jj=1:nunits1
       if isempty(allts{jj,Chans(ii)}) == 0
           temp{count} = allts{jj,Chans(ii)};
           count = count+1;
       end
   end
end
allts = temp;

strobeStart = 33;
strobeData = tsevs{1,strobeStart};

if length(strobeData) ~= numStimuli && unique(svStrobed) == 0
    strobeData = strobeData(1:2:end);
elseif length(unique(svStrobed)) == 3
    strobeData = strobeData(svStrobed==1);
end

% GATHER LFP AND MOVEMENT DATA
nonEmptyAD = ~cellfun(@isempty,allad);
inds = find(nonEmptyAD==1);
LFP = cell(length(inds),1);
for ii=1:length(inds)
   LFP{ii} = allad{inds};
end

totalTime = length(LFP{1})./adfreq;

% if isempty(allad{49}) == 0
%     movement = allad{49};
% 
%     difference = length(LFP{1})-length(movement);
%     
%     if mod(difference,2) == 0
%         addOn = difference/2;
%         movement = [zeros(addOn-1,1);movement;zeros(addOn+1,1)];
%     else
%         addOn = floor(difference/2);
%         movement = [zeros(addOn,1);movement;zeros(addOn+1,1)];
%     end
%     tempMov = conv(abs(movement),ones(adfreq/2,1),'same');
%     tempMov = tempMov-mean(tempMov);
%     stdEst = 1.4826*mad(tempMov,1);
%     movement = single(tempMov>(3*stdEst));
%     clear tempMov stdEst;
% else
%     movement = zeros(length(LFP{1}),1);
% end

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
timeMultiplier = 1000;
totalMillisecs = round(totalTime*timeMultiplier);

stimTimes = round(strobeData.*timeMultiplier);
stimTimeIndex = zeros(totalMillisecs,1,'single');

if exist('flipIntervals','var')==1
    flipInterval = mean(flipIntervals);
end

for kk=1:numStimuli
   stimTimeIndex(stimTimes(kk)+0.05*timeMultiplier:stimTimes(kk)+0.15*timeMultiplier) = kk;
end

stimOn = stimTimeIndex>0;stimOff = 1-stimOn;

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = max(1,round(allts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

reducedSpikeCount = zeros(numStimuli,totalUnits);
for ii=1:totalUnits
    for jj=1:numStimuli
        reducedSpikeCount(jj,ii) = sum(pointProcessSpikes(stimTimeIndex==jj,ii));
    end
end

clear pointProcessSpikes stimTimeIndex;

% INITIALIZE NUMBER OF PARAMETERS
numFilters = 21;
gaborParams = 9;
nonLinParams = 3;
designParams = 1;
sParams = numFilters-1;
precisionParams = numFilters+designParams+nonLinParams;
numParameters = designParams+gaborParams*numFilters+nonLinParams+precisionParams+sParams;

% FUNCTIONS THAT WILL BE USED (BUT NOT CALLED ... MATLAB IS SLOW AF WITH
%  FUNCTION CALLS)
% logPoissonPDF = @(y,mu) y.*log(mu)-mu; % assumes exp(mu) ... true is sum[i=1:N] {y.*log(mu)-mu}
% logGammaPDF = @(x,a,b) a*log(b)-log(gamma(a))+(a-1).*log(x)-x.*b;
% logGammaPDFMeanDisperse = @(x,mu,phi) (x.*(-1/mu)-log(mu)).*phi+log(phi)*phi+(phi-1).*log(x)-log(gamma(phi)); 
%                         % is 1/dispersion ... similar to precision,
%                         % high phi == low variance 
% logVonMises = @(x,k,mu) k.*cos(x-mu);% -log(besseli(0,k)); ... last part does
%                                        % not depend on the data
% logNormPDF = @(x,mu,std) -0.5*log(std*std)-(x-mu).*(x-mu)./(2*std*std);

% INTIALIZE BOUNDS
Bounds = zeros(numParameters,2);

designBounds = repmat([-1000,1000],[designParams,1]); % two baselines, history, and movement

% determine bounds for spatial frequency
spatFreqBounds = [1/150,0.5]; %cpd, for the mouse visual system,
                 % the first value is 1/150 cpd or 1 cycle per 150 degrees,
                 % this is beyond what we can even display on our screen

mmpercycle = tan(((1./spatFreqBounds)./2).*pi./180).*(DistToScreen*10*2);
cppBounds = (1./mmpercycle)*conv_factor; % cycles per pixel Bounds

gaborBounds = zeros(gaborParams,2);
gaborBounds(1,:) = [-1e4,1e4]; % B for height of gabor
gaborBounds(2,:) = [-Inf,Inf]; % A orientation of Gabor
% gaborBounds(3,:) = [min(xaxis)-50,max(xaxis)+50]; % x center
% gaborBounds(4,:) = [min(yaxis)-50,max(yaxis)+50]; % y center
gaborBounds(3,:) = [-Inf,Inf];
gaborBounds(4,:) = [-Inf,Inf];
gaborBounds(5,:) = [1,2000]; % standard deviation x
gaborBounds(6,:) = [1,2000]; % standard deviation y
gaborBounds(7,:) = [0,Inf];%cppBounds; % spatial frequency
gaborBounds(8,:) = [-Inf,Inf]; %  orientation theta
gaborBounds(9,:) = [-Inf,Inf]; % phase shift phi

nonLinBounds = zeros(nonLinParams,2);
nonLinBounds(1,:) = [-Inf,Inf]; % sigmoid base
nonLinBounds(2,:) = [0,Inf]; % sigmoid slope
nonLinBounds(3,:) = [0,Inf]; % sigmoid rise

precisionBounds = repmat([-Inf,Inf],[precisionParams,1]); % for any precision parameter

% INITIALIZE INDEX VECTORS FOR EASY CALL TO THE CORRECT PARAMETERS
designVec = 1;
bVec = designVec+1:designParams+gaborParams;
cVec = bVec(end)+1:bVec(end)+gaborParams*(numFilters-1);
bcPrecisionVec = [bVec(1),cVec(1:gaborParams:end)];
nonLinVec = cVec(end)+1:cVec(end)+nonLinParams;
precisionVec = nonLinVec(end)+1:nonLinVec(end)+precisionParams;
sVec = precisionVec(end)+1:precisionVec(end)+sParams;

Bounds(designVec,:) = designBounds;
Bounds(bVec,:) = gaborBounds;
Bounds(cVec,:) = repmat(gaborBounds,[numFilters-1,1]);
Bounds(nonLinVec,:) = nonLinBounds;
Bounds(precisionVec,:) = precisionBounds;
Bounds(sVec,:) = repmat([-Inf,Inf],[numFilters-1,1]);

clear nonLinBounds gaborBounds historyMoveBaseBounds alphaBounds designBounds;

% INITIALIZE PRIORS

designPrior = zeros(length(designVec),2);

spatFreqPrior = [2.6667,0.03]; % gamma, see Stryker 2008 J Neuro, Fig. 6A
temp = max(gamrnd(spatFreqPrior(1),spatFreqPrior(2),[5000,1]),1/150);
temp = (1./(tan(((1./temp)./2).*pi./180).*(2*DistToScreen*10)))*conv_factor;
cppPrior = mle(temp,'distribution','gamma');

gaborPrior = zeros(gaborParams,2);
gaborPrior(1,:) = [0,0]; % normal, height of gabor, variance is a precision parameter
gaborPrior(2,:) = [0,1]; % von mises(mu,k), orientation of the exponential part of gabor
gaborPrior(3,:) = [0,1000*1000]; % normal(mean,variance), xc
gaborPrior(4,:) = [0,1000*1000]; % normal, yc
gaborPrior(5,:) = [20,15]; % gamma, std x
gaborPrior(6,:) = [19,13]; % gamma, std y
gaborPrior(7,:) = cppPrior; % gamma, spatial frequency converted to units of pixels
gaborPrior(8,:) = [0,0.5]; % von mises, gabor orientation (sinusoidal part)
gaborPrior(9,:) = [0,0.5]; % von mises, phase phi

nonLinPrior = [0,0;1,0;2,0]; % base slope rise ... slope and rise must be positive
precisionPrior = [1e-3,1e-3];

sPrior = [0.5,0.5];

variances = zeros(numParameters,1);
variances(designVec) = 0.1;

count = designVec(end);
for ii=1:numFilters
   variances(count+1) = 0.1;
   variances(count+2) = pi/4;
   variances(count+3) = gaborPrior(3,2)/100;
   variances(count+4) = gaborPrior(4,2)/100;
   variances(count+5) = gaborPrior(5,1)*(gaborPrior(5,2)^2);
   variances(count+6) = gaborPrior(6,1)*(gaborPrior(6,2)^2);
   variances(count+7) = gaborPrior(7,1)*(gaborPrior(7,2)^2);
   variances(count+8) = pi/4;
   variances(count+9)  = pi/4;
   count = count+gaborParams;
end

variances(nonLinVec) = ones(nonLinParams,1);
variances(precisionVec) = 10.*ones(precisionParams,1);
% variances(sVec) = (sPrior(1)*sPrior(2))/((sPrior(1)+sPrior(2))^2*(sPrior(1)+sPrior(2)+1));
variances(sVec) = var(-10.*log(1./betarnd(sPrior(1),sPrior(2),[1000,1])-1)).*ones(sParams,1);

numIter = 16e5;burnIn = 1e5;skipRate = 1000;
fullImSize = DIM(1)*DIM(2);optimalAccept = 0.234;

pi2 = pi/2;

PosteriorSamples = zeros(totalUnits,numParameters,length(burnIn+1:skipRate:numIter));
PosteriorMean = zeros(totalUnits,numParameters);
PosteriorInterval = zeros(totalUnits,numParameters,2);
for zz=1:totalUnits
        % ORGANIZE DATA
        spikeTrain = reducedSpikeCount(:,zz);
        designMatrix = ones(numStimuli,1);
        [b,~,~] = glmfit(designMatrix,spikeTrain,'poisson','constant','off');
       
        %MCMC intialization
        designPrior(designVec,:) = [b(1),0];
        
        numStarts = 5000;
        parameterVec = zeros(numParameters,numStarts);
        posteriorProb = zeros(numStarts,1);
        
        updateMu = zeros(numParameters,1);
        updateMu(designVec) = designPrior(designVec,1)+normrnd(0,0.1,[designParams,1]);
        count = designVec(end);
        for jj=1:numFilters
            updateMu(count+1) = normrnd(0,1e-4);
            updateMu(count+2) = normrnd(0,pi/2);
            updateMu(count+3) = normrnd(0,500);
            updateMu(count+4) = normrnd(0,500);
            updateMu(count+5) = gamrnd(gaborPrior(5,1),gaborPrior(5,2));
            updateMu(count+6) = gamrnd(gaborPrior(6,1),gaborPrior(6,2));
            updateMu(count+7) = gamrnd(gaborPrior(7,1),gaborPrior(7,2));
            updateMu(count+8) = normrnd(0,pi/2);
            updateMu(count+9) = normrnd(0,pi/2);
            count = count+gaborParams;
        end
        updateMu(nonLinVec) = [normrnd(0,1);gamrnd(4,1/4);gamrnd(8,1/4)];
        updateMu(precisionVec) = log(gamrnd(1,1,[precisionParams,1]));
        updateMu(sVec) = normrnd(0,30,[sParams,1]);
        
        updateMu = min(max(updateMu,Bounds(:,1)),Bounds(:,2));
        
        sigma = diag(variances);
        halfSigma = cholcov(sigma);
        identity = eye(numParameters);
        
        updateParam = 0.1;
        
        proposalMu = zeros(numParameters,1);
        for ii=1:numStarts
            parameterVec(designVec,ii) = designPrior(designVec,1);
            count = designVec(end);
            for jj=1:numFilters
                parameterVec(count+1,ii) = normrnd(0,1e-4);
                parameterVec(count+2,ii) = normrnd(0,pi2);
                parameterVec(count+3,ii) = normrnd(gaborPrior(3,1),sqrt(gaborPrior(3,2))/2);
                parameterVec(count+4,ii) = normrnd(gaborPrior(4,1),sqrt(gaborPrior(4,2))/2);
                parameterVec(count+5,ii) = gamrnd(gaborPrior(5,1),gaborPrior(5,2));
                parameterVec(count+6,ii) = gamrnd(gaborPrior(6,1),gaborPrior(6,2));
                parameterVec(count+7,ii) = gamrnd(gaborPrior(7,1),gaborPrior(7,2));
                parameterVec(count+8,ii) = normrnd(0,pi2);
                parameterVec(count+9,ii) = normrnd(0,pi2);
                count = count+gaborParams;
            end
            parameterVec(nonLinVec,ii) = [normrnd(0,1);gamrnd(4,1/4);gamrnd(8,1/4)];
            parameterVec(precisionVec,ii) = log(gamrnd(1,1,[precisionParams,1]));
%             parameterVec(sVec,ii) = log(betarnd(sPrior(1),sPrior(2),[sParams,1]));
            parameterVec(sVec,ii) = normrnd(0,30,[sParams,1]);
            
            parameterVec(:,ii) = min(max(parameterVec(:,ii),Bounds(:,1)),Bounds(:,2));

            
            % CALCULATE POSTERIOR
            W = zeros(fullImSize,numFilters-1);
            params = parameterVec(bVec,ii);
            temp = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
            b = temp(:);
            count = 1;
            for jj=1:numFilters-1
                params = parameterVec(cVec(count:count+8),ii);
                temp = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                    ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                    .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                W(:,jj) = temp(:);
                count = count+gaborParams;
            end
            nonLins = parameterVec(nonLinVec,ii);
            sVals = 1./(1+exp(-0.1.*parameterVec(sVec,ii)));
            S = diag(2.*binornd(1,sVals)-1);

            x = ones(numStimuli,1);
            for jj=1:numStimuli
                onScreenStim = unbiasedS(jj,:)';
                temp = 0.5*onScreenStim'*W*S*W'*onScreenStim+b'*onScreenStim;
                x(jj) = nonLins(3)/(1+exp(-(temp-nonLins(1))*nonLins(2)));
            end
            mu = exp(designMatrix*parameterVec(designVec,ii)).*x;
            loglikelihood = sum(spikeTrain.*log(mu)-mu);
            
            phi1 = exp(parameterVec(precisionVec(end-1),ii));
            phi2 = exp(parameterVec(precisionVec(end),ii));
            logprior = 0.5*parameterVec(precisionVec(1),ii)-...
                0.5*exp(parameterVec(precisionVec(1),ii))*...
                (parameterVec(designVec,ii)-designPrior(designVec,1))*...
                (parameterVec(designVec,ii)-designPrior(designVec,1))+...
                sum((numFilters*0.5)*parameterVec(precisionVec(2:1+numFilters),ii)-...
                0.5*exp(parameterVec(precisionVec(2:1+numFilters),ii)).*...
                (parameterVec(bcPrecisionVec,ii).*parameterVec(bcPrecisionVec,ii)))+...
                0.5*parameterVec(precisionVec(end-2),ii)-...
                0.5*exp(parameterVec(precisionVec(end-2),ii))*(parameterVec(nonLinVec(1),ii)*parameterVec(nonLinVec(1),ii))+...
                (parameterVec(nonLinVec(2),ii).*(-1/nonLinPrior(2,1))-log(nonLinPrior(2,1))).*phi1+...
                log(phi1)*phi1+(phi1-1).*log(parameterVec(nonLinVec(2),ii))-log(gamma(phi1))+...
                (parameterVec(nonLinVec(3),ii).*(-1/nonLinPrior(3,1))-log(nonLinPrior(3,1))).*phi2+...
                log(phi2)*phi2+(phi2-1).*log(parameterVec(nonLinVec(3),ii))-log(gamma(phi2))+...
                sum((precisionPrior(1)-1).*parameterVec(precisionVec,ii)-exp(parameterVec(precisionVec,ii)).*precisionPrior(2))+...
                sum((sPrior(1)-1).*log(sVals)+(sPrior(2)-1).*log(1-sVals));
            
            count = designVec(end);
            for jj=1:numFilters
                logprior = logprior+gaborPrior(2,2).*cos(parameterVec(count+2,ii))-...
                    0.5*log(gaborPrior(3,2))-...
                    (parameterVec(count+3,ii))*(parameterVec(count+3,ii))./...
                    (2*gaborPrior(3,2))-...
                    0.5*log(gaborPrior(4,2))-...
                   (parameterVec(count+4,ii))*(parameterVec(count+4,ii))./...
                    (2*gaborPrior(4,2))+...
                    (gaborPrior(5,1)-1).*log(parameterVec(count+5,ii))-parameterVec(count+5,ii)./gaborPrior(5,2)+...
                    (gaborPrior(6,1)-1).*log(parameterVec(count+6,ii))-parameterVec(count+6,ii)./gaborPrior(6,2)+...
                    (gaborPrior(7,1)-1).*log(parameterVec(count+7,ii))-parameterVec(count+7,ii)./gaborPrior(7,2)+...
                    gaborPrior(8,2).*cos(parameterVec(count+8,ii))+...
                    gaborPrior(9,2).*cos(parameterVec(count+9,ii));
                count = count+gaborParams;
            end

            posteriorProb(ii) = loglikelihood+logprior;
            % END CALCULATE POSTERIOR
            loglambda = log(2.38*2.38/numParameters);
            for kk=1:500

                pStar = parameterVec(:,ii)+mvnrnd(proposalMu,exp(loglambda).*sigma)';

                if sum(pStar<=Bounds(:,1)) == 0 && sum(pStar>=Bounds(:,2)) == 0
                    % calculate posterior
                    W = zeros(fullImSize,numFilters-1);
                    params = pStar(bVec);
                    temp = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                        ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                        .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                    b = temp(:);
                    count = 1;
                    for jj=1:numFilters-1
                        params = pStar(cVec(count:count+8));
                        temp = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                            ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                            .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                        W(:,jj) = temp(:);
                        count = count+gaborParams;
                    end
                    nonLins = pStar(nonLinVec);
                    sVals = 1./(1+exp(-0.1.*pStar(sVec)));
                    S = diag(2.*binornd(1,sVals)-1);
                    x = ones(numStimuli,1);
                    for jj=1:numStimuli
                        onScreenStim = unbiasedS(jj,:)';
                        temp = 0.5*onScreenStim'*W*S*W'*onScreenStim+b'*onScreenStim;
                        x(jj) = nonLins(3)/(1+exp(-(temp-nonLins(1))*nonLins(2)));
                    end
                    mu = exp(designMatrix*pStar(designVec)).*x;
                    loglikelihood = sum(spikeTrain.*log(mu)-mu);
                    
                    phi1 = exp(pStar(precisionVec(end-1)));
                    phi2 = exp(pStar(precisionVec(end)));
                    logprior = 0.5*pStar(precisionVec(1))-...
                        0.5*exp(pStar(precisionVec(1)))*...
                        (pStar(designVec)-designPrior(designVec,1))*...
                        (pStar(designVec)-designPrior(designVec,1))+...
                        sum((numFilters*0.5)*pStar(precisionVec(2:1+numFilters))-...
                        0.5*exp(pStar(precisionVec(2:1+numFilters))).*...
                        (pStar(bcPrecisionVec).*pStar(bcPrecisionVec)))+...
                        0.5*pStar(precisionVec(end-2))-...
                        0.5*exp(pStar(precisionVec(end-2)))*(pStar(nonLinVec(1))*pStar(nonLinVec(1)))+...
                        (pStar(nonLinVec(2)).*(-1/nonLinPrior(2,1))-log(nonLinPrior(2,1))).*phi1+...
                        log(phi1)*phi1+(phi1-1).*log(pStar(nonLinVec(2)))-log(gamma(phi1))+...
                        (pStar(nonLinVec(3)).*(-1/nonLinPrior(3,1))-log(nonLinPrior(3,1))).*phi2+...
                        log(phi2)*phi2+(phi2-1).*log(pStar(nonLinVec(3)))-log(gamma(phi2))+...
                        sum((precisionPrior(1)-1).*pStar(precisionVec)-exp(pStar(precisionVec)).*precisionPrior(2))+...
                        sum((sPrior(1)-1).*log(sVals)+(sPrior(2)-1).*log(1-sVals));
                    
                    count = designVec(end);
                    for jj=1:numFilters
                        logprior = logprior+gaborPrior(2,2).*cos(pStar(count+2))-...
                            0.5*log(gaborPrior(3,2))-...
                            (pStar(count+3))*(pStar(count+3))./...
                            (2*gaborPrior(3,2))-...
                            0.5*log(gaborPrior(4,2))-...
                            (pStar(count+4))*(pStar(count+4))./...
                            (2*gaborPrior(4,2))+...
                            (gaborPrior(5,1)-1).*log(pStar(count+5))-pStar(count+5)./gaborPrior(5,2)+...
                            (gaborPrior(6,1)-1).*log(pStar(count+6))-pStar(count+6)./gaborPrior(6,2)+...
                            (gaborPrior(7,1)-1).*log(pStar(count+7))-pStar(count+7)./gaborPrior(7,2)+...
                            gaborPrior(8,2).*cos(pStar(count+8))+...
                            gaborPrior(9,2).*cos(pStar(count+9));
                        count = count+gaborParams;
                    end
                    
                    logA = (loglikelihood+logprior)-posteriorProb(ii);
        
                    if log(rand) < logA
                        parameterVec(:,ii) = pStar;
                        posteriorProb(ii) = loglikelihood+logprior;
                    end
                    
                    meanSubtract = parameterVec(:,ii)-updateMu;
                    updateMu = updateMu+updateParam.*meanSubtract;
                    halfSigma = halfSigma+updateParam.*(triu((inv(halfSigma))*(halfSigma'*halfSigma+meanSubtract*...
                        meanSubtract')*((inv(halfSigma))')-identity)-halfSigma);
                    sigma = halfSigma'*halfSigma;
                    loglambda = loglambda+updateParam.*(exp(min(0,logA))-optimalAccept);
                    
                else
                    display('miss');
                    display(find(pStar<=Bounds(:,1)))
                    loglambda = loglambda+updateParam.*(-optimalAccept);
                end
                scatter(kk,posteriorProb(ii));hold on;pause(0.01);
            end
        end
        
        likelihoodXprev = GetLikelihood(parameterVec(:,1));
        
        for iter = 2:burnIn
            xStar = parameterVec(:,iter-1)+normrnd(0,sigma,[numParameters,1]);
            xStar([historyParams+4,historyParams+10,historyParams+11]) = ...
                mod(xStar([historyParams+4,historyParams+10,historyParams+11]),twopi);
            temp = xStar-Bounds(:,1);
            
            if sum(temp<0) > 0 % outside bounds
                parameterVec(:,iter) = parameterVec(:,iter-1);
                rejectRate = rejectRate+1;
            else
                likelihoodXstar = GetLikelihood(xStar);

                if log(rand) < (likelihoodXstar-likelihoodXprev)
                    parameterVec(:,iter) = xStar;
                    likelihoodXprev = likelihoodXstar;
                else
                    parameterVec(:,iter) = parameterVec(:,iter-1);
                end
            end
        end
        
        parameterVec = zeros(numParameters,(numIter-burnIn)/skipRate);
        posteriorProb = zeros((numIter-burnIn)/skipRate,1);
        
        
        
        PosteriorSamples(zz,:,:) = parameterVec;
        
        alpha = 0.05;
        
        numColumns = 4;
        numRows = floor(numParameters/numColumns)+mod(numParameters,numColumns);
       figure();
       for ii=1:numParameters
           subplot(numRows,numColumns,ii);histogram(PosteriorSamples(zz,:,:));
       end
end


end

