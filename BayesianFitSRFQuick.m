function [PosteriorSamples] =...
    BayesianFitSRFQuick(AnimalName,Date,NoiseType)
%BayesianFitSRFQuick.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatial receptive field for cells in V1
%  Images are displayed in image coordinates
%  MCMC
%    reduced model from BayesianFitSRF.m, which at present would require
%    3.5 seconds / iteration to run, with an estimated run-time of 64 days
%    to get 2000 posterior samples

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
totalUnits = sum(sum(temp))-numChans;

temp = cell(totalUnits,1);
count = 1;
for ii=1:numChans
   for jj=2:nunits1
       if isempty(allts{jj,Chans(ii)}) == 0
           temp{count} = allts{jj,Chans(ii)};
           count = count+1;
       end
   end
end
allts = temp;

strobeStart = 33;
strobeData = tsevs{1,strobeStart};

if length(strobeData) ~= numStimuli && length(unique(svStrobed)) == 1
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

if exist('flipIntervals','var')==1
    flipInterval = mean(flipIntervals);
end

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = max(1,round(allts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

reducedSpikeCount = zeros(totalUnits,numStimuli,1*timeMultiplier);
for ii=1:totalUnits
    for jj=1:numStimuli
        timeInds = stimTimes(jj):stimTimes(jj)+1*timeMultiplier-1;
        reducedSpikeCount(ii,jj,:) = pointProcessSpikes(timeInds,ii)';
    end
end

clear pointProcessSpikes;

% INITIALIZE NUMBER OF PARAMETERS
numFilters = 6; % 11 filters is 0.3 seconds/iteration, 6 is 0.14
gaborRepeatParams = 5;
gaborConstantParams = 4;
gaborParams = gaborRepeatParams+gaborConstantParams;
nonLinParams = 3;
designParams = 1;
sParams = numFilters-1;
precisionParams = numFilters+designParams+nonLinParams;
% timeParams = 2;
numParameters = designParams+gaborConstantParams+gaborRepeatParams*numFilters+...
    nonLinParams+precisionParams+sParams;%+timeParams;

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

% gaborBounds(1,:) = [min(xaxis)-50,max(xaxis)+50]; % x center
% gaborBounds(2,:) = [min(yaxis)-50,max(yaxis)+50]; % y center
gaborBounds(1,:) = [-Inf,Inf];% x center
gaborBounds(2,:) = [-Inf,Inf]; % y center
gaborBounds(3,:) = [-Inf,Inf]; % standard deviation x
gaborBounds(4,:) = [-Inf,Inf]; % standard deviation y
gaborBounds(5,:) = [-Inf,Inf]; % B for height of gabor
gaborBounds(6,:) = [-Inf,Inf]; % A orientation of Gabor
gaborBounds(7,:) = [-Inf,Inf];%cppBounds; % spatial frequency
gaborBounds(8,:) = [-Inf,Inf]; %  orientation theta
gaborBounds(9,:) = [-Inf,Inf]; % phase shift phi

nonLinBounds = zeros(nonLinParams,2);
nonLinBounds(1,:) = [-Inf,Inf]; % sigmoid base
nonLinBounds(2,:) = [-Inf,Inf]; % sigmoid slope
nonLinBounds(3,:) = [-Inf,Inf]; % sigmoid rise

precisionBounds = repmat([-Inf,Inf],[precisionParams,1]); % for any precision parameter

% timeBounds = [-Inf,log(500)]; % for choosing the time window

% could add parameters that chooses the window within which we collect
% spikes ... right now, i arbitrarily pick 50:150ms, but we could add
% two parameters, a ~ gamma(25,2) and b ~ gamma(225,0.6667) ... which
% basically give us normal distributions centered at 50 and 150 with
% variance of 100 ... the distribution for a is a bit skew

% INITIALIZE INDEX VECTORS FOR EASY CALL TO THE CORRECT PARAMETERS
designVec = 1;
constantVec = designVec+1:designVec+gaborConstantParams;
bVec = constantVec(end)+1:constantVec(end)+gaborRepeatParams;
cVec = bVec(end)+1:bVec(end)+gaborRepeatParams*(numFilters-1);
bcPrecisionVec = [bVec(1),cVec(1:gaborRepeatParams:end)];
nonLinVec = cVec(end)+1:cVec(end)+nonLinParams;
precisionVec = nonLinVec(end)+1:nonLinVec(end)+precisionParams;
sVec = precisionVec(end)+1:precisionVec(end)+sParams;
% timeVec = sVec(end)+1:sVec(end)+timeParams;

Bounds(designVec,:) = designBounds;
Bounds(constantVec,:) = gaborBounds(1:4,:);
Bounds(bVec,:) = gaborBounds(5:end,:);
Bounds(cVec,:) = repmat(gaborBounds(5:end,:),[numFilters-1,1]);
Bounds(nonLinVec,:) = nonLinBounds;
Bounds(precisionVec,:) = precisionBounds;
Bounds(sVec,:) = repmat([-Inf,Inf],[numFilters-1,1]);
% Bounds(timeVec,:) = repmat(timeBounds,[timeParams,1]);

clear nonLinBounds gaborBounds historyMoveBaseBounds designBounds timeBounds ...
    precisionBounds;

% INITIALIZE PRIORS



spatFreqPrior = [2.6667,0.03]; % gamma, see Stryker 2008 J Neuro, Fig. 6A
temp = max(gamrnd(spatFreqPrior(1),spatFreqPrior(2),[5000,1]),1/150);
temp = (1./(tan(((1./temp)./2).*pi./180).*(2*DistToScreen*10)))*conv_factor;
cppPrior = mle(temp,'distribution','gamma');clear temp;

gaborPrior = zeros(gaborParams,2);

gaborPrior(1,:) = [0,1000*1000]; % normal(mean,variance), xc
gaborPrior(2,:) = [0,1000*1000]; % normal, yc
gaborPrior(3,:) = [20,15]; % gamma, std x
gaborPrior(4,:) = [19,13]; % gamma, std y
gaborPrior(5,:) = [0,0]; % normal, height of gabor, variance is a precision parameter
gaborPrior(6,:) = [0,1]; % von mises(mu,k), orientation of the exponential part of gabor
gaborPrior(7,:) = cppPrior; % gamma, spatial frequency converted to units of pixels
gaborPrior(8,:) = [0,0.5]; % von mises, gabor orientation (sinusoidal part)
gaborPrior(9,:) = [0,0.5]; % von mises, phase phi

nonLinPrior = [0,0;1,0;2,0]; % base slope rise ... slope and rise must be positive
precisionPrior = [1e-3,1e-3];

sPrior = [0.5,0.5];

% timePrior = [100,0.5;800,0.25];

variances = zeros(numParameters,1);
variances(designVec) = 0.1;

count = designVec(end);
variances(count+1) = gaborPrior(1,2)/100;
variances(count+2) = gaborPrior(2,2)/100;
variances(count+3) = var(log(gamrnd(gaborPrior(3,1),gaborPrior(3,2),[1000,1])));
variances(count+4) = var(log(gamrnd(gaborPrior(4,1),gaborPrior(4,2),[1000,1])));
count = count+4;
for ii=1:numFilters
   variances(count+1) = 0.1;
   variances(count+2) = pi/4;
   variances(count+3) = var(log(gamrnd(gaborPrior(7,1),gaborPrior(7,2),[1000,1])));
   variances(count+4) = pi/4;
   variances(count+5)  = pi/4;
   count = count+gaborRepeatParams;
end

variances(nonLinVec) = ones(nonLinParams,1);
variances(precisionVec) = 10.*ones(precisionParams,1);
% variances(sVec) = (sPrior(1)*sPrior(2))/((sPrior(1)+sPrior(2))^2*(sPrior(1)+sPrior(2)+1));
variances(sVec) = var(-10.*log(1./betarnd(sPrior(1),sPrior(2),[1000,1])-1)).*ones(sParams,1);
% variances(timeVec) = var(log(normrnd(50,10,[1000,2])))';

numIter = 17e5;burnIn = 2e5;skipRate = 1000;
fullImSize = DIM(1)*DIM(2);optimalAccept = 0.234;

pi2 = pi/2;

myCluster = parcluster('local');

if getenv('ENVIRONMENT')
   myCluster.JobStorageLocation = getenv('TMPDIR'); 
end

parpool(myCluster,4);

PosteriorSamples = zeros(totalUnits,numParameters,length(burnIn+1:skipRate:numIter));
% PosteriorMean = zeros(totalUnits,numParameters);
% PosteriorInterval = zeros(totalUnits,numParameters,2);
parfor zz=1:totalUnits
        % ORGANIZE DATA
        spikeTrain = squeeze(reducedSpikeCount(zz,:,:));
        numStimuli = size(spikeTrain,1);
        
        spikeTrain = sum(spikeTrain(:,50:500),2);
        designMatrix = ones(numStimuli,1);
        [b,~,~] = glmfit(designMatrix,spikeTrain,'poisson','constant','off');
       
        designPrior = zeros(length(designVec),2);
        
        designPrior(designVec,:) = [b(1),0];
        
        %MCMC intialization
        numStarts = 10;
        parameterVec = zeros(numParameters,numStarts);
        posteriorProb = zeros(numStarts,1);
        
        updateMu = zeros(numParameters,1);
        updateMu(designVec) = designPrior(designVec,1)+normrnd(0,0.1,[designParams,1]);
        count = designVec(end);
        updateMu(count+1) = normrnd(0,sqrt(gaborPrior(1,2))/2);
        updateMu(count+2) = normrnd(0,sqrt(gaborPrior(2,2))/2);
        updateMu(count+3) = log(gamrnd(gaborPrior(3,1),gaborPrior(3,2)));
        updateMu(count+4) = log(gamrnd(gaborPrior(4,1),gaborPrior(4,2)));
        count = count+4;
        for jj=1:numFilters
            updateMu(count+1) = normrnd(0,1e-4);
            updateMu(count+2) = normrnd(0,pi2);
            updateMu(count+3) = log(gamrnd(gaborPrior(7,1),gaborPrior(7,2)));
            updateMu(count+4) = normrnd(0,pi2);
            updateMu(count+5) = normrnd(0,pi2);
            count = count+gaborRepeatParams;
        end
        updateMu(nonLinVec) = [normrnd(0,1);log(gamrnd(4,1/4));log(gamrnd(8,1/4))];
        updateMu(precisionVec) = log(gamrnd(1,1,[precisionParams,1]));
        updateMu(sVec) = normrnd(0,30,[sParams,1]);
%         updateMu(timeVec) = [log(normrnd(50,10));log(normrnd(200,10))];
        
        updateMu = min(max(updateMu,Bounds(:,1)),Bounds(:,2));
        
        sigma = diag(variances);
        halfSigma = cholcov(sigma);
        identity = eye(numParameters);
        pcaW = eye(numParameters); %pcaW = normrnd(0,1,numParameters);
        
        updateParam = 0.1;
        
        for ii=1:numStarts
            parameterVec(designVec,ii) = designPrior(designVec,1);
            count = designVec(end);
            parameterVec(count+1,ii) = normrnd(gaborPrior(1,1),sqrt(gaborPrior(1,2))/2);
            parameterVec(count+2,ii) = normrnd(gaborPrior(2,1),sqrt(gaborPrior(2,2))/2);
            parameterVec(count+3,ii) = log(gamrnd(gaborPrior(3,1),gaborPrior(3,2)));
            parameterVec(count+4,ii) = log(gamrnd(gaborPrior(4,1),gaborPrior(4,2)));
            count = count+4;
            for jj=1:numFilters
                parameterVec(count+1,ii) = normrnd(0,1e-4);
                parameterVec(count+2,ii) = normrnd(0,pi2);
                parameterVec(count+3,ii) = log(gamrnd(gaborPrior(7,1),gaborPrior(7,2)));
                parameterVec(count+4,ii) = normrnd(0,pi2);
                parameterVec(count+5,ii) = normrnd(0,pi2);
                count = count+gaborRepeatParams;
            end
            parameterVec(nonLinVec,ii) = [normrnd(0,1);log(gamrnd(4,1/4));log(gamrnd(8,1/4))];
            parameterVec(precisionVec,ii) = log(gamrnd(1,1,[precisionParams,1]));
            parameterVec(sVec,ii) = normrnd(0,30,[sParams,1]);
%             parameterVec(timeVec,ii) = [log(normrnd(50,10));log(normrnd(200,10))];
            
            parameterVec(:,ii) = min(max(parameterVec(:,ii),Bounds(:,1)),Bounds(:,2));
            
            % CALCULATE POSTERIOR
            W = zeros(fullImSize,numFilters-1);
            params = [parameterVec(bVec(1:2),ii);parameterVec(constantVec,ii);parameterVec(bVec(3:end),ii)];
            params(5:7) = exp(params(5:7));
            temp1 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
            b = temp1(:);
            count = 1;
            for jj=1:numFilters-1
                params = [parameterVec(cVec(count:count+1),ii);parameterVec(constantVec,ii);parameterVec(cVec(count+2:count+4),ii)];
                params(5:7) = exp(params(5:7));
                temp2 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                    ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                    .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                W(:,jj) = temp2(:);
                count = count+gaborRepeatParams;
            end
            nonLins = parameterVec(nonLinVec,ii);nonLins(2:3) = exp(nonLins(2:3));
            sVals = 1./(1+exp(-0.1.*parameterVec(sVec,ii)));
            S = repmat((2.*binornd(1,sVals)-1)',[numStimuli,1]);

%             x = ones(numStimuli,1);
%             for jj=1:numStimuli
%                 onScreenStim = unbiasedS(jj,:)';
%                 x(jj) = 0.5*(onScreenStim'*W)*S*(W'*onScreenStim);%+b'*onScreenStim;
% %                 x(jj) = nonLins(3)/(1+exp(-(temp-nonLins(1))*nonLins(2)));
%             end
            temp3 = unbiasedS*W;
            x = 0.5.*sum((temp3.*S).*temp3,2)+(b'*unbiasedS')';
            x = nonLins(3)./(1+exp(-(x-nonLins(1)).*nonLins(2)));

           % x = 0.5.*diag(unbiasedS*W*S*W'*unbiasedS')+(b'*unbiasedS')';
            mu = exp(designMatrix*parameterVec(designVec,ii)).*x;
%             windowTimes = ceil(exp(parameterVec(timeVec,ii)));
%             loglikelihood = sum(sum(spikeTrain(:,windowTimes(1):windowTimes(2)),2).*log(mu)-mu);
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
                0.5*exp(parameterVec(precisionVec(end-2),ii))*(nonLins(1)*nonLins(1))+...
                (nonLins(2).*(-1/nonLinPrior(2,1))-log(nonLinPrior(2,1))).*phi1+...
                log(phi1)*phi1+(phi1-1).*log(nonLins(2))-log(gamma(phi1))+...
                (nonLins(3).*(-1/nonLinPrior(3,1))-log(nonLinPrior(3,1))).*phi2+...
                log(phi2)*phi2+(phi2-1).*log(nonLins(3))-log(gamma(phi2))+...
                sum((precisionPrior(1)-1).*parameterVec(precisionVec,ii)-exp(parameterVec(precisionVec,ii)).*precisionPrior(2))+...
                sum((sPrior(1)-1).*log(sVals)+(sPrior(2)-1).*log(1-sVals));
%                 (timePrior(1,1)-1).*parameterVec(timeVec(1),ii)-exp(parameterVec(timeVec(1),ii))./timePrior(1,2)+...
%                 (timePrior(2,1)-1).*parameterVec(timeVec(2),ii)-exp(parameterVec(timeVec(2),ii))./timePrior(2,2);

            count = designVec(end);
            logprior = logprior-...
                    (parameterVec(count+1,ii))*(parameterVec(count+1,ii))./...
                    (2*gaborPrior(1,2))-...
                   (parameterVec(count+2,ii))*(parameterVec(count+2,ii))./...
                    (2*gaborPrior(2,2))+...
                    (gaborPrior(3,1)-1).*parameterVec(count+3,ii)-exp(parameterVec(count+3,ii))./gaborPrior(3,2)+...
                    (gaborPrior(4,1)-1).*parameterVec(count+4,ii)-exp(parameterVec(count+4,ii))./gaborPrior(4,2);
            count = count+4;
            for jj=1:numFilters
                logprior = logprior+gaborPrior(6,2).*cos(parameterVec(count+2,ii))+...
                    (gaborPrior(7,1)-1).*parameterVec(count+3,ii)-exp(parameterVec(count+3,ii))./gaborPrior(7,2)+...
                    gaborPrior(8,2).*cos(parameterVec(count+4,ii))+...
                    gaborPrior(9,2).*cos(parameterVec(count+5,ii));
                count = count+gaborRepeatParams;
            end

            posteriorProb(ii) = loglikelihood+logprior;
            % END CALCULATE POSTERIOR
            loglambda = log(2.38*2.38).*ones(numParameters,1);
            eigenvals = variances;
            for kk=1:2e4
                
                index = unidrnd(numParameters);%unidrnd(numParameters);
                lambda = loglambda(index);
                stdev = sqrt(exp(lambda).*eigenvals(index));
                pStar = parameterVec(:,ii)+pcaW(:,index)*normrnd(0,stdev);

                if sum(pStar<=Bounds(:,1)) == 0 && sum(pStar>=Bounds(:,2)) == 0
                    % calculate posterior
                    params = [pStar(bVec(1:2));pStar(constantVec);pStar(bVec(3:end))];
                    params(5:7) = exp(params(5:7));
                    temp1 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2)))...
                        .*((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2)))./(2*params(5)*params(5))-...
                        ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2)))...
                        .*((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2)))./(2*params(6)*params(6)))...
                        .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                    b = temp1(:);
                  
                    count = 1;
                    for jj=1:numFilters-1
                        params = [pStar(cVec(count:count+1));pStar(constantVec);pStar(cVec(count+2:count+4))];
                        params(5:7) = exp(params(5:7));
                        temp2 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2)))...
                            .*((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2)))./(2*params(5)*params(5))-...
                            ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2)))...
                            .*((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2)))./(2*params(6)*params(6)))...
                            .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                        W(:,jj) = temp2(:);
                        count = count+gaborRepeatParams;
                    end
                    
                    nonLins = pStar(nonLinVec);nonLins(2:3) = exp(nonLins(2:3));
                    sVals = 1./(1+exp(-0.1.*pStar(sVec)));
                    S = repmat((2.*binornd(1,sVals)-1)',[numStimuli,1]);
                    
                    %             x = ones(numStimuli,1);
                    %             for jj=1:numStimuli
                    %                 onScreenStim = unbiasedS(jj,:)';
                    %                 x(jj) = 0.5*(onScreenStim'*W)*S*(W'*onScreenStim);%+b'*onScreenStim;
                    % %                 x(jj) = nonLins(3)/(1+exp(-(temp-nonLins(1))*nonLins(2)));
                    %             end
                    temp3 = unbiasedS*W;
                    x = 0.5.*sum((temp3.*S).*temp3,2)+(b'*unbiasedS')';
                    x = nonLins(3)./(1+exp(-(x-nonLins(1)).*nonLins(2)));
                    
                    mu = exp(designMatrix*pStar(designVec)).*x;
%                     windowTimes = ceil(exp(pStar(timeVec)));
%                     loglikelihood = sum(sum(spikeTrain(:,windowTimes(1):windowTimes(2)),2).*log(mu)-mu);
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
                        0.5*exp(pStar(precisionVec(end-2)))*(nonLins(1)*nonLins(1))+...
                        (nonLins(2).*(-1/nonLinPrior(2,1))-log(nonLinPrior(2,1))).*phi1+...
                        log(phi1)*phi1+(phi1-1).*log(nonLins(2))-log(gamma(phi1))+...
                        (nonLins(3).*(-1/nonLinPrior(3,1))-log(nonLinPrior(3,1))).*phi2+...
                        log(phi2)*phi2+(phi2-1).*log(nonLins(3))-log(gamma(phi2))+...
                        sum((precisionPrior(1)-1).*pStar(precisionVec)-exp(pStar(precisionVec)).*precisionPrior(2))+...
                        sum((sPrior(1)-1).*log(sVals)+(sPrior(2)-1).*log(1-sVals));
%                         (timePrior(1,1)-1).*pStar(timeVec(1))-exp(pStar(timeVec(1)))./timePrior(1,2)+...
%                         (timePrior(2,1)-1).*pStar(timeVec(2))-exp(pStar(timeVec(2)))./timePrior(2,2);
                    
                    count = designVec(end);
                    logprior = logprior-...
                            (pStar(count+1))*(pStar(count+1))./...
                            (2*gaborPrior(1,2))-...
                            (pStar(count+2))*(pStar(count+2))./...
                            (2*gaborPrior(2,2))+...
                            (gaborPrior(3,1)-1).*pStar(count+3)-exp(pStar(count+3))./gaborPrior(3,2)+...
                            (gaborPrior(4,1)-1).*pStar(count+4)-exp(pStar(count+4))./gaborPrior(4,2);
                    count = count+4;
                    for jj=1:numFilters
                        logprior = logprior+gaborPrior(6,2).*cos(pStar(count+2))+...
                            (gaborPrior(7,1)-1).*pStar(count+3)-exp(pStar(count+3))./gaborPrior(7,2)+...
                            gaborPrior(8,2).*cos(pStar(count+4))+...
                            gaborPrior(9,2).*cos(pStar(count+5));
                        count = count+gaborRepeatParams;
                    end
                    logA = (loglikelihood+logprior)-posteriorProb(ii);
        
                    if log(rand) < logA && logA < Inf
                        parameterVec(:,ii) = pStar;
                        posteriorProb(ii) = loglikelihood+logprior;
                    end
                 
                    if mod(kk,250) == 0
                        meanSubtract = parameterVec(:,ii)-updateMu;
                        updateMu = updateMu+updateParam.*meanSubtract;
                        halfSigma = halfSigma+updateParam.*(triu((inv(halfSigma))*(halfSigma'*halfSigma+meanSubtract*...
                            meanSubtract')*((inv(halfSigma))')-identity)-halfSigma);
%                         sigma = halfSigma'*halfSigma;
%                         
%                         Z = inv(tril(pcaW'*pcaW)')'*pcaW';
%                         pcaW = sigma*Z'*inv(triu(Z*sigma*Z'));
%                         pcaW = normc(pcaW);
%                         eigenvals = diag(pcaW'*sigma*pcaW);
                    end
                    lambda = lambda+updateParam.*(exp(min(0,logA))-optimalAccept);
                    
                else
                    lambda = lambda+updateParam.*(-optimalAccept);
                end
%                 scatter(kk,posteriorProb(ii));%hold on;title(sprintf('Window: %d - %d',ceil(exp(parameterVec(timeVec(1),ii))),...
%                    % ceil(exp(parameterVec(timeVec(2),ii)))));
%                 axis([0 1e4 -2000 1e4]);
%                 hold on;
                loglambda(index) = lambda;
        
            end
        end
        indeces = find(posteriorProb~=Inf & posteriorProb ~= -Inf);
        parameterVec = parameterVec(:,indeces);
        posteriorProb = posteriorProb(indeces);
        [maxPost,ind] = max(posteriorProb);
        posteriorProb = maxPost;
        parameterVec = parameterVec(:,ind);
        maxPost
        
%         eigenvals = variances;
%         
%         pcaW = normrnd(0,1,numParameters);
%         pcaW = normc(pcaW);
        updateParam = logspace(-0.3,-2,burnIn);
        for iter=1:burnIn
                
                index = unidrnd(numParameters);
                lambda = loglambda(index);
                stdev = sqrt(exp(lambda).*eigenvals(index));
                pStar = parameterVec+pcaW(:,index)*normrnd(0,stdev);

                if sum(pStar<=Bounds(:,1)) == 0 && sum(pStar>=Bounds(:,2)) == 0
                    % calculate posterior
                    
                    params = [pStar(bVec(1:2));pStar(constantVec);pStar(bVec(3:end))];
                    params(5:7) = exp(params(5:7));
                    temp1 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                        ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                        .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                    b = temp1(:);
                    count = 1;
                    for jj=1:numFilters-1
                        params = [pStar(cVec(count:count+1));pStar(constantVec);pStar(cVec(count+2:count+4))];
                        params(5:7) = exp(params(5:7));
                        temp2 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                            ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                            .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                        W(:,jj) = temp2(:);
                        count = count+gaborRepeatParams;
                    end
                    nonLins = pStar(nonLinVec);nonLins(2:3) = exp(nonLins(2:3));
                    sVals = 1./(1+exp(-0.1.*pStar(sVec)));
                    S = repmat((2.*binornd(1,sVals)-1)',[numStimuli,1]);
                    
                    temp3 = unbiasedS*W;
                    x = 0.5.*sum((temp3.*S).*temp3,2)+(b'*unbiasedS')';
                    x = nonLins(3)./(1+exp(-(x-nonLins(1)).*nonLins(2)));
                    mu = exp(designMatrix*pStar(designVec)).*x;
%                     windowTimes = ceil(exp(pStar(timeVec)));
%                     loglikelihood = sum(sum(spikeTrain(:,windowTimes(1):windowTimes(2)),2).*log(mu)-mu);
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
                        0.5*exp(pStar(precisionVec(end-2)))*(nonLins(1)*nonLins(1))+...
                        (nonLins(2).*(-1/nonLinPrior(2,1))-log(nonLinPrior(2,1))).*phi1+...
                        log(phi1)*phi1+(phi1-1).*log(nonLins(2))-log(gamma(phi1))+...
                        (nonLins(3).*(-1/nonLinPrior(3,1))-log(nonLinPrior(3,1))).*phi2+...
                        log(phi2)*phi2+(phi2-1).*log(nonLins(3))-log(gamma(phi2))+...
                        sum((precisionPrior(1)-1).*pStar(precisionVec)-exp(pStar(precisionVec)).*precisionPrior(2))+...
                        sum((sPrior(1)-1).*log(sVals)+(sPrior(2)-1).*log(1-sVals));
%                         (timePrior(1,1)-1).*pStar(timeVec(1))-exp(pStar(timeVec(1)))./timePrior(1,2)+...
%                         (timePrior(2,1)-1).*pStar(timeVec(2))-exp(pStar(timeVec(2)))./timePrior(2,2);
                    
                    count = designVec(end);
                    logprior = logprior-...
                            (pStar(count+1))*(pStar(count+1))./...
                            (2*gaborPrior(1,2))-...
                            (pStar(count+2))*(pStar(count+2))./...
                            (2*gaborPrior(2,2))+...
                            (gaborPrior(3,1)-1).*pStar(count+3)-exp(pStar(count+3))./gaborPrior(3,2)+...
                            (gaborPrior(4,1)-1).*pStar(count+4)-exp(pStar(count+4))./gaborPrior(4,2);
                    count = count+4;
                    for jj=1:numFilters
                        logprior = logprior+gaborPrior(6,2).*cos(pStar(count+2))+...
                            (gaborPrior(7,1)-1).*pStar(count+3)-exp(pStar(count+3))./gaborPrior(7,2)+...
                            gaborPrior(8,2).*cos(pStar(count+4))+...
                            gaborPrior(9,2).*cos(pStar(count+5));
                        count = count+gaborRepeatParams;
                    end
                    
                    logA = (loglikelihood+logprior)-posteriorProb;
        
                    if log(rand) < logA && logA < Inf
                        parameterVec = pStar;
                        posteriorProb = loglikelihood+logprior;
                    end
                    
                    if mod(iter,250) == 0 
                        meanSubtract = parameterVec-updateMu;
                        updateMu = updateMu+updateParam(iter).*meanSubtract;
                        halfSigma = halfSigma+updateParam(iter).*(triu((inv(halfSigma))*(halfSigma'*halfSigma+meanSubtract*...
                            meanSubtract')*((inv(halfSigma))')-identity)-halfSigma);
%                         sigma = halfSigma'*halfSigma;
%                         
%                         Z = inv(tril(pcaW'*pcaW)')'*pcaW';
%                         pcaW = sigma*Z'*inv(triu(Z*sigma*Z'));
%                         pcaW = normc(pcaW);
%                         eigenvals = diag(pcaW'*sigma*pcaW);
                        lambda = lambda+updateParam(iter).*(exp(min(0,logA))-optimalAccept);
                    end
                    
                else
                    lambda = lambda+updateParam(iter).*(-optimalAccept);
                end
                loglambda(index) = lambda;
        
        end
        clear updateParam;
        
%         tempEigs = [];tempW = [];
%         for ii=1:numParameters
%            if eigenvals(ii) > 1e-8
%               tempEigs = [tempEigs;eigenvals(ii)];
%               tempW = [tempW,pcaW(:,ii)];
%            end
%         end
%         
%         eigenvals = tempEigs;
%         pcaW = tempW;clear tempEigs tempW;
%         effectiveParams = length(eigenvals);
        
        currentParams = parameterVec;currentProb = posteriorProb;
        parameterVec = zeros(numParameters,(numIter-burnIn)/skipRate);
        
        parameterVec(:,1) = currentParams;
        figure();plot(currentParams);
        currentProb
        
        sigma = halfSigma'*halfSigma;
        lambda = 2.38*2.38/numParameters;sigma = lambda.*sigma;
        proposalMu = zeros(numParameters,1);
        bigCount = 2;
        for iter=2:(numIter-burnIn)
                
%                 index = unidrnd(effectiveParams);
%                 lambda = loglambda(index);
%                 stdev = sqrt(exp(lambda).*eigenvals(index));
%                 pStar = currentParams+pcaW(:,index)*normrnd(0,stdev);
                pStar = currentParams+mvnrnd(proposalMu,sigma)';

                if sum(pStar<=Bounds(:,1)) == 0 && sum(pStar>=Bounds(:,2)) == 0
                    % calculate posterior
                    W = zeros(fullImSize,numFilters-1);
                    params = [pStar(bVec(1:2));pStar(constantVec);pStar(bVec(3:end))];
                    params(5:7) = exp(params(5:7));
                    temp1 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                        ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                        .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                    b = temp1(:);
                    count = 1;
                    for jj=1:numFilters-1
                        params = [pStar(cVec(count:count+1));pStar(constantVec);pStar(cVec(count+2:count+4))];
                        params(5:7) = exp(params(5:7));
                        temp2 = params(1).*exp(-((X-params(3)).*cos(params(2))-(Y-params(4)).*sin(params(2))).^2./(2*params(5)*params(5))-...
                            ((X-params(3)).*sin(params(2))+(Y-params(4)).*cos(params(2))).^2/(2*params(6)*params(6)))...
                            .*sin((2*pi.*params(7)).*(cos(params(8)-pi2).*(X-params(3))-sin(params(8)-pi2).*(Y-params(4)))-params(9));
                        W(:,jj) = temp2(:);
                        count = count+gaborRepeatParams;
                    end
                    nonLins = pStar(nonLinVec);nonLins(2:3) = exp(nonLins(2:3));
                    sVals = 1./(1+exp(-0.1.*pStar(sVec)));
                    S = repmat((2.*binornd(1,sVals)-1)',[numStimuli,1]);
                    
                    temp3 = unbiasedS*W;
                    x = 0.5.*sum((temp3.*S).*temp3,2)+(b'*unbiasedS')';
                    x = nonLins(3)./(1+exp(-(x-nonLins(1)).*nonLins(2)));
                    mu = exp(designMatrix*pStar(designVec)).*x;
%                     windowTimes = ceil(exp(pStar(timeVec)));
%                     loglikelihood = sum(sum(spikeTrain(:,windowTimes(1):windowTimes(2)),2).*log(mu)-mu);
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
                        0.5*exp(pStar(precisionVec(end-2)))*(nonLins(1)*nonLins(1))+...
                        (nonLins(2).*(-1/nonLinPrior(2,1))-log(nonLinPrior(2,1))).*phi1+...
                        log(phi1)*phi1+(phi1-1).*log(nonLins(2))-log(gamma(phi1))+...
                        (nonLins(3).*(-1/nonLinPrior(3,1))-log(nonLinPrior(3,1))).*phi2+...
                        log(phi2)*phi2+(phi2-1).*log(nonLins(3))-log(gamma(phi2))+...
                        sum((precisionPrior(1)-1).*pStar(precisionVec)-exp(pStar(precisionVec)).*precisionPrior(2))+...
                        sum((sPrior(1)-1).*log(sVals)+(sPrior(2)-1).*log(1-sVals));
%                         (timePrior(1,1)-1).*pStar(timeVec(1))-exp(pStar(timeVec(1)))./timePrior(1,2)+...
%                         (timePrior(2,1)-1).*pStar(timeVec(2))-exp(pStar(timeVec(2)))./timePrior(2,2);
                    
                    count = designVec(end);
                    logprior = logprior-...
                            (pStar(count+1))*(pStar(count+1))./...
                            (2*gaborPrior(1,2))-...
                            (pStar(count+2))*(pStar(count+2))./...
                            (2*gaborPrior(2,2))+...
                            (gaborPrior(3,1)-1).*pStar(count+3)-exp(pStar(count+3))./gaborPrior(3,2)+...
                            (gaborPrior(4,1)-1).*pStar(count+4)-exp(pStar(count+4))./gaborPrior(4,2);
                    count = count+4;
                    for jj=1:numFilters
                        logprior = logprior+gaborPrior(6,2).*cos(pStar(count+2))+...
                            (gaborPrior(7,1)-1).*pStar(count+3)-exp(pStar(count+3))./gaborPrior(7,2)+...
                            gaborPrior(8,2).*cos(pStar(count+4))+...
                            gaborPrior(9,2).*cos(pStar(count+5));
                        count = count+gaborRepeatParams;
                    end
                    
                    logA = (loglikelihood+logprior)-currentProb;
        
                    if log(rand) < logA && logA < Inf
                        currentParams = pStar;
                        currentProb = loglikelihood+logprior;
                    end
                end
                if mod(iter-1,skipRate) == 0
                    parameterVec(:,bigCount) = currentParams;
                    bigCount = bigCount+1;
                end
        end
        count = designVec(end);
        parameterVec(count+3,:) = exp(parameterVec(count+3,:));
        parameterVec(count+4,:) = exp(parameterVec(count+4,:));
        count = count+4;
        for jj=1:numFilters
            parameterVec(count+2,:) = mod(parameterVec(count+2,:)+pi,2*pi)-pi;
            parameterVec(count+3,:) = exp(parameterVec(count+3,:));
            parameterVec(count+4,:) = mod(parameterVec(count+4,:)+pi,2*pi)-pi;
            parameterVec(count+5,:) = mod(parameterVec(count+5,:)+pi,2*pi)-pi;
            count = count+gaborRepeatParams;
        end
        parameterVec(nonLinVec(2:3),:) = exp(parameterVec(nonLinVec(2:3),:));
        parameterVec(precisionVec,:) = exp(parameterVec(precisionVec,:));
        parameterVec(sVec,:) = 1./(1+exp(-0.1.*parameterVec(sVec,:)));
%         parameterVec(timeVec,:) = exp(parameterVec(timeVec,:));
        
        PosteriorSamples(zz,:,:) = parameterVec;
        
        
        numColumns = 4;
        numRows = floor(numParameters/numColumns)+mod(numParameters,numColumns);
        figure();
        for ii=1:numParameters
            subplot(numRows,numColumns,ii);histogram(PosteriorSamples(zz,:,:));
        end
end
delete(gcp);
fileName = strcat(EphysFileName(1:end-9),'-Results.mat');
save(fileName,'PosteriorSamples','unbiasedS','reducedSpikeCount','totalUnits','allts','stimTimes','xaxis','yaxis');
end

