function [PosteriorSamples,PosteriorMean,PosteriorInterval,Likelihood] =...
    BayesianFitSRF(AnimalName,Date,NoiseType)
%BayesianFitSRF.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatial receptive field for cells in V1
%  MCMC

% declare global variables
global numStimuli totalMillisecs pointProcessStimTimes h ...
    unbiasedS numParameters designParams ...
    movement X Y gaborFun nonLinFun spikeTrain;
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','adfreq','allad','svStrobed','tsevs')
load(StimulusFileName)

gaborFun = @(x,y,B,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
    B.*exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
    ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)-sin(theta-pi/2).*(y-yc))-phi);

%nonLinFun = @(x,baseX,scale) exp((x-baseX)/scale);
nonLinFun = @(x,base,slope,rise) rise./(1+exp(-(x-base).*slope));
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
    
xaxis = linspace(-screenPix_to_effPix*DIM(2)/2,...
    screenPix_to_effPix*DIM(2)/2,DIM(2));
yaxis = linspace(-screenPix_to_effPix*DIM(1)/4,...
    3*screenPix_to_effPix*DIM(1)/4,DIM(1));

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

if length(strobeData) ~= numStimuli
    strobeData = strobeData(1:2:end);
end

% GATHER LFP AND MOVEMENT DATA
nonEmptyAD = ~cellfun(@isempty,allad);
inds = find(nonEmptyAD==1);
LFP = cell(length(inds),1);
for ii=1:length(inds)
   LFP{ii} = allad{inds};
end

totalTime = length(LFP{1})./adfreq;

if isempty(allad{49}) == 0
    movement = allad{49};

    difference = length(LFP{1})-length(movement);
    
    if mod(difference,2) == 0
        addOn = difference/2;
        movement = [zeros(addOn-1,1);movement;zeros(addOn+1,1)];
    else
        addOn = floor(difference/2);
        movement = [zeros(addOn,1);movement;zeros(addOn+1,1)];
    end
    tempMov = conv(abs(movement),ones(adfreq/2,1),'same');
    tempMov = tempMov-mean(tempMov);
    stdEst = 1.4826*mad(tempMov,1);
    movement = single(tempMov>(3*stdEst));
    clear tempMov stdEst;
else
    movement = zeros(length(LFP{1}),1);
end

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
timeMultiplier = 1000;
totalMillisecs = round(totalTime*timeMultiplier);

stimTimes = round(strobeData.*timeMultiplier);
pointProcessStimTimes = ones(totalMillisecs,numStimuli,'single');

if exist('flipIntervals','var')==1
    flipInterval = mean(flipIntervals);
end

for kk=1:numStimuli
   pointProcessStimTimes(stimTimes(kk)+0.05*timeMultiplier:stimTimes(kk)+0.15*timeMultiplier,kk) = 1;
end

temp = sum(pointProcessStimTimes,2);
stimOn = temp>0;stimOff = 1-stimOn;

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = max(1,round(allts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

% INITIALIZE NUMBER OF PARAMETERS
numFilters = 21;
gaborParams = 9;
nonLinParams = 3;
historyParams = 100;
moveParams = 1;
baseParams = 2;
precisionParams = numFilters+moveParams+baseParams+nonLinParams;
designParams = historyParams+moveParams+baseParams;
numParameters = designParams+gaborParams*numFilters+nonLinParams+precisionParams;

% FUNCTIONS THAT WILL CALLED
logPoissonPDF = @(y,mu) y.*log(mu)-mu; % assumes exp(mu) ... true is sum[i=1:N] {y.*log(mu)-mu}
logGammaPDF = @(x,a,b) -a*log(b)-log(gamma(a))+(a-1).*log(x)-x./b;
logVonMises = @(x,k,mu) k.*cos(x-mu)-log(besseli(0,k));

% INTIALIZE BOUNDS
Bounds = zeros(numParameters,2);

designBounds = repmat([-200,200],[designParams,1]); % two baselines, history, and movement

gaborBounds = zeros(gaborParams,2);
gaborBounds(1,:) = [-5000,5000]; % B for height of gabor
gaborBounds(2,:) = [-pi,pi]; % A orientation of Gabor
gaborBounds(3,:) = [min(xaxis)-50,max(xaxis)+50]; % x center
gaborBounds(4,:) = [min(yaxis)-50,max(yaxis)+50]; % y center
gaborBounds(5,:) = [1,2000]; % standard deviation x
gaborBounds(6,:) = [1,2000]; % standard deviation y
gaborBounds(7,:) = [0,100]; % spatial frequency
gaborBounds(8,:) = [-pi,pi]; %  orientation theta
gaborBounds(9,:) = [-pi,pi]; % phase shift phi

nonLinBounds = zeros(nonLinParams,2);
nonLinBounds(1,:) = [-1000,1000]; % sigmoid base
nonLinBounds(2,:) = [0,200]; % sigmoid slope
nonLinBounds(3,:) = [0,200]; % sigmoid rise

alphaBounds = repmat([-Inf,Inf],[precisionParams,1]); % for any precision parameter

% INITIALIZE INDEX VECTORS FOR EASY CALL TO THE CORRECT PARAMETERS
stimOffVec = 1;stimOnVec = 2;moveVec = 3;
designVec = 1:designParams;
bVec = designParams+1:designParams+gaborParams;
cVec = bVec(end)+1:bVec(end)+gaborParams*(numFilters-1);
nonLinVec = cVec(end)+1:cVec(end)+nonLinParams;
precisionVec = nonLinVec(end)+1:nonLinVec(end)+precisionParams;

Bounds(designVec,:) = designBounds;
Bounds(bVec,:) = gaborBounds;
Bounds(cVec,:) = repmat(gaborBounds,[numFilters-1,1]);
Bounds(nonLinVec,:) = nonLinBounds;
Bounds(precisionVec,:) = alphaBounds;

clear nonLinBounds gaborBounds historyMoveBaseBounds alphaBounds designBounds;

% INITIALIZE PRIORS

designPrior = zeros(length(designVec),2);
gaborPrior = zeros(gaborParams,2);
gaborPrior(1,:) = [0,0];
gaborPrior(2,:) = [0,1];
gaborPrior(3,:) = [0,1000]; % normal
gaborPrior(4,:) = [0,1000]; % normal
gaborPrior(5,:) = [20,15]; % gamma
gaborPrior(6,:) = [19,13]; % gamma


nonLinPrior = [0,0;1,0;2,0]; % base slope rise
precisionPrior = repmat([1e-3,1e3],[precisionParams,1]);

numIter = 16e5;burnIn = 1e5;skipRate = 1000;

PosteriorSamples = zeros(totalUnits,numParameters,length(burnIn+1:skipRate:numIter));
PosteriorMean = zeros(totalUnits,numParameters);
PosteriorInterval = zeros(totalUnits,numParameters,2);
for zz=1:totalUnits
        % ORGANIZE DATA
        spikeTrain = pointProcessSpikes(:,zz);
        historyDesign = zeros(length(spikeTrain),historyParams+2);
        for kk=1:historyParams-1
            temp = y;shift = zeros(kk,1);
            history = [shift;temp];
            historyDesign(:,kk+1) = history(1:(end-kk));
        end
        historyDesign(:,1) = stimOff;
        historyDesign(:,end-1) = stimOn;
        historyDesign(:,end) = movement;
        
        spikeTrain = spikeTrain(historyParams:end);
        historyDesign = historyDesign(historyParams:end,:);
        
       
        %MCMC intialization
        designPrior(stimOffVec,:) = [sum(spikeTrain)/length(spikeTrain),0];
        designPrior(stimOnVec,:) = [sum(spikeTrain)/length(spikeTrain),0];

        parameterVec = zeros(numParameters,numIter);
        posteriorProb = zeros(numIter,1);
        
        parameterVec(:,1) = max(Bounds(:,1),min(parameterVec(:,1),Bounds(:,2)));
        
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
                    rejectRate = rejectRate+1;
                end
            end
        end
        
        
        PosteriorSamples(zz,:,:) = parameterVec(:,burnIn:skipRate:numIter);
        
        alpha = 0.05;
        
        numColumns = 4;
        numRows = floor(numParameters/numColumns)+mod(numParameters,numColumns);
        for ii=1:numParameters
           subplot(numRows,numColumns,ii);
           histogram(squeeeze(PosteriorSamples(zz,ii,:)));
           title('Marginal Posterior, Parameter %d',ii);
           PosteriorMean(zz,ii) = mean(PosteriorSamples(ii,:));
           PosteriorInterval(zz,ii,:) = quantile(squeeze(PosteriorSamples(zz,ii,:)),...
               [alpha/2,1-alpha/2]);
        end
end


end

function [loglikelihood] = GetLikelihood(parameterVec)
global unbiasedS pointProcessStimTimes historyParams historyDesign ...
    stimulusTime X Y T gaborFun nonLinFun spikeTrain ...
    movement trainMillisecs;

baseMu = historyDesign(1:trainMillisecs,:)*parameterVec(1:historyParams);
moveMu = parameterVec(end).*movement(1:trainMillisecs);

gabor = gaborFun(X,Y,T,parameterVec(historyParams+1),parameterVec(historyParams+2),...
    parameterVec(historyParams+3),parameterVec(historyParams+4),parameterVec(historyParams+5),...
    parameterVec(historyParams+6),parameterVec(historyParams+7),...
    parameterVec(historyParams+8),parameterVec(historyParams+9),...
    parameterVec(historyParams+10),parameterVec(historyParams+11));
%filterEnergy = sum(gabor(:).*gabor(:));

filterOutput = zeros(trainMillisecs,1);
%stimLen = 29000:stimulusTime:trainMillisecs;

for kk=30000:stimulusTime:trainMillisecs
%     onScreenStim = squeeze(newS(pointProcessStimTimes(zz(kk)),:,:));
%     if mod(kk,stimulusTime) == 1
%         filterOutput = sum(newS(pointProcessStimTimes(zz(kk),:))'.*gaborFilter)./numPixels;
%     end
    filterOutput(kk-stimulusTime/2:kk+stimulusTime/2-1) = ...
        sum(sum(sum(unbiasedS(:,:,pointProcessStimTimes(kk-201:stimulusTime:kk))...
        .*gabor)));
end

filterOutput(1:30000) = sum(sum(sum(unbiasedS(:,:,pointProcessStimTimes(1:stimulusTime:201))...
        .*gabor)));

% for kk=1:numStimuli
%    filterOutput = sum(newS(kk,:)'.*gaborFilter)./numPixels;
%    loglikelihood(pointProcessStimTimes==kk) = Y(pointProcessStimTimes==kk).*log(baseMu*exp(filterOutput))-...
%        baseMu*exp(filterOutput);
% end

%figure();plot(filterOutput);

filterOutput = nonLinFun(filterOutput,parameterVec(historyParams+12),...
    parameterVec(historyParams+13),parameterVec(historyParams+14),...
    parameterVec(historyParams+15));


loglikelihood = sum(spikeTrain(1:trainMillisecs).*(baseMu+moveMu+filterOutput)-...
    exp(baseMu+moveMu+filterOutput));
end

