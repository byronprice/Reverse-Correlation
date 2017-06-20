function [PosteriorSamples,PosteriorMean,PosteriorInterval,Likelihood] =...
    BayesianFitSRF(AnimalName,Date,NoiseType)
%BayesianFitSRF.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatial receptive field for cells in V1
%  MCMC

% declare global variables
global numStimuli totalMillisecs pointProcessStimTimes h ...
    unbiasedS numParameters historyParams ...
    reducedMovement X Y gaborFun nonLinFun spikeTrain;
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','adfreq','allad','svStrobed','tsevs')
load(StimulusFileName)

gaborFun = @(x,y,A,xc,yc,sigmax,sigmay,spatFreq,theta,phi) ...
    exp(-((x-xc).*cos(A)-(y-yc).*sin(A)).^2./(2*sigmax*sigmax)-...
    ((x-xc).*sin(A)+(y-yc).*cos(A)).^2/(2*sigmay*sigmay))...
    .*sin((2*pi.*spatFreq).*(cos(theta-pi/2).*(x-xc)+sin(theta-pi/2).*(y-yc))-phi);

%nonLinFun = @(x,baseX,scale) exp((x-baseX)/scale);
nonLinFun = @(x,base,slope,rise,drop) rise./(1+exp((x-base).*slope));
%nonLinFun = @(x,slope,intercept) max(0,slope.*x+intercept);

% horzDegrees = atan((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
% vertDegrees = atan((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);

DIM = zeros(2,1);
if size(effectivePixels,2) == 1
    DIM(1) = sqrt(effectivePixels);
    DIM(2) = sqrt(effectivePixels);
elseif size(effectivePixels,2) == 2
    DIM(1) = effectivePixels(1);
    DIM(2) = effectivePixels(2);
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

unbiasedS = zeros(size(S));
for ii=1:numStimuli
    unbiasedS(ii,:) = reshape(real(ifftn(fftn(double(reshape(S(ii,:),[DIM(1),DIM(2)])))./S_f)),[DIM(1)*DIM(2),1])';
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
totalMillisecs = round(totalTime*1000);

stimTimes = round(strobeData.*1000);
pointProcessStimTimes = ones(totalMillisecs,numStimuli+1);

if exist('flipIntervals','var')==1
    flipInterval = mean(flipIntervals);
end

for kk=1:numStimuli-1
   pointProcessStimTimes(stimTimes(kk)+50:stimTimes(kk)+120,kk+1) = 1;
end
pointProcessStimTimes(stimTimes(numStimuli):(stimTimes(numStimuli)+1000*flipInterval)) = numStimuli+1;

temp = sum(pointProcessStimTimes(:,2:end),2);
stimOn = temp>0;

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = round(allts{ii}.*1000);
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

gaborParams = 8*21;
nonLinParams = 3;
historyParams = 101;
numParameters = historyParams+gaborParams+nonLinParams+3;

Bounds = zeros(numParameters,2);
Bounds(1:historyParams,:) = [-200,200]; % baseline & history
Bounds(historyParams+1,:) = [-pi,pi]; % A orientation of Gabor
Bounds(historyParams+2,:) = [min(xaxis)-50,max(xaxis)+50]; % x center
Bounds(historyParams+3,:) = [min(yaxis)-50,max(yaxis)+50]; % y center
Bounds(historyParams+4,:) = [1,500]; % standard deviation x
Bounds(historyParams+5,:) = [1,500]; % standard deviation y
Bounds(historyParams+6,:) = [0,100]; % spatial frequency
Bounds(historyParams+7,:) = [-pi,pi]; %  orientation theta
Bounds(historyParams+8,:) = [-pi,pi]; % phase shift phi
Bounds(historyParams+9,:) = [-1000,1000]; % sigmoid base
Bounds(historyParams+10,:) = [0,200]; % sigmoid slope
Bounds(historyParams+11,:) = [0,200]; % sigmoid rise
Bounds(historyParams+12,:) = [-200,200]; % parameter for stimulus on screen (DC gain)
Bounds(end,:) = [-200,200]; % parameter for movement modulation

numIter = 11e5;burnIn = 1e5;skipRate = 500;

PosteriorSamples = zeros(totalUnits,numParameters,length(burnIn+1:skipRate:numIter));
PosteriorMean = zeros(totalUnits,numParameters);
PosteriorInterval = zeros(totalUnits,numParameters,2);
for zz=1:totalUnits
        % ORGANIZE DATA
        spikeTrain = pointProcessSpikes(:,zz);
        historyDesign = zeros(length(spikeTrain),historyParams-1);
        for kk=1:historyParams
            temp = y;shift = zeros(kk,1);
            history = [shift;temp];
            historyDesign(:,kk) = history(1:(end-kk));
        end
        
        
        
        % RUN MCMC
        h = ones(numParameters,1)./1000;

        parameterVec = zeros(numParameters,numIter);
        posteriorProb = zeros(numIter,1);
        
        parameterVec(1,1) = b(1);
        parameterVec(end,1) = b(2);
        
        parameterVec(historyParams+1,1) = normrnd(0.5,0.1);
        parameterVec(historyParams+2,1) = normrnd(15,5);
        parameterVec(historyParams+3,1) = normrnd(0,0.1);
        parameterVec(historyParams+4,1) = normrnd(pi,pi/4);
        parameterVec(historyParams+5,1) = normrnd(0,300);
        parameterVec(historyParams+6,1) = normrnd(0,200);
        parameterVec(historyParams+7,1) = normrnd(150,50);
        parameterVec(historyParams+8,1) = normrnd(150,50);
        parameterVec(historyParams+9,1) = normrnd(0.4,0.1);
        parameterVec(historyParams+10,1) = normrnd(pi,pi/4);
        parameterVec(historyParams+11,1) = normrnd(pi,pi/4);
        parameterVec(historyParams+12,1) = normrnd(2,1);
        parameterVec(historyParams+13,1) = normrnd(5,0.5);
        parameterVec(historyParams+14,1) = normrnd(1,0.5);
        parameterVec(historyParams+15,1) = normrnd(0,1);
        
        parameterVec(:,1) = max(Bounds(:,1),min(parameterVec(:,1),Bounds(:,2)));
        
        likelihoodXprev = GetLikelihood(parameterVec(:,1));
        
        for iter = 2:numIter
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
