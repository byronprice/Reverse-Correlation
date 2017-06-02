function [] = BayesianReceptiveFields(AnimalName,Date,NoiseType,startUnit)
%BayesianReceptiveFields.m
%   %   Analysis of single unit recording data in response to a white
%   or pink noise movie (see Noise_Movie.m for Psychtoolbox
%   stimulus)
%    See Bayesian Spike-Triggered Covariance from Jonathan Pillow
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%       NoiseType - 'white' or 'pink' or 'brown'

%OUTPUT: saves a file with the spatiotemporal receptive field estimate
%
% Created: 2016/06/01, 24 Cummington Mall, Boston
%  Byron Price
% Updated: 2017/06/01
% By: Byron Price

%cd('~/CloudStation/ByronExp/NoiseRetino');
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseMovieData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort2.mat');

% if exist(EphysFileName,'file') ~= 2
%    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
% end

StimulusFileName = strcat('NoiseMovieStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'allad','allts','adfreq','tsevs','nunits1')
load(StimulusFileName)

% CREATE UNBIASED VERSION OF MOVIE BY DIVIDING OUT POWER SPECTRUM
%  USED TO GENERATE THE MOVIE
DIM = [effectivePixels(1),effectivePixels(2),numStimuli];

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
t = [(0:floor(DIM(3)/2)) -(ceil(DIM(3)/2)-1:-1:1)]'/(DIM(3));
[V,U,T] = meshgrid(v,u,t);
S_f = single((U.^spaceExp+V.^spaceExp+T.^timeExp).^(beta/2));

S_f(S_f==inf) = 1;

unbiasedS = real(ifftn(fftn(double(S))./S_f));

unbiasedS = reshape(unbiasedS,[DIM(2)*DIM(1),numStimuli]);

clear S S_f U V T u v t;

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

totalTime = (length(LFP{1})./adfreq)/10;
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
    movement = movement(1:10:end);
    clear tempMov stdEst;
end

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
timeMultiplier = 100;
totalCentisecs = round(totalTime*timeMultiplier);
timeVector = linspace(1/timeMultiplier,totalTime,totalCentisecs);

stimTimes = round(strobeData.*timeMultiplier);
pointProcessStimTimes = zeros(totalCentisecs,1);

for kk=1:numStimuli-1
   pointProcessStimTimes(stimTimes(kk):stimTimes(kk+1)-1) = kk;
end
pointProcessStimTimes(stimTimes(numStimuli):(stimTimes(numStimuli)+timeMultiplier/movie_FrameRate)) = numStimuli;

pointProcessSpikes = zeros(totalCentisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = round(allts{ii}.*timeMultiplier);
   for jj=1:length(spikeTimes)
       if spikeTimes(jj) > 0 
          pointProcessSpikes(spikeTimes(jj),ii) = pointProcessSpikes(spikeTimes(jj),ii)+1;
       end
   end
end

filterTime = 0.2;filterSteps = 1/timeMultiplier;
filterLen = filterTime/filterSteps+1;

onScreenInds = zeros(totalCentisecs,filterLen,'single');
for kk=filterLen+1:totalCentisecs
    stimOffset = timeVector(kk);
    temp = pointProcessStimTimes(round((stimOffset-filterTime)*timeMultiplier):...
            round(filterSteps*timeMultiplier):round(stimOffset*timeMultiplier));
    onScreenInds(kk,:) = temp;
end

fprintf('Gathered movie indices\n');

fileName = strcat('NoiseMovieConvertedData',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
save(fileName,'unbiasedS','Response','allts','totalUnits',...
    'movieTime_Seconds','LFP','movement','DIM','pointProcessStimTimes',...
    'pointProcessSpikes','totalTime','totalStims','filterTime','filterLen',...
    'totalCentisecs','beta','filterSteps','movie_FrameRate');

clear beta;
% MCMC
for ii=1:totalUnits
    y = pointProcessSpikes(:,ii);
    historyParams = 10;
    historyDesign = zeros(totalCentisecs,historyParams);
    
    for kk=1:historyParams
        temp = y;shift = zeros(kk,1);
        history = [shift;temp];
        historyDesign(:,kk) = history(1:(end-kk));
    end
    y = y(filterLen+1:end);
    historyDesign = historyDesign(filterLen+1:end,:);
    
    X = [ones(length(historyDesign),1),historyDesign];
    onScreenInds = onScreenInds(filterLen+1:end,:);
    
    logPoissonPDF = @(y,mu) y.*mu-exp(mu);
    filterSize = DIM(1)*DIM(2)*filterLen;
    Q = 10;
    
    baseFiring = sum(y)/length(y);
    numIter = 1.5e6;burnIn = 5e5;numParams = historyParams+1+(Q+1)*filterSize+2*(Q+1);
    skipRate = 1000;totalSamples = (numIter-burnIn)/skipRate;
    params = zeros(numParams,totalSamples);
    updateMu = zeros(numParams,1);
    updateSigma = eye(numParams);
    posteriorProb = zeros(totalSamples,1);
    
    bVec = (historyParams+2):(history+2+filterSize);
    designVec = 1:(historyParams+1);
    
    baseMu = log(baseFiring);
    historyPriorD = (historyParams)/2;
    filterPriorD = filterSize/2;
    
    abprior = [1e-3,1e3];
    
    % intialize params and recursive mean estimation vector
    params(1,1) = normrnd(baseMu,1);
    params(2:historyParams+1,1) = normrnd(0,1,[historyParams,1]);
    
    params(historyParams+2:historyParams+2+(Q+1)*filterSize,1) = normrnd(0,1,[(Q+1)*filterSize,1]);
    params(historyParams+2+(Q+1)*filterSize+1:end,1) = normrnd(0,1,[2*(Q+1),1]);
    
    updateMu(1,1) = normrnd(baseMu,1);
    updateMu(2:historyParams+1,1) = normrnd(0,1,[historyParams,1]);
    
    updateMu(historyParams+2:historyParams+2+(Q+1)*filterSize,1) = normrnd(0,1,[(Q+1)*filterSize,1]);
    updateMu(historyParams+2+(Q+1)*filterSize+1:end,1) = normrnd(0,1,[2*(Q+1),1]);
    
    b = params(historyParams+2:historyParams+2+filterSize,1);
    W = params(historyParams+3+filterSize);
end

horzDegrees = atand((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
vertDegrees = atand((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);
xaxis = linspace(-horzDegrees/2,horzDegrees/2,DIM(1));
yaxis = linspace(-vertDegrees/4,3*vertDegrees/4,DIM(2));
taxis = linspace(-filterTime*1000,0,filterLen);


fileName = sprintf('NoiseMovieResults%s%d_%d.mat',NoiseType,Date,AnimalName);
save(fileName,'F','Response','bigLambda','totalUnits','spikeCountLen',...
    'xaxis','yaxis','taxis','DIM','kernelLenFull','DistToScreen',...
    'screenPix_to_effPix','kernelLen','onScreenInds','RMS','totalStims');

%cd('~/Documents/Current-Projects/Reverse-Correlation');

end
