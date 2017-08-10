function [] = RevCorrMov(AnimalName,Date,NoiseType)
%RevCorrMov.m
%   %   Analysis of single unit recording data in response to a white
%   or pink noise movie (see Noise_Movie.m for Psychtoolbox
%   stimulus)
%    See Smyth et al. 2003 Receptive Field Organization ... We'll solve a
%    matrix equation using a regularized pseudo-inverse technique.
%
% Briefly, if the receptive field for a given neuron is described by a
%  p-by-1 vector (p being the number of pixels), f, and the stimuli described
%  by a s-by-p (s being the number of stimuli presented) matrix, S, and the
%  response of that neuron described by an s-by-1 vector, r, then we can write
%  the responses in terms of the stimuli and the receptive field as
%  r = S * f .  We want to infer the spatial receptive field structure
%  f, so we need S^-1 * r = S^-1 * S * f => f = S^-1 * r ... So, the
%  following code will get the spatial receptive fields f for a series of
%  neurons whose responses were measured with single-unit recordings during
%  the display of the white/pink noise stimuli in S.
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%       NoiseType - 'white' or 'pink' or 'brown'

%OUTPUT: saves a file with the spatiotemporal receptive field estimate
%
% Created: 2016/04/29, Commuter Rail Boston to Providence
%  Byron Price
% Updated: 2017/08/09
% By: Byron Price

%cd('~/CloudStation/ByronExp/NoiseRetino');
% read in the .plx file
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseMovieData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseMovieStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
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
    
xaxis = linspace(-round(screenPix_to_effPix*DIM(2)/2)+1,...
    round(screenPix_to_effPix*DIM(2)/2),DIM(2));
yaxis = linspace(round(3*screenPix_to_effPix*DIM(1)/4),...
    -round(screenPix_to_effPix*DIM(1)/4)+1,DIM(1));
taxis = 300:16:-40;

[X,Y,T] = meshgrid(xaxis,yaxis,taxis);

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

strobeData = strobeData(strobeData>0);

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
    movement = logical(tempMov>(3*stdEst));
    clear tempMov stdEst;
else
    movement = zeros(length(LFP{1}),1);
end

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
timeMultiplier = 1000;
totalMillisecs = round(totalTime*timeMultiplier);

stimTimes = round(strobeData.*timeMultiplier);

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = max(1,round(allts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

reducedMov = zeros(numStimuli,1*timeMultiplier);
reducedSpikeCount = zeros(totalUnits,numStimuli,1*timeMultiplier);
for ii=1:totalUnits
    for jj=1:numStimuli
        timeInds = stimTimes(jj):stimTimes(jj)+1*timeMultiplier-1;
        reducedSpikeCount(ii,jj,:) = pointProcessSpikes(timeInds,ii)';
        reducedMov(jj,:) = movement(timeInds);
    end
end

clearvars -except EphysFileName totalUnits numStimuli ...
    reducedSpikeCount DIM unbiasedS allts strobeData xaxis yaxis ...
    X Y totalMillisecs reducedMov movement;

fullSize = DIM(1)*DIM(2);
basisStdDevs = [200,300,400,450,500,550,600,650,700,750,1000];
numStdDevs = length(basisStdDevs);

finalResultsPoissonB = cell(totalUnits,numStdDevs);
finalResultsPoissonDev = cell(totalUnits,numStdDevs);
finalResultsPoissonSE = cell(totalUnits,numStdDevs);
finalResultsNormalB = cell(totalUnits,numStdDevs);
finalResultsNormalDev = cell(totalUnits,numStdDevs);
finalResultsNormalSE = cell(totalUnits,numStdDevs);
for zz=1:totalUnits
   spikeTrain = squeeze(reducedSpikeCount(zz,:,:));
   spikeTrain = sum(spikeTrain(:,50:500),2);
   movDesign = sum(reducedMov(:,50:500),2);
   
   %    r = spikeTrain;
   %    fhat = unbiasedS\r;
   gaussFun = @(x,y,xc,yc,std) exp(-((x-xc).*(x-xc))./(2*std*std)-...
       ((y-yc).*(y-yc))./(2*std*std));
   
   center1 = xaxis(1:3:end);
   center2 = yaxis(1:3:end);
   numBasis1 = length(center1);
   numBasis2 = length(center2);
   totalParams = numBasis1*numBasis2;
   
   basisFuns = zeros(fullSize,totalParams);
   for stddev = 1:numStdDevs
       count = 1;
       for ii=1:numBasis1
           for jj=1:numBasis2
               temp = gaussFun(X,Y,center1(ii),center2(jj),basisStdDevs(stddev));
               basisFuns(:,count) = temp(:)./max(temp(:));
               count = count+1;
           end
       end
       design = [movDesign,unbiasedS*basisFuns];
       [bPoiss,devPoiss,statsPoiss] = glmfit(design,spikeTrain,'poisson');
       finalResultsPoissonB{zz,stddev} = bPoiss;
       finalResultsPoissonDev{zz,stddev} = devPoiss;
       finalResultsPoissonSE{zz,stddev} = statsPoiss.se;
       
       [b,dev,stats] = glmfit(design,spikeTrain);
       finalResultsNormalB{zz,stddev} = b;
       finalResultsNormalDev{zz,stddev} = dev;
       finalResultsNormalSE{zz,stddev} = stats.se;
       
       clear bPoiss devPoiss statsPoiss b dev stats count design;
   end
   clear spikeTrain movDesign basisFuns center1 center2 numBasis1 numBasis2 totalParams;
end

fileName = strcat(EphysFileName(1:end-9),'-GLMResults.mat');
save(fileName,'finalResultsNormalB','finalResultsPoissonB',...
    'finalResultsNormalDev','finalResultsPoissonDev',...
    'finalResultsNormalSE','finalResultsPoissonSE','allts','totalUnits',...
    'reducedSpikeCount','DIM','unbiasedS','movement',...
    'xaxis','yaxis','basisStdDevs','reducedMov');
%cd('~/Documents/Current-Projects/Reverse-Correlation');

end
