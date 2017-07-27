function [] = RF_PlaybackAnalysis(AnimalName,Date,NoiseType)
%RF_Playback_Analysis.m     
%    
% Created: 2016/07/26, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/07/26
% By: Byron Price

% read in the .plx file
beta = 0;

EphysFileName = strcat('RFPlaybackData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('RFPlaybackStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','adfreq','allad','svStrobed','tsevs')
load(StimulusFileName)


stimNums = unique(stimulusNumber(stimulusNumber~=0));

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

strobeData = strobeData(svStrobed~=0);
svStrobed = svStrobed(svStrobed~=0);

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


end