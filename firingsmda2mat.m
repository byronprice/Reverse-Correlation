function [] = firingsmda2mat(directory)
% convert the firings.mda file output from MountainSort back to MATLAB
%  format for data analysis

% fileID = fopen('fileName.txt','r');
% fileName1 = fscanf(fileID,'%s');

cd directory

files = dir('firing*.mda');
numFiles = length(files);

spikeTimes = cell(numFiles,1);
totalUnits = 0;
for ii=1:numFiles
    firingFile = files(ii).name;
    A = readmda(firingFile);
    
    channels = unique(A(1,:));
    numChannels = length(channels);
    
    units = unique(A(3,:));
    numUnits = length(units);
    
    spikeTimes{ii} = cell(numUnits,1);
    
    count = 1;
    for jj=1:numChannels
        pseudo_event_times = A(2,A(1,:)==channels(jj));
        unitcode = A(3,A(1,:)==channels(jj));
        unitIDs = unique(unitcode);
        nunits = length(unitIDs);
        for kk=1:nunits
            spikeTimes{ii}{count} = pseudo_event_times(unitcode==unitIDs(kk))+1;
            totalUnits = totalUnits+1;
            count = count+1;
        end
    end
end

clearvars -except spikeTimes numFiles

!cd ..

fileName1 = strcat(directory,'-mda.mat');
load(fileName1)

% convert from pseudo event times to experimental time
newts = cell(totalUnits,1);
newwaves = cell(totalUnits,1);

count = 1;
for ii=1:numFiles
   spikeTimeTempArray = spikeTimes{ii};
   numUnits = size(spikeTimeTempArray,1);
   
   trueEventTimes = allts{ii};
   trueWaves = allwaves{ii};
   correspondingEventIndices = double(allEventTimes{ii}');
   for jj=1:numUnits 
       pseudo_event_times = unique(spikeTimeTempArray{jj});
       indices = zeros(length(pseudo_event_times),1);
       for kk=1:length(pseudo_event_times)
          [difference,ind] = min(abs(correspondingEventIndices-pseudo_event_times(kk)));
          if abs(difference)<5
            indices(kk) = ind;
          end
       end
       indices = indices(indices>0);
       newwaves{count} = trueWaves(indices,:);
       newts{count} = trueEventTimes(indices);count=count+1;
   end
end
clear difference ind trueWaves count;

% include and exclude units
%  criterion for inclusion: 
%     1) <=1% of spikes within a 2ms refractory period
%     2) <0.9 correlation with all other recorded neurons
%     3) neurons with significant negative correlation may be merged

timeMultiplier = 1000;
nonemptyad = ~cellfun(@isempty,allad);
temp = allad(nonemptyad);temp2 = adfreqs(nonemptyad);
totalTime = length(temp{1})/temp2(1);

pointProcessSpikes = zeros(round(totalTime*timeMultiplier),totalUnits);
for ii=1:totalUnits
   spikeTimes = max(1,round(newts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = 1;
   end
end

refractory_cutoff = 2/1000;
refractory_inclusion = 0.01;
spikeHz_cutoff = 0.1;spikeNum_cutoff = spikeHz_cutoff*totalTime;
correlation_inclusion = [-0.1,0.8];
toInclude = ones(totalUnits,1);

for ii=1:totalUnits
   spikeTimes = newts{ii};
   for jj=ii+1:totalUnits
       [r,~] = corrcoef(pointProcessSpikes(:,ii),pointProcessSpikes(:,jj));
       if r(1,2) < correlation_inclusion(1)
          toInclude(jj) = 0; 
          pointProcessSpikes(:,jj) = 0;
          temp = newts{jj};
          newts{ii} = [spikeTimes;temp];
          newts{jj} = 0;
       end
   end
end

for ii=1:totalUnits
    spikeTimes = newts{ii};
    if length(spikeTimes) < spikeNum_cutoff
        toInclude(ii) = 0;
    end
    isi = diff([0;spikeTimes]);
%     figure();subplot(2,1,1);plot(spikeTimes);
%     subplot(2,1,2);histogram(isi);
    criterion1 = sum(isi<refractory_cutoff)/length(isi);
    fprintf('\nProportion refractory violations: %3.2e\n',criterion1);
    if criterion1 > refractory_inclusion
        toInclude(ii) = 0;
    end
    for jj=ii+1:totalUnits
       [r,~] = corrcoef(pointProcessSpikes(:,ii),pointProcessSpikes(:,jj));
%        fprintf('Neuronal Correlation: %3.2e\n',r(1,2));
       if r(1,2) > correlation_inclusion(2)
          toInclude(jj) = 0; 
       end
    end
end

totalUnits = sum(toInclude);
allts = cell(totalUnits,1);
allwaves = cell(totalUnits,1);
meanWaves = cell(totalUnits,1);
quantile95Waves = cell(totalUnits,1);
alpha = 0.05;

inds = find(toInclude==1);
for ii=1:totalUnits
   allts{ii} = newts{inds(ii)};
   allwaves{ii} = newwaves{inds(ii)};
   
   temp = allwaves{ii};
   meanWaves{ii} = mean(temp,1);
   quantile95Waves{ii} = quantile(temp,[alpha/2,1-alpha/2]);
end

clear pointProcessSpikes temp temp2 ii jj timeMultiplier newts spikeTimes ...
    trueEventTimes spikeTimeTempArray correspondingEventIndices toInclude r ...
    criterion1 allEventTimes inds isi kk index indices nonemptyad pseudo_event_times ...
    newwaves;

fprintf('\nTotal Units: %d\n',totalUnits);

pause(1);
newFileName = sprintf('%s-mounsort.mat',fileName1(1:end-8));
save(newFileName);


cd('/home/byron/CloudStation/ByronExp/NoiseRetino_MountainSort');

exit;
end