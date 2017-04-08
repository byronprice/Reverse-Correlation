function [ output_args ] = FitReceptiveField_SimpleCell(AnimalName,Date)
%FitReceptiveField_SimpleCell.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatial receptive field for a simple cell in L4 of V1.

% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName)
load(StimulusFileName)

% generate theoretical stimulus power spectrum
N = sqrt(effectivePixels);
u = 0:(N-1);
[U,V] = meshgrid(u,u);
S_f = (U.^spaceExp+V.^spaceExp).^(beta/2);
S_f(S_f==inf) = 1;


temp = ~cellfun(@isempty,allts);
Chans = find(sum(temp,1));numChans = length(Chans);

% tsevs are the strobed times of stimulus onset, then offset
%  Onset at tsevs{1,33}(2), offset at tsevs{1,33}(3), onset at
%  tsevs{1,33}(4), offset at 5, etc.

%x = find(~cellfun(@isempty,tsevs));
strobeStart = 3;

temp = cell(numChans,nunits1);
for ii=1:numChans
    for jj=1:nunits1
        if isempty(allts{jj,Chans(ii)}) == 0
            temp{ii,jj} = allts{jj,Chans(ii)};
        end
    end
end
allts = temp;
fullSpots = 1-cellfun(@isempty,allts);
Neurons = sum(fullSpots,2);

% ASSUME THAT A RESPONSE TO A STIMULUS OFFSET IS THE SAME AS A RESPONSE TO
%  THE NEGATIVE OF THAT IMAGE, image created with values from 0 to 255
Grey = 127;
numStimuli = numStimuli*2;
newS = zeros(numStimuli,size(S,2),'single');
for ii=1:numStimuli
    if mod(ii,2) == 1
        newS(ii,:) = S(floor(ii/2)+1,:);
    elseif mod(ii,2) == 0
        newS(ii,:) = 255-S(ii/2,:); %255-S(ii/2,:)
%         temp = reshape(S(ii/2,:),[N,N]);
%         temp2 = reshape(255-S(ii/2,:),[N,N]);
%         figure();imagesc(temp);
%         figure();imagesc(temp2);
%         display('blah');
    end
end
newS = Grey-newS;
clear S;

strobeData = tsevs{1,strobeStart};

historyParams = 25;
gaborParams = 8;
numParameters = historyParams+gaborParams;
finalParameters = zeros(numChans,nunits1,numParameters);
fisherInfo = zeros(numChans,nunits1,numParameters,numParameters);
ninetyfiveErrors = zeros(numChans,nunits1,numParameters);

maxITER = 100;
likelyTolerance = 1e-6;
gradientTolerance = 1e-6;
numRepeats = 5e4;

for zz=1:numChans
    numNeurons = Neurons(zz);
    for yy=1:numNeurons
        totalTime = strobeData(end)+2;
        totalMillisecs = round(totalTime*1000);
        
    end
end


end

