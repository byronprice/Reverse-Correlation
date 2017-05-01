function [] = RevCorrMovies2(AnimalName,Date,NoiseType)
%RevCorrMovies.m
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
% Updated: 2017/04/29
% By: Byron Price

cd('~/CloudStation/ByronExp/NoiseRetino');

% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseMovieData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseMovieStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'allad','allts','adfreq','tsevs','nunits1')
load(StimulusFileName)

% CREATE UNBIASED VERSION OF MOVIE BY DIVIDING OUT POWER SPECTRUM
%  USED TO GENERATE THE MOVIE
DIM = [effectivePixels(1),effectivePixels(2),numStimuli];

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]'/DIM(2);
t = [(0:floor(DIM(3)/2)) -(ceil(DIM(3)/2)-1:-1:1)]'/(DIM(3));
[U,V,T] = meshgrid(u,v,t);
S_f = single((U.^spaceExp+V.^spaceExp+T.^timeExp).^(beta/2));

S_f(S_f==inf) = 1;

newS = real(ifftn(fftn(S)./S_f));

newS = single(reshape(newS,[DIM(1)*DIM(2),numStimuli]))';

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
LFP_Movement = cell(length(inds),1);
for ii=1:length(inds)
    LFP_Movement{ii} = allad{inds};
end

totalTime = length(LFP_Movement{1})./adfreq;
movement = LFP_Movement{end};
difference = length(LFP_Movement{1})-length(LFP_Movement{end});

if mod(difference,2) == 0
   addOn = difference/2;
   movement = [zeros(addOn,1);movement;zeros(addOn,1)];
else
   addOn = floor(difference/2); 
   movement = [zeros(addOn,1);movement;zeros(addOn+1,1)];
end

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
totalMillisecs = round(totalTime*1000);

stimTimes = round(strobeData.*1000);
pointProcessStimTimes = zeros(totalMillisecs,1);

for kk=1:numStimuli-1
    pointProcessStimTimes(stimTimes(kk):stimTimes(kk+1)-1) = kk;
end
pointProcessStimTimes(stimTimes(numStimuli):(stimTimes(numStimuli)+1000/movie_FrameRate)) = numStimuli;

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
    spikeTimes = round(allts{ii}.*1000);
    for jj=1:length(spikeTimes)
       pointProcessSpikes(spikeTimes(jj),ii) = 1;
    end
end

kernelLen = 0.25;
stimLen = 0.03;

Response = cell(totalUnits,2);
for ii=1:totalUnits
        firstStim = strobeData(1)+kernelLen;
        lastStim = strobeData(end)+kernelLen;
        
        totalStims = round((lastStim-firstStim)/stimLen);
        stimOffsets = linspace(firstStim,lastStim,totalStims);
        Response{ii,1} = zeros(totalStims,1,'single');
        for kk=1:totalStims
            stimOffset = stimOffsets(kk);
            temp = sum(pointProcessSpikes(round((stimOffset)*1000):round((stimOffset+stimLen)*1000),ii));
            
            Response{ii,1}(kk,1) = temp;
            Response{ii,2}(kk,2) = pointProcessStimTimes(round((stimOffset-kernelLen)*1000):...
                round((stimOffset)*1000));
        end
end

clear pointProcessStimTimes pointProcessSpikes totalMillisecs strobeData;


% CREATE LAPLACIAN MATRIX
L = zeros(DIM(1)*DIM(2),DIM(1)*DIM(2),'single');

%operator = [0,-1,0;-1,4,-1;0,-1,0];
bigCount = 1;
for jj=1:DIM(1)
    for ii=1:DIM(2)
        tempMat = zeros(DIM(2),DIM(1));
        tempMat(ii,jj) = 4;
        if ii > 1
            tempMat(ii-1,jj) = -1;
        end
        if ii < DIM(2)
            tempMat(ii+1,jj) = -1;
        end
        if jj > 1
            tempMat(ii,jj-1) = -1;
        end
        if jj < DIM(1)
            tempMat(ii,jj+1) = -1;
        end
        L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
    end
end

% REGULARIZED PSEUDO-INVERSE SOLUTION
horzDegrees = atand((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
vertDegrees = atand((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);
xaxis = linspace(-horzDegrees/2,horzDegrees/2,DIM(1));
yaxis = linspace(-vertDegrees/4,3*vertDegrees/4,DIM(2));

bigLambda = logspace(3,6,10);
F = zeros(totalUnits,numLags,DIM(1)*DIM(2));

parfor ii=1:totalUnits
    for jj=1:numLags
        r = squeeze(Response{ii,jj}(:,1));
        constraints = [r;zeros(DIM(1)*DIM(2),1,'single')];
        onScreenMovie = single(newS(Response{ii,jj}(:,2),:));
        
        tempF = zeros(length(bigLambda),DIM(1)*DIM(2));
        RMS =  zeros(length(bigLambda),1);
        for lambda = 1:length(bigLambda)
            % VERTICALLY CONCATENATE S and L
            % onScreenMovie is size numStimuli X effectivePixels
            A = [onScreenMovie;bigLambda(lambda).*L];
            fhat = pinv(A)*constraints;
            
            tempF(lambda,:) = fhat;
            RMS(lambda) = norm(r-onScreenMovie*fhat)./sqrt(DIM(1)*DIM(2));
        end

        rmsDiff = diff(RMS');
        lambdaDiff = diff(bigLambda);
        deltaRMSdeltaLambda = rmsDiff./lambdaDiff;
        [maxVal,index] = max(abs(deltaRMSdeltaLambda));
        onepercent = 0.01.*maxVal;
        firstBelow = find(abs(deltaRMSdeltaLambda(index:end)<onepercent,1));
        bestMap = index+firstBelow-1;
        F(ii,jj,:) = tempF(bestMap,:);
     end
end


for ii=1:totalUnits
    h(ii) = figure();
    for jj=1:numLags
        subplot(numLags/2,2,jj);
        imagesc(xaxis,yaxis,reshape(F(ii,jj,:),[DIM(2),DIM(1)]));set(gca,'YDir','normal');
        title(sprintf('Unbiased STRF Unit %d, Lag %3.0f ms',ii,lagTimes(jj)*1000));
        xlabel('Azimuth (dva)');ylabel('Altitude (dva)');
    end
end

FileName = strcat('NoiseMovieResults',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
save(FileName,'F','newS','Response','allts','bigLambda',...
    'beta','totalUnits','h','stimLen');

end