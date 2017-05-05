function [] = RevCorrMov2(AnimalName,Date,NoiseType)
%RevCorrMov2.m
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
% Updated: 2017/05/03
% By: Byron Price

%cd('~/CloudStation/ByronExp/NoiseRetino');
tic;
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseMovieData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

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
       if spikeTimes(jj) > 0 
          pointProcessSpikes(spikeTimes(jj),ii) = 1;
       end
   end
end

numStimuli = 10000;
kernelLen = 0.2;
spikeCountLen = 0.02;
kernelSteps = 0.01;
kernelLenFull = kernelLen/kernelSteps+1;

firstStim = strobeData(1)+kernelLen;
lastStim = strobeData(end);

stimOffsets = unifrnd(firstStim,lastStim,[numStimuli,1]);

if exist('movement','var')
    
    Response = cell(totalUnits,1);
    for ii=1:totalUnits
        Response{ii,1} = zeros(numStimuli,2,'single');
        for kk=1:numStimuli
            stimOffset = stimOffsets(kk);
            Response{ii,1}(kk,1) = sum(pointProcessSpikes(round(stimOffset*1000):...
                round((stimOffset+spikeCountLen)*1000),ii));
            temp = movement(round(stimOffset*1000):round((stimOffset+spikeCountLen)*1000));
            Response{ii,1}(kk,2) = sum(temp)<=((spikeCountLen*1000)/10);
        end
    end
    
else
    Response = cell(totalUnits,1);
    for ii=1:totalUnits
        
        Response{ii,1} = zeros(numStimuli,2,'single');
        for kk=1:numStimuli
            stimOffset = stimOffsets(kk);
            Response{ii,1}(kk,1) = sum(pointProcessSpikes(round(stimOffset*1000):...
                round((stimOffset+spikeCountLen)*1000),ii));
            Response{ii,1}(kk,2) = 1;
        end
    end
end

onScreenInds = zeros(sum(Response{1}(:,2)),kernelLenFull,'single');
count = 1;
for kk=1:numStimuli
    if Response{1}(kk,2) == 1
        stimOffset = stimOffsets(kk);
        temp = pointProcessStimTimes(round((stimOffset-kernelLen)*1000):...
            round(kernelSteps*1000):round(stimOffset*1000));
        onScreenInds(count,:) = temp;
        count = count+1;
    end
end

totalStims = sum(Response{1}(:,2));
fprintf('Gathered movie indices\n');

fileName = strcat('NoiseMovieConvertedData',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
save(fileName,'unbiasedS','Response','allts','totalUnits',...
    'movieTime_Seconds','LFP','movement','DIM','pointProcessStimTimes',...
    'pointProcessSpikes','totalTime','totalStims','kernelLen','kernelLenFull',...
    'totalMillisecs','beta','kernelSteps','movie_FrameRate','stimOffsets');

clearvars -except unbiasedS onScreenInds Response DIM kernelLenFull screenPix_to_effPix ...
    DistToScreen totalStims conv_factor kernelLen totalUnits Date AnimalName NoiseType ...
    spikeCountLen;


% CREATE LAPLACIAN MATRIX
%L = zeros(DIM(1)*DIM(2),DIM(1)*DIM(2),'single');
L = zeros(DIM(1)*DIM(2),DIM(1)*DIM(2),'single');

%2D operator = [0,-1,0;-1,4,-1;0,-1,0];
bigCount = 1;
for jj=1:DIM(2)
    for ii=1:DIM(1)
        tempMat = zeros(DIM(1),DIM(2));
        tempMat(ii,jj) = 4;
        if ii > 1
            tempMat(ii-1,jj) = -1;
        end
        if ii < DIM(1)
            tempMat(ii+1,jj) = -1;
        end
        if jj > 1
            tempMat(ii,jj-1) = -1;
        end
        if jj < DIM(2)
            tempMat(ii,jj+1) = -1;
        end
        L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
    end
end

fprintf('Populated matrix L\n');
clear tempMat ii jj kk;

% REGULARIZED PSEUDO-INVERSE SOLUTION

bigLambda = [0,5,10];
F = zeros(totalUnits,kernelLenFull,DIM(1)*DIM(2));
train = round(totalStims*0.7);
for ii=1:totalUnits
    fprintf('Running unit: %d\n',ii);
    r = squeeze(Response{ii}(Response{ii}(:,2)==1,1));
    constraints = [r(1:train);zeros(DIM(1)*DIM(2),1)];
    
    tempF = zeros(length(bigLambda),kernelLenFull,DIM(1)*DIM(2));
    RMS =  zeros(length(bigLambda),1);
    for lambda = 1:length(bigLambda)
        fprintf('Lambda: %3.2e\n',bigLambda(lambda));
            % VERTICALLY CONCATENATE S and L
            % onScreenMovie is size numStimuli X effectivePixels
        tic;
        for kk=1:kernelLenFull
            inds = onScreenInds(1:train,kk);
            tempF(lambda,kk,:) = double([unbiasedS(:,inds)';bigLambda(lambda).*L])\constraints;
        end
        toc;
        fhat = squeeze(tempF(lambda,:,:))';
        inds = onScreenInds(train+1:totalStims,:);
        A = zeros(length(inds),DIM(1)*DIM(2)*kernelLenFull);
        for ll=1:length(inds)
            tempInds = inds(ll,:);
            tempMov = unbiasedS(:,tempInds);
            A(ll,:) = tempMov(:);
        end
        RMS(lambda) = norm(r(train+1:totalStims)-A*fhat(:)./kernelLenFull)...
                ./sqrt(DIM(1)*DIM(2)*kernelLenFull);
        clear A tempMov tempInds;
    end
    [~,bestMap] = min(RMS);
    F(ii,:,:) = tempF(bestMap,:,:);
end

horzDegrees = atand((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
vertDegrees = atand((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);
xaxis = linspace(-horzDegrees/2,horzDegrees/2,DIM(1));
yaxis = linspace(-vertDegrees/4,3*vertDegrees/4,DIM(2));
taxis = linspace(-kernelLen*1000,0,kernelLenFull);


% for ii=1:totalUnits
%     h(ii) = figure();
%     for jj=1:kernelLenFull
%         subplot(numLags/2,2,jj);
% imagesc(xaxis,yaxis,reshape(F(ii,jj,:),[DIM(2),DIM(1)]));set(gca,'YDir','normal');
%         title(sprintf('Unbiased STRF Unit %d, Lag %3.0f ms',ii,taxis(jj)*1000));
%         xlabel('Azimuth (dva)');ylabel('Altitude (dva)');
%     end
% end
runTime = toc;

fileName = strcat('NoiseMovieResults',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
save(fileName,'F','Response','bigLambda','totalUnits','spikeCountLen',...
    'xaxis','yaxis','taxis','DIM','kernelLenFull','DistToScreen',...
    'screenPix_to_effPix','kernelLen','runTime','onScreenInds');

%cd('~/Documents/Current-Projects/Reverse-Correlation');

end
