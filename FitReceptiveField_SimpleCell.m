function [PosteriorSamples,PosteriorMean,PosteriorInterval] = FitReceptiveField_SimpleCell(AnimalName,Date,NoiseType)
%FitReceptiveField_SimpleCell.m
%  Use data from a receptive-field mapping experiment and fit a Gabor model
%   to obtain the spatial receptive field for a simple cell in L4 of V1.

% declare global variables
global xmesh ymesh numPixels numStimuli totalMillisecs Y pointProcessStimTimes h ...
    historyDesign newS numParameters historyParams gaborParams stimulusTime;
% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'nunits1','allts','svStrobed','tsevs')
load(StimulusFileName)

% generate theoretical stimulus power spectrum
N = sqrt(effectivePixels);
u = 0:(N-1);v = 0:(N-1);
[U,V] = meshgrid(u,v);
S_f = (U.^spaceExp+V.^spaceExp).^(beta/2);
S_f(S_f==inf) = 1;

xaxis = linspace(-N/2,N/2,N);
yaxis = linspace(-N/4,3*N/4,N);
[xmesh,ymesh] = meshgrid(xaxis,yaxis);

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
numPixels = N*N;
newS = zeros(numStimuli,numPixels,'single');
for ii=1:numStimuli
    if mod(ii,2) == 1
        newS(ii,:) = S(floor(ii/2)+1,:);
    elseif mod(ii,2) == 0
        newS(ii,:) = 255-S(ii/2,:); %255-S(ii/2,:)
    end
end
newS = Grey-newS;
clear S;

strobeData = tsevs{1,strobeStart};

historyParams = 1;
gaborParams = 8;
numParameters = historyParams+gaborParams;

twopi = 2*pi;
Bounds = zeros(numParameters,2);
Bounds(1:historyParams,:) = ones(historyParams,2).*[-250,50];
Bounds(historyParams+1,:) = [0,1]; % height of Gabor parameter
Bounds(historyParams+2,:) = [min(xaxis)-5,max(xaxis)+5]; % x center
Bounds(historyParams+3,:) = [min(yaxis)-5,max(yaxis)+5]; % y center
Bounds(historyParams+4,:) = [0,20]; % standard deviation x
Bounds(historyParams+5,:) = [0,20]; % standard deviation y
Bounds(historyParams+6,:) = [0,10]; % spatial frequency
Bounds(historyParams+7,:) = [0,twopi]; %  orientation theta
Bounds(historyParams+8,:) = [0,twopi]; % phase shift phi

% convert timeStamps to point process
totalTime = strobeData(end)+2;
totalMillisecs = round(totalTime*1000);

stimTimes = round(strobeData.*1000);
pointProcessStimTimes = zeros(totalMillisecs,1);
stimOnset = 60;
stimOffset = 120;
stimulusTime = stimOffset-stimOnset+1;
for kk=1:numStimuli
    pointProcessStimTimes((stimTimes(kk)+stimOnset):(stimTimes(kk)+stimOffset)) = kk;
end

for zz=1:numChans
    numNeurons = Neurons(zz);
    for yy=1:numNeurons
        % ORGANIZE DATA
        
        Y = zeros(totalMillisecs,1);
        historyDesign = zeros(totalMillisecs,historyParams);
        
        % convert spike times to point process
        temp = round(allts{zz,yy+1}.*1000);
        for kk=1:length(temp)
            Y(temp(kk)) = 1;
        end
        
        % get history dependence (past 25ms)
        historyDesign(:,1) = ones(totalMillisecs,1,'single');
%         for kk=2:historyParams
%             temp = Y;shift = zeros(kk-1,1);
%             history = [shift;temp];
%             historyDesign(:,kk) = history(1:(end-(kk-1)));
%         end
        
        
        % RUN MCMC
        h = ones(numParameters,1)./1000;
        % repeat gradient ascent from a number of different starting
        %  positions
        numIter = 5e5;burnIn = 1e5;
        parameterVec = zeros(numParameters,numIter);
        
        parameterVec(1,1) = normrnd(-50,50);
        for ii=2:historyParams
            parameterVec(ii,1) = normrnd(0,10);
        end
        for ii=1:gaborParams
            parameterVec(historyParams+ii,1) = unifrnd(Bounds(historyParams+ii,1),Bounds(historyParams+ii,2));
        end
        
        parameterVec(:,1) = max(Bounds(1:numParameters,1),min(parameterVec(:,1),Bounds(1:numParameters,2)));
        
        % currently requires 0.045 seconds ... 0.0864 seconds
        %  just to get it to run in 1 day
        likelihoodXprev = GetLikelihood(parameterVec(:,1));

        
        sigma = zeros(numParameters,1);
        sigma(1:historyParams) = ones(historyParams,1).*10;
        sigma(historyParams+1) = 0.1;
        sigma(historyParams+2) = 10;
        sigma(historyParams+3) = 10;
        sigma(historyParams+4) = 10;
        sigma(historyParams+5) = 10;
        sigma(historyParams+6) = 0.1;
        sigma(historyParams+7) = 0.1;
        sigma(historyParams+8) = 0.1;
        
        rejectRate = 0;
        for iter = 2:numIter
            xStar = parameterVec(:,iter)+normrnd(0,sigma,[numParameters,1]);
            xStar(historyParams+7:historyParams+8) = mod(xStar(historyParams+7:historyParams+8),twopi);
            temp = xStar(1:historyParams+6)-Bounds(1:historyParams+6,1);
            
            if sum(temp<0) > 0 % outside bounds
                parameterVec(:,iter) = parameterVec(:,iter-1);
                rejectRate = rejectRate+1;
            else
                likelihoodXstar = GetLikelihood(xStar);
                %         A = min(1,(likelihoodXstar*priorXstar)./(likelihoodXprev*priorXprev));
                %             A = min(1,likelihoodXprev/likelihoodXstar);
                logA = likelihoodXstar-likelihoodXprev;
                if log(rand) < logA
                    parameterVec(:,iter) = xStar;
                    likelihoodXprev = likelihoodXstar;
                else
                    parameterVec(:,iter) = parameterVec(:,iter-1);
                    rejectRate = rejectRate+1;
                end
            end
                
        end
        
        fprintf('Rejection Rate: %3.2f',rejectRate/numIter);
        
        PosteriorSamples = parameterVec(:,burnIn:40:numIter);
        
        PosteriorMean = zeros(numParameters,1);
        PosteriorInterval = zeros(numParameters,2);
        alpha = 0.05;
        
        numColumns = 4;
        numRows = floor(numParameters/numColumns)+mod(numParameters,numColumns);
        for ii=1:numParameters
           subplot(numRows,numColumns,ii);
           histogram(PosteriorSamples(ii,:));
           title('Marginal Posterior, Parameter %d',ii);
           PosteriorMean(ii) = mean(PosteriorSamples(ii,:));
           PosteriorInterval(ii,:) = quantile(PosteriorSamples(ii,:),[alpha/2,1-alpha/2]);
        end
       
    end
end


end

function [loglikelihood] = GetLikelihood(parameterVec)
global newS pointProcessStimTimes historyParams xmesh ymesh numPixels ...
    historyDesign Y stimulusTime;

baseMu = exp(historyDesign(1,:)*parameterVec(1:historyParams));
loglikelihood = sum(Y(pointProcessStimTimes==0).*log(baseMu)-baseMu);



gaborFilter = parameterVec(historyParams+1)*exp(-(xmesh-parameterVec(historyParams+2)).^2 ...
        ./(2*parameterVec(historyParams+4).^2)-(ymesh-parameterVec(historyParams+3)).^2 ...
        ./(2*parameterVec(historyParams+5).^2)).*sin(2*pi*parameterVec(historyParams+6)*...
        ((xmesh-parameterVec(historyParams+2))*cos(parameterVec(historyParams+7)-pi/2)+...
        (ymesh-parameterVec(historyParams+3))*sin(parameterVec(historyParams+7)-pi/2))+...
        parameterVec(historyParams+8));
gaborFilter = gaborFilter(:);

zz = find(pointProcessStimTimes~=0);
loglikelihoodTwo = zeros(length(zz),1);

for kk=1:stimulusTime:length(zz)
%     onScreenStim = squeeze(newS(pointProcessStimTimes(zz(kk)),:,:));
%     if mod(kk,stimulusTime) == 1
%         filterOutput = sum(newS(pointProcessStimTimes(zz(kk),:))'.*gaborFilter)./numPixels;
%     end
    filterOutput = sum(newS(pointProcessStimTimes(zz(kk),:))'.*gaborFilter)./numPixels;
    
    loglikelihoodTwo(kk:kk+stimulusTime-1) = Y(zz(kk:kk+stimulusTime-1))*log(baseMu*exp(filterOutput))-baseMu*exp(filterOutput);
end

% for kk=1:numStimuli
%    filterOutput = sum(newS(kk,:)'.*gaborFilter)./numPixels;
%    loglikelihood(pointProcessStimTimes==kk) = Y(pointProcessStimTimes==kk).*log(baseMu*exp(filterOutput))-...
%        baseMu*exp(filterOutput);
% end

loglikelihood = loglikelihood+sum(loglikelihoodTwo);
end

