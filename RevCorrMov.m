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

forMovStrobed = svStrobed(svStrobed>0 & svStrobed<254);

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

temp = strobeData(svStrobed>0 & svStrobed<254);
stimTimes = round(temp.*timeMultiplier);

pointProcessSpikes = zeros(totalMillisecs,totalUnits);

for ii=1:totalUnits
   spikeTimes = max(1,round(allts{ii}.*timeMultiplier));
   for jj=1:length(spikeTimes)
      pointProcessSpikes(spikeTimes(jj),ii) = pointProcessSpikes(spikeTimes(jj),ii)+1;
   end
end


movieLen = movie_FrameRate*movieTime_Seconds;

temp = randperm(numMoviesToDisplay);
divide = round(0.7*numMoviesToDisplay);
train = temp(1:divide);test = temp(divide+1:end);clear temp divide;

load(sprintf('%sMovie1.mat',movieType),'DIM');

xaxis = linspace(-round(screenPix_to_effPix*DIM(2)/2)+1,...
    round(screenPix_to_effPix*DIM(2)/2),DIM(2));
yaxis = linspace(round(3*screenPix_to_effPix*DIM(1)/4),...
    -round(screenPix_to_effPix*DIM(1)/4)+1,DIM(1));
taxis = linspace(500,50,28);


movieIndices = zeros(movieLen*numMoviesToDisplay,length(taxis));
movieFrames = zeros(movieLen*numMoviesToDisplay,DIM(1)*DIM(2),'uint8');
reducedSpikeCount = zeros(movieLen*numMoviesToDisplay,totalUnits);
reducedMov = zeros(movieLen*numMoviesToDisplay,1);

count = 1;
for ii=1:numMoviesToDisplay
    index = movieNums(ii);
    fileName = sprintf('%sMovie%d.mat',movieType,index);
    load(fileName,'S');
    currentMovStrobes = forMovStrobed(forMovStrobed==index);
    for jj=1:movieLen
        temp = S(:,:,jj);
        movieFrames(count,:) = temp(:);
        temp2 = count-length(taxis):count-1;
        movieIndices(count,:) = max(temp2,0);
        onsetTime = stimTimes(currentMovStrobes(jj));
        offsetTime = onsetTime+round(1/movie_FrameRate);
        for kk=1:totalUnits
            reducedSpikeCount(count,kk) = sum(pointProcessSpikes(onsetTime:offsetTime,kk));
            reducedMov(count) = sum(movement(onsetTime:offsetTime));
        end
        count = count+1;
    end
end

movieFrames = [Grey.*ones(1,DIM(1)*DIM(2),'uint8');movieFrames];
movieIndices = movieIndices+1;

clearvars -except EphysFileName totalUnits numStimuli ...
    reducedSpikeCount DIM unbiasedS allts strobeData xaxis yaxis taxis...
    totalMillisecs reducedMov movement AnimalName numChans movieIndices ...
    movieFrames pointProcessSpikes svStrobed forMovStrobed stimTimes;

% REGULARIZED PSEUDO-INVERSE SOLUTION, CREATE CONVOLUTION MATRIX L
[centerPositions,rfInds,newDims] = GetRetinoMap(AnimalName,xaxis,yaxis,numChans);

warning('off','all');

numLambda = 50;
loglambda = logspace(1,7,numLambda);
F = cell(totalUnits,1);
STA = zeros(totalUnits,DIM(1)*DIM(2));

bestLambda = zeros(totalUnits,1);
heldOutDeviance = zeros(totalUnits,5);
heldOutExplainedVariance = zeros(totalUnits,1);
sigmoidNonlin = zeros(totalUnits,4);
visualResponsiveness = zeros(totalUnits,2);
for zz=1:totalUnits
   fprintf('Running unit %d ...\n',zz);
   spikeTrain = squeeze(reducedSpikeCount(:,zz));
   
   baseRate = sum(pointProcessSpikes
   
   y = [sum(spikeTrain(:,51:150),2);sum(spikeTrain(:,901:1000),2)];
   design1 = ones(2*numStimuli,1);
   design2 = [[ones(numStimuli,1);zeros(numStimuli,1)],...
       [zeros(numStimuli,1);ones(numStimuli,1)]];
   [~,dev1,~] = glmfit(design1,y,'poisson','constant','off');
   [~,dev2,~] = glmfit(design2,y,'poisson','constant','off');
   dfDiff = 1;
   pVal = chi2cdf(dev1-dev2,dfDiff,'upper');
   visualResponsiveness(zz,1) = pVal;
   visualResponsiveness(zz,2) = pVal<0.01;
   
   if visualResponsiveness(zz,2) == 1
       spikeTrain = sum(spikeTrain(:,50:150),2);
       
       % PCA solution
       %    [V,D] = eig(cov(S));
       %    eigenvals = diag(D);
       %    start = find(eigenvals>10,1,'first');
       %    q = size(eigenvals(start:end),1);
       %    meanEig = mean(eigenvals(1:start-1));
       %    W = V(:,start:end)*sqrtm(D(start:end,start:end)-meanEig.*eye(q));
       %    W = fliplr(W);
       %    Winv = pinv(W);
       %    x = zeros(q,numStimuli);
       %    for ii=1:numStimuli
       %       x(:,ii) = Winv*(S(ii,:)'-127);
       %    end
       %    x = x';
       %    fhat = x(train,:)\spikeTrain(train); % then proceed as usual, with no smoothing
       
       % GET SMALLER SECTION OF IMAGE DISPLAYED ON THE SCREEN
       fullSize = newDims(unitChannel(zz),1)*newDims(unitChannel(zz),2);
       L = sparse(fullSize,fullSize);
       
       %operator = [0,-1,0;-1,4,-1;0,-1,0];
       bigCount = 1;
       for jj=1:newDims(unitChannel(zz),2)
           for ii=1:newDims(unitChannel(zz),1)

               tempMat = zeros(newDims(unitChannel(zz),1),newDims(unitChannel(zz),2));
               
               if ii==1 && jj==1
                   tempMat(ii,jj) = 4;
                   tempMat(ii+1,jj) = -4/2;
                   tempMat(ii,jj+1) = -4/2;
               elseif ii==newDims(unitChannel(zz),1) && jj==1
                   tempMat(ii,jj) = 4;
                   tempMat(ii-1,jj) = -4/2;
                   tempMat(ii,jj+1) = -4/2;
               elseif ii==1 && jj==newDims(unitChannel(zz),2)
                   tempMat(ii,jj) = 4;
                   tempMat(ii,jj-1) = -4/2;
                   tempMat(ii+1,jj) = -4/2;
               elseif ii == newDims(unitChannel(zz),1) && jj == newDims(unitChannel(zz),2)
                   tempMat(ii,jj) = 4;
                   tempMat(ii-1,jj) = -4/2;
                   tempMat(ii,jj-1) = -4/2;
               elseif ii==1
                   tempMat(ii,jj) = 4;
                   tempMat(ii,jj-1) = -4/3;
                   tempMat(ii+1,jj) = -4/3;
                   tempMat(ii,jj+1) = -4/3;
               elseif jj==1
                   tempMat(ii,jj) = 4;
                   tempMat(ii-1,jj) = -4/3;
                   tempMat(ii,jj+1) = -4/3;
                   tempMat(ii+1,jj) = -4/3;
               elseif ii==newDims(unitChannel(zz),1)
                   tempMat(ii,jj) = 4;
                   tempMat(ii-1,jj) = -4/3;
                   tempMat(ii,jj-1) = -4/3;
                   tempMat(ii,jj+1) = -4/3;
               elseif jj==newDims(unitChannel(zz),2)
                   tempMat(ii,jj) = 4;
                   tempMat(ii-1,jj) = -4/3;
                   tempMat(ii+1,jj) = -4/3;
                   tempMat(ii,jj-1) = -4/3;
               else
                   tempMat(ii,jj) = 4;
                   tempMat(ii-1,jj) = -1;
                   tempMat(ii+1,jj) = -1;
                   tempMat(ii,jj+1) = -1;
                   tempMat(ii,jj-1) = -1;
               end
               L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
               %         imagesc(tempMat);caxis([-1 4]);pause(0.1);
           end
       end
       
       r = [spikeTrain(train);sparse(fullSize,1)];
       
       tempF = zeros(numLambda,fullSize);
       tempSigmoid = zeros(numLambda,4);
       tempDev = zeros(numLambda,4);
       %    allOnesTrain = ones(length(train),1);
       %    allOnesTest = ones(length(test),1);
       for jj=1:numLambda
           % calculate regularized pseudoinverse solution
           constraints = [S(train,rfInds{unitChannel(zz)});loglambda(jj).*L];
           fhat = constraints\r;fhat = full(fhat);
           
           % fit sigmoid nonlinearity
           result = S(train,rfInds{unitChannel(zz)})*fhat;
           myFun = @(x) x(1)./(1+exp(-(result-x(2)).*x(3)))+x(4)-spikeTrain(train);
           x0 = [5,0.5,10,baseRate];
           lb = [0,0,0,0];ub = [5e2,Inf,Inf,Inf];
           %       options.Algorithm = 'levenberg-marquardt';
           %options.MaxFunctionEvaluations = 1000;
           options = optimoptions(@lsqnonlin,'MaxIterations',1000);
           options.Display = 'off';
           sigmoidParams = lsqnonlin(myFun,x0,lb,ub,options);
           
           temp = sigmoidParams(1)./(1+exp(-(result-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
           initialDev = spikeTrain(train).*log(spikeTrain(train)./temp)-(spikeTrain(train)-temp);
           initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
           sigmoidParams = FitSigmoid(sigmoidParams,2*sum(initialDev),spikeTrain(train),result);
           
           % save parameters
           tempF(jj,:) = fhat;
           tempSigmoid(jj,:) = sigmoidParams;
           
           % calculate held-out Poisson deviance
           temp = S(test,rfInds{unitChannel(zz)})*fhat;
           temp = sigmoidParams(1)./(1+exp(-(temp-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
           
           initialDev = spikeTrain(test).*log(spikeTrain(test)./temp)-(spikeTrain(test)-temp);
           initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
           tempDev(jj,1) = 2*sum(initialDev);
           
           % rotate the receptive field and recalculate the held-out deviance
           rf = reshape(fhat,[newDims(unitChannel(zz),1),newDims(unitChannel(zz),2)]);
           for kk=2:4
               rf = imrotate(rf,90);
               tempfhat = rf(:);
               temp = S(test,rfInds{unitChannel(zz)})*tempfhat;
               temp = sigmoidParams(1)./(1+exp(-(temp-sigmoidParams(2)).*sigmoidParams(3)))+sigmoidParams(4);
               
               initialDev = spikeTrain(test).*log(spikeTrain(test)./temp)-(spikeTrain(test)-temp);
               initialDev(isnan(initialDev) | isinf(initialDev)) = temp(isnan(initialDev) | isinf(initialDev));
               tempDev(jj,kk) = 2*sum(initialDev);
           end
       end
       
       [~,bestMap] = min(tempDev(:,1));
       fprintf('Best Map: %d\n',bestMap);
       F{zz} = tempF(bestMap,:);
       bestLambda(zz) = loglambda(bestMap);
       heldOutDeviance(zz,1:4) = tempDev(bestMap,:);
       sigmoidNonlin(zz,:) = tempSigmoid(bestMap,:);
       
       % get held-out deviance of model with single regressor
       [~,dev,~] = glmfit(ones(length(test),1),spikeTrain(test),'poisson','constant','off');
       heldOutDeviance(zz,5) = dev;
       
       % get held-out explained variance
       %  see R-squared measures for Count Data Regression Models
       %   A. Colin Cameron & Frank A.G. Windmeijer April 1995
       %   Journal of Business and Economic Statistics
       heldOutExplainedVariance(zz,1) = 1-tempDev(bestMap,1)/dev;
       fprintf('Fraction of Explained Variance: %3.2f\n\n',1-tempDev(bestMap,1)/dev);
       
       if heldOutExplainedVariance(zz,1) >= 0.05
           figure();imagesc(reshape(fhat,[newDims(unitChannel(zz),1),newDims(unitChannel(zz),2)]));
           title(sprintf('%s',EphysFileName(1:end-9)));
       end
       
       totalSpikes = sum(spikeTrain);
       for ii=1:numStimuli
           if spikeTrain(ii) > 0
               STA(zz,:) = STA(zz,:)+S(ii,:).*(spikeTrain(ii)/totalSpikes);
           end
       end
   end
end

% save the results
fileName = strcat(EphysFileName(1:end-9),'-PseudoInvResults.mat');
save(fileName,'F','totalUnits','bestLambda',...
    'reducedSpikeCount','DIM','S','movement',...
    'xaxis','yaxis','reducedMov','allts','strobeData','totalMillisecs',...
    'svStrobed','heldOutDeviance','numStimuli','sigmoidNonlin',...
    'heldOutExplainedVariance','STA','beta','centerPositions','rfInds',...
    'newDims','unitChannel','visualResponsiveness');

% REVERSE CORRELATION SOLUTION
% for ii=1:numChans
%     numNeurons = Neurons(ii,1);
%     figure();plotRows = ceil(numNeurons/2);
%     for jj=1:numNeurons
%         r = squeeze(Response(ii,jj,:));
% % %         r = newS*f
%         fhat = newS\r;
%         F(ii,jj,:) = fhat;
%         subplot(plotRows,2,jj);imagesc(reshape(fhat,[N,N]));
%         title(sprintf('RevCorr: Chan %d, Unit %d, Animal %d',ii,jj,AnimalName));
%     end
% end
% Img = reshape(S(tt,:),[minPix/screenPix_to_effPix,minPix/screenPix_to_effPix]);
% Img = kron(double(Img),ones(screenPix_to_effPix));

end

function [sigmoidParams] = FitSigmoid(sigmoidParams,initialDev,spikeTrain,result)
N = 1e5;
numParams = length(sigmoidParams);
paramVec = zeros(numParams,N);
devVec = zeros(N,1);

variance = 0.01;

paramVec(:,1) = sigmoidParams';
devVec(1) = initialDev;

W = diag(ones(numParams,1));
loglambda = log(2.38^2)*ones(numParams,1);
optimalAccept = 0.44;
for ii=2:N
    index = random('Discrete Uniform',numParams);
    
    pStar = paramVec(:,ii-1)+W(:,index).*normrnd(0,sqrt(exp(loglambda(index))*variance));
    
    if sum(pStar<=0) == 0
        temp = pStar(1)./(1+exp(-(result-pStar(2)).*pStar(3)))+pStar(4);
        
        pDev = spikeTrain.*log(spikeTrain./temp)-(spikeTrain-temp);
        pDev(isnan(pDev) | isinf(pDev)) = temp(isnan(pDev) | isinf(pDev));
        
        logA = devVec(ii-1) - 2*sum(pDev);
        if log(rand) < logA
           paramVec(:,ii) = pStar;
           devVec(ii) = 2*sum(pDev);
        else
            paramVec(:,ii) = paramVec(:,ii-1);
            devVec(ii) = devVec(ii-1);
            
        end
        loglambda(index) = loglambda(index)+0.01*(exp(min(0,logA))-optimalAccept);
    else
       loglambda(index) = loglambda(index)+0.01*(-optimalAccept);
       paramVec(:,ii) = paramVec(:,ii-1);
       devVec(ii) = devVec(ii-1);      
    end
end

[~,ind] = min(devVec);
sigmoidParams = paramVec(:,ind);

end

function [centerPositions,rfInds,newDims] = GetRetinoMap(AnimalName,xaxis,yaxis,numChans)
cd ~/CloudStation/ByronExp/Retino/
fileName = strcat('RetinoMapBayes*',num2str(AnimalName),'.mat');
files = dir(fileName);

if isempty(files)==1
    centerPositions = zeros(numChans,2);
    rfInds = cell(numChans,1);
    newDims = zeros(numChans,2);
    for ii=1:numChans
        tempIm = ones(length(yaxis),length(xaxis));
        rfInds{ii} = find(tempIm==1);
        newDims(ii,1) = length(yaxis);
        newDims(ii,2) = length(xaxis);
    end
else
    load(files(end).name);
    
    [numChans,~,numSamples] = size(posteriorSample);
    
    x = 1:w_pixels;y=1:h_pixels;
    [X,Y] = meshgrid(x,y);
    
    xpix = linspace(-round(w_pixels/2)+1,...
        round(w_pixels/2),w_pixels);
    ypix = linspace(round(3*h_pixels/4),...
        -round(h_pixels/4)+1,h_pixels);
    
    % we want final RF to be ~1200 by 1200 screen pixels, about 50 degrees
    %  of visual arc on a side
    
    centerPositions = zeros(numChans,2);
    rfInds = cell(numChans,1);
    newDims = zeros(numChans,2);
    for ii=1:numChans
        finalIm = zeros(length(y),length(x));
        samples = squeeze(posteriorSample(ii,:,:));
        N = 1000;
        for ll=1:N
            index = random('Discrete Uniform',numSamples);
            parameterVec = samples(:,index);
            b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(6)];
            distX = X-parameterVec(2);distY = Y-parameterVec(3);
            finalIm = finalIm+b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
                (distY.^2)./(2*b(3)*b(3)))+b(4);
            %         for jj=1:length(x)
            %             for kk=1:length(y)
            %                 distX = x(jj)-parameterVec(2);
            %                 distY = y(kk)-parameterVec(3);
            %
            %                 finalIm(jj,kk) = finalIm(jj,kk)+b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
            %                     (distY.^2)./(2*b(3)*b(3)))+b(4);
            %             end
            %         end
        end
        finalIm = finalIm./N;
        [~,maxInd] = max(finalIm(:));
        [row,col] = ind2sub(size(finalIm),maxInd);
        centerPositions(ii,1) = col;
        centerPositions(ii,2) = row;
        
        centerX = xpix(col);centerY = ypix(length(ypix)-row+1);
        xLow = centerX-600;xHigh = centerX+600;
        yLow = centerY-600;yHigh = centerY+600;
        
        [~,xLow] = min(abs(xLow-xaxis));
        [~,xHigh] = min(abs(xHigh-xaxis));
        [~,yLow] = min(abs(yLow-yaxis));
        [~,yHigh] = min(abs(yHigh-yaxis));
        tempIm = zeros(length(yaxis),length(xaxis));
        tempIm(yHigh:yLow,xLow:xHigh) = 1;
        rfInds{ii} = find(tempIm==1);
        newDims(ii,1) = yLow-yHigh+1;
        newDims(ii,2) = xHigh-xLow+1;
    end
end

cd ~/CloudStation/ByronExp/NoiseRetino/

end