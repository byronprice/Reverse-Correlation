function [] = RevCorrUnits(AnimalName,Date,NoiseType)
%RevCorrUnits.m
%   %   Analysis of single unit recording data in response to a series of white
%   or pink noise stimuli (see Noise_RevCorr.m for Psychtoolbox
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
%       Chans - channel numbers, input as [6,8], defaults to 6 and 8
%    
%      Input through saved file 'NoiseStimData_AnimalName.mat'
%      S - number of stimuli-by-number of pixels matrix of white noise
%        stimuli used by the function Noise_RevCorr.m
%       R - matrix of size number of frames-by-number of neurons that
%        contains neuronal data which has already been processed into a
%        point process and which will be converted to a new matrix of responses, 
%        number of stimuli-by-number of neurons matrix of neuronal
%        numbers obtained from calcium imaging (count threshold
%        crossings as spikes ... we will want to know the number of spikes
%        in response to the onset and offset of each stimuli.

%OUTPUT: F - number of pixels-by-number of neurons matrix representing the
%         receptive field for each neuron recorded from in R
%            if N = sqrt(size(F,1));
%            then the receptive field for a given neuron can be 
%            visualized by:
%            figure();imagesc(reshape(F(:,ii),[N,N]);      
%    
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/04/05
% By: Byron Price

% CREATE GIANT CONVOLUTION MATRIX TO TAKE DISCRETE
%  LAPLACIAN OF THE SPATIAL RECEPTIVE FIELD f for each
%  neuron

% assume we sample at some frequency, Fs, and have a point process
%  representation of the spiking behavior of a given neuron, we 
%  have the set of stimuli used in the input matrix S and the 
%  timeStamps for when each was displayed
%  If neuron i fired after the onset of stimulus j, then we'll take 
%  the number of spikes in a certain window after the onset of that
%  stimulus as the response, with one caveat: we'll take the spikes 
%  at the offset of the stimulus and count those as negative responses.
%
%   Note: we could record from these neurons for a period with no visual
%   stimulus present and get a background firing rate for each and then
%   normalize the response by that background firing rate, though this may
%   not add anything

% assume we want the period from 50msec to 120msec after stimulus
%  onset and that the calcium imaging recording started just in time 
%  with the first stimulus (which was a grey screen)

% read in the .plx file
beta = 0;

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'-sort.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName(1:end-4),'.plx'));pause(1);
end

StimulusFileName = strcat('NoiseStim',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
load(EphysFileName,'allts','allad','adfreq','tsevs','svStrobed','nunits1');
load(StimulusFileName,'numStimuli','minPix','maxPix','conv_factor','screenPix_to_effPix','effectivePixels',...
    'DistToScreen','beta','S');

% generate theoretical stimulus power spectrum
if length(effectivePixels) == 1
    N = sqrt(effectivePixels);
    DIM = [N,N];
else
   temp = effectivePixels;
   effectivePixels(1) = temp(2);effectivePixels(2) = temp(1);
   DIM = effectivePixels;
end

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
% Reproduce these frequencies along ever row
U = repmat(u,1,DIM(2)); 
% v is the set of frequencies along the second dimension.  For a square
% region it will be the transpose of u
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]/DIM(2);
% Reproduce these frequencies along ever column
V = repmat(v,DIM(1),1);

% [U,V] = meshgrid(u,v); U = U'; V = V';
% Generate the power spectrum
S_f = (U.^2 + V.^2).^(beta/2);

% Set any infinities to zero
S_f(S_f==inf) = 0;

tempSf = S_f; tempSf(tempSf==0) = 1;

% correct power spectrum of the images to remove bias in RF estimate
%  I control the image generation so we know the theoretical power spectrum
parfor ii=1:numStimuli
    tempS = reshape(S(ii,:),[DIM(1),DIM(2)]);
    tempFFT = fft2(double(tempS));
    correctedIm = ifft2(tempFFT./tempSf);
    S(ii,:) = real(correctedIm(:));
end

temp = ~cellfun(@isempty,allts);
Chans = find(sum(temp,1));numChans = length(Chans);
totalUnits = sum(sum(temp));

% tsevs are the strobed times of stimulus onset, then offset
%  Onset at tsevs{1,33}(2), offset at tsevs{1,33}(3), onset at
%  tsevs{1,33}(4), offset at 5, etc.

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

% ASSUME THAT A RESPONSE TO A STIMULUS OFFSET IS THE SAME AS A RESPONSE TO
%  THE NEGATIVE OF THAT IMAGE, image created with values from 0 to 255
% Grey = 127;
numStimuli = numStimuli*2;
newS = zeros(numStimuli,size(S,2),'single');
for ii=1:numStimuli
    if mod(ii,2) == 1
        newS(ii,:) = S(floor(ii/2)+1,:);
    elseif mod(ii,2) == 0
        newS(ii,:) = 255-S(ii/2,:);
    end
end
%newS = newS-Grey;
clear S;

% newS = S;

strobeStart = 33;
strobeData = tsevs{1,strobeStart};

stimLen = 0.07;

nonEmptyAD = ~cellfun(@isempty,allad);
LFP_Movement = allad{nonEmptyAD};

if iscell(LFP_Movement)==1
    totalTime = length(LFP_Movement{1})./adfreq;
else
    totalTime = length(LFP_Movement)./adfreq;
end
% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = zeros(totalUnits,numStimuli);
baseRate = zeros(totalUnits,1);
for ii=1:totalUnits
    stimStart = 0.05;%stimStartTimes(ii,jj+1);
    baseRate(ii) = length(allts{ii})./totalTime;
    %             display(baseRate(ii,jj));
    %         figure();plot(0,0);axis([0 totalTime -10 10]);hold on;
    for kk=1:numStimuli
        stimOnset = strobeData(kk);
        high = find(allts{ii} > (stimOnset+stimStart));
        low = find(allts{ii} < (stimOnset+stimStart+stimLen));
        temp = intersect(low,high);
        
        Response(ii,kk) = length(temp)./stimLen;
%         if mod(kk,2) == 1
%             Response(ii,kk) = (length(temp)./stimLen)./baseRate(ii,jj+1);%-1
%         elseif mod(kk,2) == 0
%             Response(ii,kk) = -(length(temp)./stimLen)./baseRate(ii,jj+1); %+1
%         end
        
    end
end


% REGULARIZED PSEUDO-INVERSE SOLUTION
horzDegrees = atand((screenPix_to_effPix*DIM(2)*conv_factor/10)/DistToScreen);
vertDegrees = atand((screenPix_to_effPix*DIM(1)*conv_factor/10)/DistToScreen);
xaxis = linspace(-horzDegrees/2,horzDegrees/2,DIM(2));
yaxis = linspace(-vertDegrees/4,3*vertDegrees/4,DIM(1));

fullSize = DIM(1)*DIM(2);
numParams = fullSize+1+1+1;

Design = zeros(numStimuli,numParams-2);
Design(:,1) = ones(numStimuli,1);
Design(:,2:end) = newS;clear newS;
logPoissonPDF = @(y,mu) y.*mu-exp(mu);
logGammaPDF = @(x,a,b) -a*log(b)-log(gamma(a))+(a-1).*log(x)-x./b;

for unit = 1:totalUnits
    y = Response(unit,:)';
    baseFiring = baseRate(unit);
    numIter = 7e5;burnIn = 2e5;skipRate = 500;
    params = zeros(numParams,(numIter-burnIn)/skipRate);
    posteriorProb = zeros((numIter-burnIn)/skipRate,1);
    prevParams = zeros(numParams,1);
    
    priorMu = log(baseFiring);
    
    alpha = 1;beta = 1;
    abprior1=1e-3;abprior2=1e3;
    
    
    prevPosterior = -Inf;
    while prevPosterior < -1e10 || isnan(prevPosterior)==1
        prevParams(1:end-2,1) = mvnrnd([priorMu;zeros(fullSize,1)],eye(numParams-2))';
        prevParams(end-1,1) = log(alpha);
        prevParams(end,1) = log(beta);
        
        
        tempMu = Design*prevParams(1:end-2,1);
        prevLogLikelihood = sum(logPoissonPDF(y,tempMu));
        
        smoothPrior = del2(reshape(prevParams(2:end-2,1),[DIM(1),DIM(2)]));
        prevLogPrior = ((fullSize)/2)*log(alpha)-0.5*alpha*(prevParams(2:end-2,1)'*prevParams(2:end-2,1))+...
            0.5*log(alpha)-0.5*alpha*(prevParams(1,1)-priorMu)^2+...
            ((fullSize)/2)*log(beta)-0.5*beta*(smoothPrior(:)'*smoothPrior(:))+...
            sum(logGammaPDF([alpha,beta],abprior1,abprior2));
        
        prevPosterior = prevLogLikelihood+prevLogPrior;
    end
    
    
%     FOR AUTOMATIC CREATION OF UPDATE MATRIX
    updateParam = logspace(-0.3,-3,burnIn);
    loglambda = log(2.38^2/numParams);
    updateMu = zeros(numParams,1);
    updateMu(1:end-2) = mvnrnd([priorMu;zeros(fullSize,1)],eye(numParams-2))';
    updateMu(end-1) = log(1.5);updateMu(end) = log(1.5);
    optimalAccept = 0.234;
    
    sigma = eye(numParams);halfSigma = cholcov(sigma);
    identity = eye(numParams);
    
    proposalMu = zeros(numParams,1)';
    for ii=2:burnIn
        tic;
        pStar = prevParams+mvnrnd(proposalMu,exp(loglambda).*sigma)';
        
        if sum(pStar(end-1:end)<=-30) == 0
            tempMu = Design*pStar(1:end-2);
            pStarLogLikelihood = sum(logPoissonPDF(y,tempMu));
            smoothPrior = del2(reshape(pStar(2:end-2),[DIM(1),DIM(2)]));
            
            beta = exp(pStar(end));alpha = exp(pStar(end-1));
            pStarLogPrior = ((fullSize)/2)*log(alpha)-0.5*alpha*(pStar(2:end-2)'*pStar(2:end-2))+...
                0.5*log(alpha)-0.5*alpha*(pStar(1)-priorMu)^2+...
                ((fullSize)/2)*log(beta)-0.5*beta*(smoothPrior(:)'*smoothPrior(:))+...
                sum(logGammaPDF([alpha,beta],abprior1,abprior2));
            logA = (pStarLogLikelihood+pStarLogPrior)-prevPosterior;

            if log(rand) < logA
                prevParams = pStar;
                prevPosterior = pStarLogLikelihood+pStarLogPrior;
            end
            
            meanSubtract = prevParams-updateMu;
            updateMu = updateMu+updateParam(ii).*meanSubtract;
            halfSigma = halfSigma+updateParam(ii).*(triu((halfSigma^-1)*(halfSigma'*halfSigma+meanSubtract*...
                meanSubtract')*((halfSigma^-1)')-identity)-halfSigma);
            sigma = halfSigma'*halfSigma;
            
            loglambda = loglambda+updateParam(ii).*(exp(min(0,logA))-optimalAccept);
        else
            loglambda = loglambda+updateParam(ii).*(-optimalAccept);
        end
        toc;
    end
    
    [V,D] = eig(sigma);
    W = V*sqrtm(D);
    eigenvals = diag(W'*W);
    
    tempW = [];
    tempEigs = [];
    for jj=1:numParams
        if eigenvals(jj) > 1e-6
            tempW = [tempW,W(:,jj)];
            tempEigs = [tempEigs,eigenvals(jj)];
        end
    end
    
    W = fliplr(tempW);
    eigenvals = fliplr(tempEigs);
    q = length(eigenvals);
    p = ones(q,1)./q;
    updateParam = 1e-2;
    loglambda = ones(q,1).*loglambda;
    params(:,1) = prevParams;
    posteriorProb(1) = prevPosterior;
    
    count = 2;
    for ii=burnIn+1:numIter
        index = find(mnrnd(1,p)==1);
        lambda = loglambda(index);
        stdev = sqrt(exp(lambda).*eigenvals(index));
        pStar = prevParams+W(:,index)*normrnd(0,stdev);
        
        if sum(pStar(end-1:end)<=-30) == 0
            tempMu = Design*pStar(1:end-2);
            pStarLogLikelihood = sum(logPoissonPDF(y,tempMu));
            smoothPrior = del2(reshape(pStar(2:end-2),[DIM(1),DIM(2)]));
            
            beta = exp(pStar(end));alpha = exp(pStar(end-1));
            pStarLogPrior = ((fullSize)/2)*log(alpha)-0.5*alpha*(pStar(2:end-2)'*pStar(2:end-2))+...
                0.5*log(alpha)-0.5*alpha*(pStar(1)-priorMu)^2+...
                ((fullSize)/2)*log(beta)-0.5*beta*(smoothPrior(:)'*smoothPrior(:))+...
                sum(logGammaPDF([alpha,beta],abprior1,abprior2));
            logA = (pStarLogLikelihood+pStarLogPrior)-prevPosterior;

            if log(rand) < logA
                prevParams = pStar;
                prevPosterior = pStarLogLikelihood+pStarLogPrior;
            end
            
            lambda = lambda+updateParam.*(exp(min(0,logA))-optimalAccept);
        else
            lambda = lambda+updateParam.*(-optimalAccept);
        end
        loglambda(index) = lambda;
        if mod(count,skipRate) == 0
            params(:,count) = prevParams;
            posteriorProb(count) = prevPosterior;
        end
        count = count+1;
    end
    fileName = strcat('NoiseResults',NoiseType,num2str(Date),'_',num2str(AnimalName),num2str(unit),'.mat');
    save(fileName,'Design','params','DIM','numParams');

end

% FileName = strcat('NoiseResults',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
% save(FileName,'F','newS','Response','allts','bigLambda','numChans','nunits1',...
%     'beta','spaceExp','N','RMS','bestMaps','totalUnits');

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

% SIMULATE DATA and test algorithm
% numStimuli = 5000;
% N = 48;
% newS = randn([numStimuli,N*N]);
% gaborFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3)).*sin(2*pi*0.05.*(x.*cos(pi/4)+y.*sin(pi/4)));
% gaussFilter = @(x,y) exp(-x.^2./(2*3*3)-y.^2./(2*3*3));
% x = linspace(-20,20,N);y = linspace(-20,20,N);
% [X,Y] = meshgrid(x,y);
% gabor = gaborFilter(X,Y);
% gauss = gaussFilter(X,Y);
% r = zeros(numStimuli,1);
% for ii=1:numStimuli
%     temp = newS(ii,:);
%     temp = temp-min(temp);temp = (temp./max(temp)).*255;
%     temp = temp-(mean(temp)-127);
%     newS(ii,:) = temp-127;
%     tempIm = reshape(temp,[N,N]);
%     gaussOutput = 1;%sum(sum(conv2(tempIm,gauss)));
%     gaborOutput = sum(newS(ii,:)'.*gabor(:));
%     lambda = exp((gaborOutput/gaussOutput)./(N*N));
%     r(ii) = poissrnd(lambda);
% end
% 
% L = zeros(N*N,N*N,'single');
% 
% %operator = [0,-1,0;-1,4,-1;0,-1,0];
% bigCount = 1;
% for jj=1:N
%     for ii=1:N
%         tempMat = zeros(N,N);
%         tempMat(ii,jj) = 4;
%         if ii > 1
%             tempMat(ii-1,jj) = -1;
%         end
%         if ii < N
%             tempMat(ii+1,jj) = -1;
%         end
%         if jj > 1
%             tempMat(ii,jj-1) = -1;
%         end
%         if jj < N
%             tempMat(ii,jj+1) = -1;
%         end
%         L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
%     end
% end
% bigLambda = [0,5e1,1e2,1e3,1e4,5e4,1e5,1e6,1e7,1e8];
% RMS = zeros(length(bigLambda),1);
% parfor ii=1:length(bigLambda)
%     A = [newS;bigLambda(ii).*L];
%     constraints = [r;zeros(N*N,1)];
%     fhat = pinv(A)*constraints;
%     RMS(ii) = (N*N)^(-0.5)*norm(r-newS*fhat);
% end
% rmsDiff = diff(RMS);lambdaDiff = diff(bigLambda)';
% deltaRMSdeltaLambda = rmsDiff./lambdaDiff;
% [maxVal,index] = max(abs(deltaRMSdeltaLambda));
% onepercent = 0.01*maxVal;
% firstBelow = find(abs(deltaRMSdeltaLambda(index:end))<onepercent,1);
% bestMap = index+firstBelow-1;
% 
% A = [newS;bigLambda(bestMap).*L];
% fhat = pinv(A)*constraints;
% figure();subplot(2,1,1);imagesc(gabor);
% subplot(2,1,2);imagesc(reshape(fhat,[N,N]));
% title(sprintf('Lambda: %3.0e',bigLambda(bestMap)));