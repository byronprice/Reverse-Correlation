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

EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');

if exist(EphysFileName,'file') ~= 2
    readall(strcat(EphysFileName,'.plx'));pause(1);
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
numStimuli = numStimuli*2;
newS = zeros(numStimuli,size(S,2),'single');
for ii=1:numStimuli
    if mod(ii,2) == 1
        newS(ii,:) = S(floor(ii/2)+1,:);
    elseif mod(ii,2) == 0
        newS(ii,:) = 255-S(ii/2,:); %S(ii/2,:)
%         temp = reshape(S(ii/2,:),[N,N]);
%         temp2 = reshape(255-S(ii/2,:),[N,N]);
%         figure();imagesc(temp);
%         figure();imagesc(temp2);
%         display('blah');
    end
end
clear S;

strobeData = tsevs{1,strobeStart};

% GLM to check for visual responsiveness
significantVisualResponsiveness = zeros(size(allts));
stimStartTimes = zeros(size(allts));
for ii=1:numChans
    numNeurons = Neurons(ii,1);
    for jj=1:numNeurons
        totalTime = strobeData(end)+2;
        totalMillisecs = round(totalTime*1000);
        histDependentCoeffs = 150;
        timeAfterStimOnsetCoeffs = 150;
        
        % create Gaussian basis functions with 2ms standard deviation
        stdev = 5;
        numBasis = 25;
        Shist = zeros(histDependentCoeffs,numBasis);
        SstimOnset = zeros(timeAfterStimOnsetCoeffs,numBasis);
        
        timeHist = 0:(histDependentCoeffs-1);
        timeStimOnset = 0:(timeAfterStimOnsetCoeffs-1);
        for kk=1:numBasis
            Shist(:,kk) = exp(-(timeHist'-(kk-1)*histDependentCoeffs/numBasis).^2./(2*stdev*stdev));
            SstimOnset(:,kk) = exp(-(timeStimOnset'-(kk-1)*timeAfterStimOnsetCoeffs/numBasis).^2./(2*stdev*stdev));
        end

        Y = zeros(totalMillisecs,1);
        visualResponseDesign = zeros(totalMillisecs,timeAfterStimOnsetCoeffs);
        histDependentDesign = zeros(totalMillisecs,histDependentCoeffs);
        
        % convert spike times to point process
        temp = round(allts{ii,jj+1}.*1000);
        for kk=1:length(temp)
            Y(temp(kk)) = 1;
        end
        
        % convert timeStamps to point process
        stimTimes = round(strobeData.*1000);
        pointProcessStimOnsets = zeros(totalMillisecs,1);
        for kk=1:numStimuli
            pointProcessStimOnsets(stimTimes(kk)) = 1;
        end
        
        % get history dependence (past 20ms)
        for kk=1:histDependentCoeffs
            temp = Y;shift = zeros(kk,1);
            history = [shift;temp];
            histDependentDesign(:,kk) = history(1:(end-kk));
        end
        
        % collect time period from 1 to 150 ms
        for kk=1:timeAfterStimOnsetCoeffs
            temp = pointProcessStimOnsets;shift = zeros(kk,1);
            history = [shift;temp];
            visualResponseDesign(:,kk) = history(1:(end-kk));
        end
        
        Design = histDependentDesign*Shist;
        [b1,dev1,stats1] = glmfit(Design,Y,'poisson');
        AIC1 = dev1+2*length(b1);
        
        Design = [histDependentDesign*Shist,visualResponseDesign*SstimOnset];
        [b2,dev2,stats2] = glmfit(Design,Y,'poisson');
        AIC2 = dev2+2*length(b2);
        
        if AIC2 < AIC1
            significantVisualResponsiveness(ii,jj+1) = 1;
            fprintf('Chan %d, Neuron %d Visually Responsive\n',ii,jj);
        else
            fprintf('Chan %d, Neuron %d NOT Visually Responsive\n',ii,jj); 
        end
        
        fullS = [[Shist;zeros(size(SstimOnset))],[zeros(size(Shist));SstimOnset]];
        
        [yhat,dylo,dyhi] = glmval(b2,fullS,'log',stats2);
        
        [~,index] = max(abs(yhat(histDependentCoeffs+1:end)-...
            mean(yhat(histDependentCoeffs:histDependentCoeffs+45))));
        stimStartTimes(ii,jj+1) = index-30;

        figure();boundedline(1:size(fullS,1),yhat,[dylo,dyhi]);
        set(gca,'XTick',[linspace(0,histDependentCoeffs,10),...
            linspace(histDependentCoeffs+1,histDependentCoeffs+timeAfterStimOnsetCoeffs,10)]);
        set(gca,'XTickLabels',[round(linspace(0,histDependentCoeffs,10)),...
        round(linspace(0,timeAfterStimOnsetCoeffs,10))]);
        xlabel('Lag (ms), then Time After Stim Onset (ms)');
        ylabel('Intensity');
        title(sprintf('Chan %d, Neuron %d- %d',ii,jj,AIC2<AIC1));
    end
end

% CREATE LAPLACIAN MATRIX
L = zeros(effectivePixels,effectivePixels,'single');
for ii=1:effectivePixels
    if ii == 1  % top left corner
        L(ii,1) = 2;
        L(ii,2) = -1;
        L(ii,N+1) = -1;
    elseif ii == N % bottom left corner
        L(ii,N) = 2;
        L(ii,N-1) = -1;
        L(ii,2*N) = -1;
    elseif ii == ((N-1)*N+1) % top right corner
        L(ii,(N-1)*N+1) = 2;
        L(ii,(N-2)*N+1) = -1;
        L(ii,(N-1)*N+2) = -1;
    elseif ii == N*N % bottom right corner
        L(ii,N*N) = 2;
        L(ii,N*N-1) = -1;
        L(ii,(N-1)*N) = -1;
    elseif ismember(ii,2:1:N-1) % left edge
        L(ii,ii) = 3;
        L(ii,ii-1) = -1;
        L(ii,ii+1) = -1;
        L(ii,ii+N) = -1;
    elseif ismember(ii,2*N:N:(N-1)*N) % bottom edge
        L(ii,ii) = 3;
        L(ii,ii-1) = -1;
        L(ii,ii-N) = -1;
        L(ii,ii+N) = -1;
    elseif ismember(ii,(N-1)*N+2:1:N*N-1) % right edge
        L(ii,ii) = 3;
        L(ii,ii-1) = -1;
        L(ii,ii+1) = -1;
        L(ii,ii-N) = -1;
    elseif ismember(ii,(N+1):N:(N-2)*N+1) % top edge
        L(ii,ii) = 3;
        L(ii,ii+1) = -1;
        L(ii,ii-N) = -1;
        L(ii,ii+N) = -1;
    else        % interior points
        L(ii,ii) = 4;
        L(ii,ii-1) = -1;
        L(ii,ii+1) = -1;
        L(ii,ii-N) = -1;
        L(ii,ii+N) = -1;
    end
end

% could simplify by ignoring edge conditions and doing ...
%   diagFours = 4.*diag(effectivePixels);
%   diagUpOnes = -1.*diag(effectivePixels,1);
%   diagDownOnes = -1.*diag(effectivePixels,-1);
%   L = diagFours+diagUpOnes+diagDownOnes;


stimLen = 0.06;

totalTime = strobeData(end)+2;
% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = zeros(numChans,nunits1,numStimuli);
baseRate = zeros(numChans,nunits1);
for ii=1:numChans
    numNeurons = Neurons(ii,1);
    for jj=1:numNeurons
        stimStart = stimStartTimes(ii,jj+1);
        baseRate(ii,jj+1) = length(allts{ii,jj+1})./totalTime;
        %             display(baseRate(ii,jj));
        %         figure();plot(0,0);axis([0 totalTime -10 10]);hold on;
        for kk=1:numStimuli
            stimOnset = strobeData(kk);
            high = find(allts{ii,jj+1} > (stimOnset+stimStart));
            low = find(allts{ii,jj+1} < (stimOnset+stimStart+stimLen));
            temp = intersect(low,high);
            
            Response(ii,jj+1,kk) = (length(temp)./stimLen)./baseRate(ii,jj+1);
            
            clear temp;
        end
    end
end


% REGULARIZED PSEUDO-INVERSE SOLUTION
degreesVisualSpace = degPerPix*screenPix_to_effPix*N;
xaxis = linspace(-degreesVisualSpace/2,degreesVisualSpace/2,N);
yaxis = linspace(-degreesVisualSpace/4,3*degreesVisualSpace/4,N);
maxNeurons = max(Neurons);

bigLambda = [5e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6,5e6];
F = zeros(length(bigLambda),numChans,nunits1,effectivePixels);
RMS = zeros(length(bigLambda),numChans,nunits1);
for lambda = 1:length(bigLambda)
    % VERTICALLY CONCATENATE S and L
    % S is size numStimuli X effectivePixels
    A = [newS;bigLambda(lambda).*L];
    for ii=1:numChans
        numNeurons = Neurons(ii,1);
        for jj=1:numNeurons
            
            r = squeeze(Response(ii,jj+1,:));
            constraints = [r;zeros(effectivePixels,1)];
            % %         r = newS*f ... constraints = A*f
            %         [fhat,~,~] = glmfit(A,constraints,'normal','constant','off');
            %         [rhat,lBound,uBound] = glmval(fhat,S,'identity',stats,'confidence',1-alpha,'constant','off');
            
            %         fhat = pinv(A)*constraints;
            %         fhat = pinv(newS)*r;
            % alternatively
            %         fhat = newS\r;
            
            fhat = A\constraints;
            %         fhat = newS(1:2:end,:)'*r(1:2:end)/sum(r(1:2:end));
            %         fhat = newS'*r./sum(r);
            
            tempfhat = reshape(fhat,[N,N]);
            tempFFT = fft2(tempfhat);
            tempFFT = numStimuli.*tempFFT./S_f;
            tempIFFT = ifft2(tempFFT);%tempIFFT = sqrt(tempIFFT.*conj(tempIFFT));
            F(lambda,ii,jj+1,:) = tempIFFT(:);
            RMS(lambda,ii,jj+1) = (effectivePixels)^(-0.5)*norm(r-newS*tempIFFT(:),2);
%             subplot(numChans,maxNeurons,plotCount);
%             imagesc(xaxis,yaxis,tempIFFT);set(gca,'YDir','normal');
%             title(sprintf('Lambda %3.1e',bigLambda(lambda)));
%             xlabel('Azimuth (degrees of visual arc)');
%             ylabel('Altitude(degrees of visual arc)');
%             plotCount = plotCount+1;
        end
    end
end

% for ii=1:numChans
%     numNeurons = Neurons(ii,1);
%     for jj=1:numNeurons
%         
%     end
% end

FileName = strcat('NoiseResults',NoiseType,num2str(Date),'_',num2str(AnimalName),'.mat');
save(FileName,'F','newS','Response','allts','bigLambda','numChans','nunits1',...
    'Neurons','beta','spaceExp','N','significantVisualResponsiveness','RMS',...
    'stimStartTimes');

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






