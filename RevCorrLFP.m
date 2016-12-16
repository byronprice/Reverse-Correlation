function [] = RevCorrLFP(AnimalName,Date,NoiseType)
%RevCorrLFP.m
%   %   Analysis of LFP recording data in response to a series of white
%   noise stimuli (see Noise_RevCorr.m for Psychtoolbox
%   stimulus)
%    See Smyth et al. 2003 Receptive Field Organization ... We'll solve a
%    matrix equation using a regularized pseudo-inverse technique.
%
% Briefly, if the receptive field for a given neuron is described by a
%  p-by-1 vector (p being the number of pixels), f, and the stimuli described
%  by a s-by-p (s being the number of stimuli presented) matrix, S, and the 
%  response of that neuron described by an s-by-1 vector, then we can write
%  the responses in terms of the stimuli and the receptive field as 
%  r = S * f .  We want to infer the spatial receptive field structure
%  f, so we need S^-1 * r = S^-1 * S * f => f = S^-1 * r ... So, the
%  following code will get the spatial receptive fields f for a series of 
%  neurons in R whose responses were measured with calcium imaging during
%  the display of the white noise stimuli in S.
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
% Created: 2016/07/25, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/25
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
EphysFileName = strcat('NoiseData',NoiseType,num2str(Date),'_',num2str(AnimalName));

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

StimulusFileName = strcat('NoiseStim',num2str(Date),'_',num2str(AnimalName),'.mat');
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)
load(StimulusFileName)

sampleFreq = adfreq;
% newFs = 5;
% downSampleRate = sampleFreq/newFs;


% tsevs are the strobed times of stimulus onset, then offset
%  Onset at tsevs{1,33}(2), offset at tsevs{1,33}(3), onset at
%  tsevs{1,33}(4), offset at 5, etc.
% allad contains the continuous data from each channel, which appear to be
%  recorded at 1000 Hz rather than 40,000

%totalAD = size(allad,2);
%totalSEVS = size(tsevs,2);

Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);
strobeStart = 33;

% lowpass filter the data
dataLength = length(allad{1,Chans(1)});

ChanData = zeros(dataLength,numChans);
preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    temp = smooth(voltage,0.013*sampleFreq);
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,temp);
end


timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeData = tsevs{1,strobeStart};


% ASSUME THAT A RESPONSE TO A STIMULUS OFFSET IS THE SAME AS A RESPONSE TO
%  THE NEGATIVE OF THAT IMAGE, image created with values from 0 to 255
numStimuli = numStimuli*2;
newS = zeros(numStimuli,size(S,2),'single');
for ii=1:numStimuli
    if mod(ii,2) == 1
        newS(ii,:) = S(floor(ii/2)+1,:);
    elseif mod(ii,2) == 0
        newS(ii,:) = 255-S(ii/2,:);
    end
end
clear S;

strobeData = tsevs{1,strobeStart};
stimLen = round(0.05*sampleFreq);
stimStart = round(0.05*sampleFreq);

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = zeros(numChans,numStimuli);
for ii=1:numChans
    for jj=1:numStimuli
        stimOnset = strobeData(jj);
        [~,index] = min(abs(timeStamps-stimOnset));
        temp = ChanData(index+stimStart:index+stimStart+stimLen-1,ii);
        Response(ii,jj) = abs(min(temp));
    end
end

% CREATE LAPLACIAN MATRIX
N = sqrt(effectivePixels);
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

lambda = 1e4;
% VERTICALLY CONCATENATE S and L
% S is size numStimuli X effectivePixels
A = [newS;lambda.*L];
F = zeros(numChans,effectivePixels);

% REGULARIZED PSEUDO-INVERSE SOLUTION
figure();
for ii=1:numChans
    r = squeeze(Response(ii,:)); 
    constraints = [r';zeros(effectivePixels,1)];
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
    F(ii,:) = fhat;
    subplot(numChans,1,ii);imagesc(reshape(fhat,[N,N]));
    title(sprintf('Laplacian Map for Channel %d , Animal %d',ii,AnimalName));
    %         subplot(plotRows,2,jj);histogram(r);
end

% REVERSE CORRELATION SOLUTION
figure();
for ii=1:numChans
    r = squeeze(Response(ii,:));
    % %         r = newS*f ... constraints = A*f
    fhat = newS\r';
    F(ii,:) = fhat;
    subplot(numChans,1,ii);imagesc(reshape(fhat,[N,N]));
    title(sprintf('RevCorr Map for Channel %d , Animal %d',ii,AnimalName));
end
% Img = reshape(S(tt,:),[minPix/screenPix_to_effPix,minPix/screenPix_to_effPix]);
% Img = kron(double(Img),ones(screenPix_to_effPix));

end






