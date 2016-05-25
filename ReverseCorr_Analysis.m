function [F] = ReverseCorr_Analysis(S,R,timeStamps,Fs)
%ReverseCorr_Analysis.m
%   %   Analysis of calcium-imaging data in response to a series of white
%   noise stimuli (see WhiteNoise_ReverseCorrelation.m for Psychtoolbox
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
%INPUT: S - number of stimuli-by-number of pixels matrix of white noise
%        stimuli used by the function WhiteNoise_ReverseCorrelation.m
%       R - matrix of size number of frames-by-number of neurons that
%        contains neuronal data which has already been processed into a
%        point process and which will be converted to a new matrix of responses, 
%        number of stimuli-by-number of neurons matrix of neuronal
%        numbers obtained from calcium imaging (count threshold
%        crossings as spikes ... we will want to know the number of spikes
%        in response to the onset and offset of each stimuli.
%       timeStamps - timeStamps output from WhiteNoise_ReverseCorrelation.m
%       Fs - sampling frequency of the calcium imaging recording in Hz 
%OUTPUT: F - number of pixels-by-number of neurons matrix representing the
%         receptive field for each neuron recorded from in R
%            if N = sqrt(size(F,1));
%            then the receptive field for a given neuron can be 
%            visualized by:
%            figure();imagesc(reshape(F(:,ii),[N,N]);      
%    
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/03/05
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

% assume we want the period from 30msec to 60msec after stimulus
%  onset and that the calcium imaging recording started just in time 
%  with the first stimulus (which was a grey screen)

period = [30,60]; % period after stimulus to include for spike-triggered
       % average stimulus (in milliseconds)
startRecord = round(period(1)/1000*Fs);
endRecord = round(period(2)/1000*Fs);
numNeurons = size(R,2);
numStimuli = size(S,1)*2;

RR = zeros(numStimuli/2,numNeurons);
for ii=1:numNeurons
    Trace = square(R(:,ii));
    currentIndex = 1;
    for tt=2:numStimuli
        interval = timeStamps(tt)-timeStamps(tt-1);
        numFrames = interval*Fs;
        currentIndex = currentIndex+numFrames;
        response = sum(Trace(currentIndex+startRecord:currentIndex+endRecord));
        if mod(tt,2) == 0
            RR(tt/2,ii) = response;
        elseif mod(tt,2) == 1
            RR((tt-1)/2,ii) = -response;
        end
    end
end


effectivePixels = size(S,2);
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

lambda = 10;
% VERTICALLY CONCATENATE S and L
A = [S;lambda.*L];
F = zeros(effectivePixels,numNeurons);

% REGULARIZED LEAST-SQUARES SOLUTION
for ii=1:numNeurons
    r = RR(:,ii);
    constraints = [r;zeros(effectivePixels,1)];
    
    %Singular Value Decomposition to solve
    % linear system of equations for the receptive
    % field of each neuron, f
    % given a matrix A, solve
    % Af = constraints for f
    % first require U,S,V such that A = U*S*V'
    [U,S,V] = svd(A);
    % pseudo-inverse is then
    Sinv = inv(S);
    Ainv = V'*Sinv*U;
    f = Ainv*constraints;
    % alternatively 
    % f = A\constraints;
    F(:,ii) = f;
end
end

