function [] = Noise_RevCorr(AnimalName,NoiseType,DistToScreen,flipInterval,WaitTime)
%Noise_ReverseCorrelation.m
%   Display a series of white noise stimuli to infer the receptive fields
%    of neurons using reverse correlation.
%    See Smyth et al. 2003 Receptive Field Organization ...
%
% Briefly, if the receptive field for a given neuron is described by a
%  p-by-1 vector (p being the number of pixels), f, and the stimuli described
%  by a s-by-p (s being the number of stimuli presented) matrix, S, and the 
%  response of that neuron described by an s-by-1 vector, then we can write
%  the responses in terms of the stimuli and the receptive field as 
%  r = S * f .  We want to infer the spatial receptive field structure
%  f, so we need S^-1 * r = S^-1 * S * f => f = S^-1 * r ... So, the
%  following code will present a series of white noise stimuli with each
%  stimulus lasting flipInterval msec, followed by flipInterval msec of flat 
%  grey.  The stimuli will be output as the matrix S, along with the 
%  timeStamps for when they were generated.
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Optional Inputs
%       NoiseType - 'white' or 'pink' ... defaults to white
%       DistToScreen - distance from the mouse to the screen, in
%        centimeters
%       flipInterval - time (milliseconds) to display the noise, then the
%          grey screen ... the screen will flip from grey to noise to grey
%          and each will display for flipInterval milliseconds
%       WaitTime - time (milliseconds) during which microscope will record
%          the visually-evoked response
%
%OUTPUT: file named 'NoiseStimDate_AnimalName.mat' , e.g. 
%          NoiseStim20160718_12345.mat
%       S - matrix sized numStimuli-by-numPixels that represent each of
%          the stimuli presented, try 
%          image = reshape(S(1,:),[width,height]); to view one of the white
%          noise stimuli
%        effectivePixels - effective width (and height) of the display window 
%          in pixels, first chosen as the minimum of the width and height 
%          of the current display ... the display is forced to be square
%          and the width in effectivePixels is 1/4 the width in true screen
%          pixels
%        DistToScreen - as above
%
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/18
% By: Byron Price

switch nargin
    case 1
        NoiseType = 'white';
        DistToScreen = 25;
        flipInterval = 200;
        WaitTime = 1000;
    case 2
        DistToScreen = 25;
        flipInterval = 200;
        WaitTime = 1000;
    case 3
        flipInterval = 200;
        WaitTime = 1000;
    case 4
        WaitTime = 1000;
end
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;
usb = usb1208FSPlusClass;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

flipInterval = flipInterval/1000;
WaitTime = WaitTime/1000;
numStimuli = 1500;
TimeEstimate = numStimuli*(flipInterval+WaitTime)/60;
display(sprintf('Estimated time is %.2f minutes.',TimeEstimate))
WaitSecs(5);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display,
[win,~] = Screen('OpenWindow', screenid,0);

% Query window size in pixels
[width, height] = Screen('WindowSize', win);
minPix = min(width,height);

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to 4x4 squares of on-screen pixels
screenPix_to_effPix = 4;
minPix = minPix-mod(minPix,screenPix_to_effPix);
numPixels = minPix*minPix;

effectivePixels = numPixels/(screenPix_to_effPix*screenPix_to_effPix);
% uniformly-distributed noise (Is Gaussian-distributed noise white?)
if strcmp(NoiseType,'white') == 1
    S = randi([0,255],[numStimuli,effectivePixels],'uint8');
elseif strcmp(NoiseType,'pink') == 1
    S = randi([0,255],[numStimuli,effectivePixels],'uint8');
    N = sqrt(effectivePixels);
    for ii=1:numStimuli
        stim = reshape(S(ii,:),[N,N]);
        y = fft2(stim);
        xfreq = (1./(1:N))./40;
        xfreq = bsxfun(@times,xfreq,ones(length(xfreq),1));
        yfreq = xfreq';
        mask = 1./sqrt(xfreq.^2+yfreq.^2);
        mask = rot90(mask,2);
        y = y.*mask;
        stim = real(ifft2(y));
        stim = reshape(stim,[N*N,1]);
        stim = stim-min(stim);
        stim = round(stim.*(255/max(stim)));
        S(ii,:) = reshape(uint8(stim),[1,N*N]);
    end
else 
    display('NoiseType must be ''white'' or ''pink'' as a string.')
    return;
end

Grey = 128*ones(minPix,minPix);
timeStamps = zeros(numStimuli*2,1);

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

% usb.startRecording;
WaitSecs(2);

vbl = Screen('Flip', win);
for tt=1:numStimuli
    % Convert it to a texture 'tex':
    Img = reshape(S(tt,:),[minPix/screenPix_to_effPix,minPix/screenPix_to_effPix]);
    Img = kron(double(Img),ones(screenPix_to_effPix));
    tex = Screen('MakeTexture',win,Img);
    clear Img;
    Screen('DrawTexture',win, tex);
    vbl = Screen('Flip',win, vbl + ifi/2);%
%     usb.strobe;
    WaitSecs(flipInterval);vbl = vbl+flipInterval;

    Img = Grey;
    tex = Screen('MakeTexture',win,Img);
    clear Img;
    Screen('DrawTexture',win, tex);
    vbl = Screen('Flip',win, vbl + ifi/2);
%     usb.strobe;
    WaitSecs(WaitTime);vbl = vbl+WaitTime;
    Screen('Close', tex);
end
% usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

cd('~/Documents/MATLAB/Byron/RetinoExp')
Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');
filename = strcat('NoiseStim',Date,'_',num2str(AnimalName),'.mat');
save(filename,'S','numStimuli','effectivePixels','DistToScreen');
end
