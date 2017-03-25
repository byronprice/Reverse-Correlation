function [] = Noise_RevCorr(AnimalName,NoiseType)
%Noise_RevCorr.m
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
%       NoiseType - 'white' or 'pink' or 'brown' ... defaults to pink
%
%       see file NoiseVars.mat for more changeable presets
%
%OUTPUT: file named 'NoiseStimDate_AnimalName.mat' , e.g. 
%          NoiseStim20160718_12345.mat
%
%        S - matrix sized numStimuli-by-numPixels that represent each of
%          the stimuli presented, try 
%          image = reshape(S(1,:),[width,height]); to view one of the white
%          noise stimuli
%        effectivePixels - effective width (and height) of the display window 
%          in pixels, first chosen as the minimum of the width and height 
%          of the current display ... the display is forced to be square
%          and the width in effectivePixels is 1/4 the width in true screen
%          pixels
%        DistToScreen - 25 cm for now
%
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/02/03
% By: Byron Price

cd('~/CloudStation/ByronExp/NoiseRetino')
load('NoiseVars.mat');

switch nargin
    case 1
        NoiseType = 'pink';
end
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;
usb = usb1208FSPlusClass;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

TimeEstimate = numStimuli*(flipInterval+0.1+WaitTime+0.2)/60;
display(sprintf('\nEstimated time is %3.2f minutes.',TimeEstimate))
WaitSecs(10);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
background = 127;
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

% Query window size in pixels
[w_pixels,h_pixels] = Screen('WindowSize', win);
minPix = min(w_pixels,h_pixels);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to larger squares of on-screen pixels
screenPix_to_effPix = 25;
minPix = minPix-mod(minPix,screenPix_to_effPix);
numPixels = minPix*minPix;

effectivePixels = numPixels/(screenPix_to_effPix*screenPix_to_effPix);

% GENERATION OF NOISE
if strcmp(NoiseType,'white') == 1
    beta = 0;
elseif strcmp(NoiseType,'pink') == 1
    beta = -2;
elseif strcmp(NoiseType,'brown') == 1
    beta = -4;
else 
    display('NoiseType must be ''white'', ''pink'' or ''brown'' as a string.')
    return;
end

Grey = 127;
S = zeros(numStimuli,effectivePixels,'uint8');
N = sqrt(effectivePixels);
% perform unit conversions
degPerPix = atand((1*conv_factor)/(DistToScreen*10));
% below pink noise from Jon Yearsley, 1/f noise generate spatial data
DIM = [N,N];
for ii=1:numStimuli
    u = 0:(N-1);
    [U,V] = meshgrid(u,u);
    S_f = (U.^2+V.^2).^(beta/2);
    S_f(S_f==inf) = 1;
    noise = wgn(DIM(1),DIM(2),0);
    Y = fft2(noise);
    Y = Y.*S_f;
    X = ifft2(Y);
    X = sqrt(X.*conj(X)); % real(X)
    X = X-min(min(X));
    X = (X./max(max(X))).*255;
    y = reshape(X,[1,effectivePixels]);
    meanVal = mean(y);difference = meanVal-Grey;
    S(ii,:) = y-difference;
end
S = uint8(S);

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);
flipIntervals = flipInterval+exprnd(0.1,[numStimuli,1]);
WaitTimes = WaitTime+exprnd(0.2,[numStimuli,1]);

usb.startRecording;
WaitSecs(5);
vbl = Screen('Flip',win);
for tt=1:numStimuli
    % Convert it to a texture 'tex':
    Img = reshape(S(tt,:),[minPix/screenPix_to_effPix,minPix/screenPix_to_effPix]);
    Img = kron(double(Img),ones(screenPix_to_effPix));
    tex = Screen('MakeTexture',win,Img);
    Screen('DrawTexture',win, tex);
    vbl = Screen('Flip',win);usb.strobe;
    vbl = Screen('Flip',win,vbl-ifi/2+flipIntervals(tt));usb.strobe;
    vbl = Screen('Flip',win,vbl-ifi/2+WaitTimes(tt));
    Screen('Close',tex);
end
pause(1);
usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('NoiseStim%s%d_%d.mat',NoiseType,Date,AnimalName);
save(filename,'S','numStimuli','flipIntervals','effectivePixels',...
    'DistToScreen','screenPix_to_effPix','minPix','NoiseType','degPerPix','WaitTimes');
end

function gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
% Generates a 256x3 gamma lookup table suitable for use with the
% psychtoolbox Screen('LoadNormalizedGammaTable',win,gammaTable) command
% 
% gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
%
%   gamma defines the level of gamma correction (1.8 or 2.2 common)
%   blackSetPoint should be the highest value that results in a non-unique
%   luminance value on the monitor being used (sometimes values 0,1,2, all
%   produce the same black pixel value; set to zero if this is not a
%   concern)
%   whiteSetPoint should be the lowest value that returns a non-unique
%   luminance value (deal with any saturation at the high end)
% 
%   Both black and white set points should be defined on a 0:255 scale

gamma = max([gamma 1e-4]); % handle zero gamma case
gammaVals = linspace(blackSetPoint/255,whiteSetPoint/255,256).^(1./gamma);
gammaTable = repmat(gammaVals(:),1,3);
end

