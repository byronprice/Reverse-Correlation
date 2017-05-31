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
% Updated: 2017/05/18
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

TimeEstimate = numStimuli*(flipInterval+WaitTime+0.2)/60;
fprintf('\nEstimated time is %3.2f minutes.',TimeEstimate);
WaitSecs(5);

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
maxPix = max(w_pixels,h_pixels)-200;

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to larger squares of on-screen pixels
screenPix_to_effPix = 10;
maxPix = maxPix-mod(maxPix,screenPix_to_effPix);
minPix = minPix-mod(minPix,screenPix_to_effPix);

effectivePixels = [maxPix/screenPix_to_effPix,minPix/screenPix_to_effPix];

gridSpacing = screenPix_to_effPix*conv_factor/10;
spatialSampleFreq = 1/gridSpacing; % in units of 1/cm

% GENERATION OF NOISE
if strcmp(NoiseType,'white') == 1
    beta = 0;
    spaceExp = 2;
elseif strcmp(NoiseType,'pink') == 1
    beta = -2;
    spaceExp = 2;
elseif strcmp(NoiseType,'brown') == 1
    beta = -4;
    spaceExp = 2;
else 
    display('NoiseType must be ''white'', ''pink'' or ''brown'' as a string.')
    return;
end

Grey = 127;
DIM = [effectivePixels(1),effectivePixels(2)];
S = zeros(numStimuli,DIM(1)*DIM(2),'uint8');
% below pink noise from Jon Yearsley, 1/f noise generate spatial data

for ii=1:numStimuli
    X = spatialPattern([DIM(2),DIM(1)],beta);
    X = X-min(min(X));
    X = (X./max(max(X))).*255;
    Y = X(:);
    meanVal = mean(Y);difference = meanVal-Grey;
    %figure();imagesc(reshape(Y-difference,DIM));
    S(ii,:) = Y-difference;
end
S = uint8(S);

wLow = round((w_pixels-maxPix)/2);
wHigh = round(w_pixels-wLow);
hLow = round((h_pixels-minPix)/2);
hHigh = round(h_pixels-hLow);
destRect = [wLow hLow wHigh hHigh];

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

WaitTimes = WaitTime+exprnd(0.2,[numStimuli,1]);

usb.startRecording;usb.strobeEventWord(0);
WaitSecs(30);
tt = 1;
vbl = Screen('Flip',win);
while tt <= numStimuli
    % Convert it to a texture 'tex':
    Img = reshape(S(tt,:),[DIM(2),DIM(1)]);
    tex = Screen('MakeTexture',win,Img);
    Screen('DrawTexture',win, tex,[],destRect,[],0); % 0 is nearest neighbor
                                        % 1 is bilinear filter
    vbl = Screen('Flip',win);usb.strobeEventWord(1);
    vbl = Screen('Flip',win,vbl-ifi/2+flipInterval);usb.strobeEventWord(2);
    vbl = Screen('Flip',win,vbl-ifi/2+WaitTimes(tt));
    Screen('Close',tex);
    tt = tt+1;
end
Screen('Flip',win);usb.strobeEventWord(0);
WaitSecs(30);
usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('NoiseStim%s%d_%d.mat',NoiseType,Date,AnimalName);
save(filename,'S','numStimuli','flipInterval','effectivePixels',...
    'DistToScreen','screenPix_to_effPix','minPix','NoiseType','degPerPix',...
    'conv_factor','WaitTimes','beta','spaceExp','DIM','spatialSampleFreq');
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

