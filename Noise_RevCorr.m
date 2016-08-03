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
%       NoiseType - 'white' or 'pink' or 'brown' ... defaults to brown
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
% Updated: 2016/08/03
% By: Byron Price

cd('~/CloudStation/ByronExp/RetinoExp')
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

flipInterval = flipInterval/1000;
WaitTime = WaitTime/1000;
numStimuli = 1500;
TimeEstimate = numStimuli*(flipInterval+WaitTime)/60;
display(sprintf('\nEstimated time is %3.2f minutes.',TimeEstimate))
WaitSecs(10);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display,
[win,~] = Screen('OpenWindow', screenid,128);

% Query window size in pixels
[w_pixels,h_pixels] = Screen('WindowSize', win);
minPix = min(w_pixels,h_pixels);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to larger squares of on-screen pixels
screenPix_to_effPix = 20;
minPix = minPix-mod(minPix,screenPix_to_effPix);
numPixels = minPix*minPix;

effectivePixels = numPixels/(screenPix_to_effPix*screenPix_to_effPix);

% GENERATION OF NOISE
if strcmp(NoiseType,'white') == 1
    beta = 0;
elseif strcmp(NoiseType,'pink') == 1
    beta = -1;
elseif strcmp(NoiseType,'brown') == 1
    beta = -2;
else 
    display('NoiseType must be ''white'' or ''pink'' as a string.')
    return;
end

S = zeros(numStimuli,effectivePixels);
N = sqrt(effectivePixels);
% perform unit conversions
degPerPix = atan((minPix*conv_factor)/(DistToScreen*10))./N;
% below pink noise from Jon Yearsley, 1/f noise generate spatial data
DIM = [N,N];
for ii=1:numStimuli
    u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
    u = repmat(u,1,DIM(2));
    v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]/DIM(2);
    v = repmat(v,DIM(1),1);
    S_f = (u.^2 + v.^2).^(beta);
    S_f(S_f==inf) = 0;
    phi = rand(DIM);
    x = ifft2(S_f.^0.5 .* (cos(2*pi*phi)+1i*sin(2*pi*phi)));
    x = real(x);
    x = x+abs(min(min(x)));
    x = round((x./max(max(x))).*255);
    S(ii,:) = reshape(x,[1,effectivePixels]);
end


Grey = 128*ones(minPix,minPix);

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

usb.startRecording;
WaitSecs(5);
for tt=1:numStimuli
    vbl = Screen('Flip', win);
    % Convert it to a texture 'tex':
    Img = reshape(S(tt,:),[minPix/screenPix_to_effPix,minPix/screenPix_to_effPix]);
    Img = kron(double(Img),ones(screenPix_to_effPix));
    tex = Screen('MakeTexture',win,Img);
    clear Img;
    Screen('DrawTexture',win, tex);
    vbl = Screen('Flip',win);usb.strobe;
    vbl = Screen('Flip',win,vbl-ifi/2+flipInterval);usb.strobe;
    vbl = Screen('Flip',win,vbl-ifi/2+WaitTime);
end
usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('NoiseStim%d_%d.mat',Date,AnimalName);
save(filename,'S','numStimuli','flipInterval','effectivePixels','DistToScreen','screenPix_to_effPix','minPix','NoiseType');
end
