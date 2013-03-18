function c_noise = add_observation_noise(c, noise)
%add_observation_noise adds artificial Gaussian distributed noise to
% synthetic receptor observations.
%
% INPUTS:
%   c         The receptor observations matrix.
%
%   noise     Noise mode parameter which controls the type of noise 
%               which is added the receptor observations .
%
%               noise > 0   the slope of the line which relates the 
%                             observation concentration to the standard
%                             deviation of the Gaussian distributed 
%                             noise (0.1 = 10% noise).
%
%               noise = 0   no noise is added to the receptor 
%                             observations.
%
%               noise = -1  the standard deviation of the added noise
%                             follows the estimated TDLAS sensitivity
%                             curve.
%
% OUTPUT:
%   c_noise   The modified receptor observation matrix with Gaussian 
%               distributed observational noise and has the same
%               dimension as 'c'.

assert(isnumeric(c));
assert(isnumeric(noise));
assert(isscalar(noise));
assert( noise >= 0 || noise == -1 );

% conversion factor obtained from:
%http://www.lenntech.com/calculators/ppm/converter-parts-per-million.htm
kg_m3_per_ppm = 0.706E-6; 

dqII = [[      -9E99, 27.25700089]; [5.000000000, 27.25700089]; ...
        [55.00000000, 27.58277605]; [105.0000000, 27.93364826]; ...
        [155.0000000, 28.30884459]; [205.0000000, 28.70757301]; ...
        [255.0000000, 29.12903028]; [305.0000000, 29.57240881]; ...
        [355.0000000, 30.03690350]; [405.0000000, 30.52171671]; ...
        [455.0000000, 31.02606363]; [505.0000000, 31.54917604]; ...
        [555.0000000, 32.09030604]; [605.0000000, 32.64872808]; ...
        [655.0000000, 33.22374174]; [705.0000000, 33.81467295]; ...
        [755.0000000, 34.42087499]; [805.0000000, 35.04172937]; ...
        [855.0000000, 35.67664589]; [905.0000000, 36.32506299]; ...
        [955.0000000, 36.98644688]; [1005.000000, 37.66029193]; ...
        [1055.000000, 38.34611920]; [1105.000000, 39.04347620]; ...
        [1155.000000, 39.75193600]; [1205.000000, 40.47109625]; ...
        [1255.000000, 41.20057794]; [1305.000000, 41.94002512]; ...
        [1355.000000, 42.68910327]; [1405.000000, 43.44749841]; ...
        [1455.000000, 44.21491645]; [1505.000000, 44.99108199]; ...
        [1555.000000, 45.77573711]; [1605.000000, 46.56864102]; ...
        [1655.000000, 47.36956883]; [1705.000000, 48.17831088]; ...
        [1755.000000, 48.99467163]; [1805.000000, 49.81846955]; ...
        [1855.000000, 50.64953529]; [1905.000000, 51.48771219]; ...
        [1955.000000, 52.33285506]; [2005.000000, 53.18482946]; ...
        [9E99, 53.18482946] ] * kg_m3_per_ppm;
 

if noise > 0  % use a percentage of the obsevation
  min_noise = dqII(1,2) % 27.25700089*kg_m3_per_ppm;
  dc = noise*c + min_noise;
  dc(dc < min_noise) = min_noise;
  c_error = dc.*randn(size(c));
elseif noise == -1  % use dqII curve
  c_error = interp1(dqII(:,1),dqII(:,2),c).*randn(size(c));
elseif  noise == 0  % no noise
  c_error = zeros(size(c));
else
  error('ERROR: Invalid observation noise parameter');
end

c_noise = c + c_error;