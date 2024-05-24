
function out = applyCT_1D( f, domTraj, rangeTraj, kC, C )
  % out = applyCT_1D( f, domTraj, rangeTraj, kC, C )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   f - An nTraj array representing the values of the function evaluated at each point in domTraj
  %   domTraj - An nTraj element array specifying the coordinates of the domain points
  %   rangeTraj - An nNew element array specifying the coordinates of the new points
  %   kC - array of convolution kernel domain values
  %   C  - array of convolution kernel values
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  % Adjoint is just applyC with the trajectories reversed
  out = applyC_1D( f, rangeTraj, domTraj, kC, C );
end

