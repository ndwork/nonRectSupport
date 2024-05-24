
function out = applyC_1D( F, domTraj, rangeTraj, kC, C )
  % out = applyC_1D( F, domTraj, rangeTraj, N, kC, C )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F - An nTraj x nFs array representing the values of the function evaluated at each point in domTraj
  %     nFs is the number of independent values for each of the trajectory points
  %   domTraj - An nTraj x 1 array specifying the k-space coordinates of the domain trajectory points
  %   rangeTraj - An nNew x 1 array specifying the k-space coordinates of the new points
  %   kC - array of convolution kernel domain values
  %   C  - array of convolution kernel range values
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  kDistThresh = max( kC );

  nDomTraj = numel( domTraj );
  nRangeTraj = numel( rangeTraj );

  nFs = size( F, 2 );
  out = zeros( nRangeTraj, nFs );
  for domTrajIndx = 1 : nDomTraj

    dTraj_m_rangeTraj = domTraj(domTrajIndx) - rangeTraj;

    kDists = min( abs( dTraj_m_rangeTraj         ), ...
                  abs( dTraj_m_rangeTraj - 1.0 ) );
    kDists = min( kDists, ...
                  abs( dTraj_m_rangeTraj + 1.0 ) );

    shortIndxs = find( kDists < kDistThresh );
    if numel( shortIndxs ) == 0, continue; end

    CVals = interp1( kC, C, kDists( shortIndxs ), 'linear', 0 );

    FCVals = CVals * F( domTrajIndx, : );

    out( shortIndxs, : ) = out( shortIndxs, : ) + FCVals;
  end

end
