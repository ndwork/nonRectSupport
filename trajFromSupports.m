
function [ trajOdd, trajEven, outerRows ] = trajFromSupports( supports )
  % [ trajOdd, trajEven, outerRows ] = trajFromSupports( supports )
  %
  % Inputs:
  % supports - and MxNxC array of 1s and 0s that indicates support.
  %   Note: C is the number of coils.
  %
  % Outputs:
  % trajOdd - an array with each row corresponding to (ky/kx) for the odd columns of data
  % trajEven - an array with each row corresponding to (ky/kx) for the even columns of data
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It is offered without any warranty expressed
  % or implied, including the implied warranties of merchantability or fitness for a particular purpose.

  outerRows = outerRowsFromSupports( supports );
  innerRows = 1 - outerRows;
  nInnerRows = sum( innerRows );
  if numel( nInnerRows ) > 1, nInnerRows = max( nInnerRows ); end

  kyOdd = size2fftCoordinates( nInnerRows );

  ks = size2fftCoordinates( size( supports ) );
  ky = ks{1};  kx = ks{2};

  [  kxOdd,  kyOdd ] = meshgrid( kx(1:2:end), kyOdd );
  [ kxEven, kyEven ] = meshgrid( kx(2:2:end), ky );

  trajEven = [ kyEven(:), kxEven(:) ];
  trajOdd  = [  kyOdd(:),  kxOdd(:) ];

end
