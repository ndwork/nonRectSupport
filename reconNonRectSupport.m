
function recon = reconNonRectSupport( support, kDataEvenCols, kDataOddCols, varargin )
  % recon = reconNonRectSupport( support, kDataEvenCols, kDataOddCols )
  %
  % Inputs:
  % support - a 2D binary mask with 1 indicating which pixels are in the support
  % kDataEvenCols - a 2D array representing the Fourier values of even columns
  % kDataOddCols - a 2D array representing the Fourier values of odd columns
  %
  % Outputs:
  % recon - a two dimensional image
  %
  % Written by Nicholas Dwork - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nRows = size( support, 1 );
  nCols = size( support, 2 );

  sFull = size( kDataEvenCols );
  sFull(2) = sFull(2) * 2;
  outerRows = outerRowsFromSupport( support );
  innerRows = 1 - outerRows;
  nInnerRows = sum( innerRows );

  %-- Find the relevant trajectories

  kxFull = size2fftCoordinates( nCols );
  kxInnerMissing = kxFull(2:2:end);
  kyInner = size2fftCoordinates( nInnerRows );
  [ kxInnerMissingMesh, kyInnerMissingMesh ] = meshgrid( kxInnerMissing, kyInner );
  trajMissing = [ kyInnerMissingMesh(:) kxInnerMissingMesh(:) ];

  [ kxInnerMesh, kyInnerMesh ] = meshgrid( kxFull, kyInner );
  trajInner = [ kyInnerMesh(:) kxInnerMesh(:) ];


  %-- Perform the reconstruction

  % 1) Reconstruct with even columns; this results in aliased images
  kDataEvenZF = zeros( sFull );  % ZF - zero filled
  kDataEvenZF(:,2:2:end,:) = kDataEvenCols;
  reconEvenCols = fftshift2( ifft2( ifftshift2( kDataEvenZF ) ) );

  % 2) Interpolate aliased images onto odd column trajectory with nInnerRows in each column
  kMissing = iGrid_2D( reconEvenCols, trajMissing );
  nxInnerMissing = numel( kxInnerMissing );
  kMissing = reshape( kMissing, [ nInnerRows nxInnerMissing ] );

  kInner = zeros( nInnerRows, nCols );
  kInner( :, 1:2:end ) = kDataOddCols;
  kInner( :, 2:2:end ) = kMissing;

  % 3) Reconstruct outer image
  outerSupport = bsxfun( @times, support, outerRows );
  reconsOuter = 2 * bsxfun( @times, reconEvenCols, outerSupport );

  % 4) Subtract outer image from inner trajectory
  tmp = iGrid_2D( reconsOuter, trajInner );
  kRemaining = kInner - reshape( tmp, [ nInnerRows nCols ] );

  % 5) Reconstruct the inner image
  reconsInner = grid_2D( kRemaining(:), trajInner, [ nRows nCols ] );
  reconsInner = bsxfun( @times, reconsInner, innerRows );

  % 6) Sum to create the whole image
  recon = reconsInner + reconsOuter;
end

