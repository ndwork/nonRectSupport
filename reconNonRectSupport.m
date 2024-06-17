
function recon = reconNonRectSupport( support, kDataEvenCols, kDataOddCols )
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
  nCoils = size( kDataEvenCols, 3 );
  kxFull = size2fftCoordinates( nCols );

  [ trajOdd, trajEven, outerRows ] = trajFromSupports( support );
  nEven = size( trajEven, 1 );
  nOdd = size( trajOdd, 1 );

  %-- Perform the reconstruction

  % Create the trajectory of the inner grid
  nInner = sum( 1 - outerRows );
  kInner = size2fftCoordinates([ nInner nCols ]);
  kyInner = kInner{1};  kxInner = kInner{2};
  [ kxInnerMesh, kyInnerMesh ] = meshgrid( kxInner, kyInner );
  trajInner = [ kyInnerMesh(:) kxInnerMesh(:) ];

  % 1) Reconstruct with even columns; this results in aliased images
  kDataEvenZF = zeros( [ nRows nCols nCoils ] );  % ZF - zero filled
  kDataEvenZF(:,2:2:end,:) = reshape( kDataEvenCols, size( kDataEvenZF(:,2:2:end,:) ) );
  reconEvenCols = fftshift2( ifft2( ifftshift2( kDataEvenZF ) ) );

  % 2) Reconstruct outer image
  outerSupport = bsxfun( @times, support, outerRows );
  reconsOuter = 2 * bsxfun( @times, reconEvenCols, outerSupport );

  % ) Interpolate onto the inner grid
  kInner = zeros( nInner, nCols, nCoils );
  kInner( :, 1:2:end, : ) = kDataOddCols;
  [ kxMissing, kyMissing ] = meshgrid( kxFull(2:2:end), kyInner );  % missing from inner grid
  trajMissing = [ kyMissing(:) kxMissing(:) ];
  kMissing = iGrid_2D( reconEvenCols, trajMissing );
  kInner( :, 2:2:end, : ) = reshape( kMissing, nInner, floor(nCols/2), nCoils );

  % 4) Subtract outer Fourier values from full Fourier values
  kOuter = iGrid_2D( reconsOuter, trajInner );
  kOuter = reshape( kOuter, nInner, nCols, nCoils );
  kRemaining = kInner - kOuter;

  % 5) Reconstruct the inner region
  reconsInner = grid_2D( reshape( kRemaining, [], nCoils ), trajInner, [ nRows nCols ] );
  %reconsInner = grid_2D( [ kRemainingEven; kRemainingOdd; ], [ trajEven; trajOdd; ], [ nRows nCols ] );
  reconsInner = bsxfun( @times, reconsInner, 1 - outerRows );

  % 6) Sum to create the whole image
  recons = reconsInner + reconsOuter;

  if nCoils > 1
    recon = mri_reconRoemer( recons ) .* support;
  else
    recon = recons .* support;
  end
end

