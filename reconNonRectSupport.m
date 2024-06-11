
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

  [ trajOdd, trajEven, outerRows ] = trajFromSupports( support );
  nEven = size( trajEven, 1 );

  %-- Perform the reconstruction

  % 1) Reconstruct with even columns; this results in aliased images
  kDataEvenZF = zeros( [ nRows nCols nCoils ] );  % ZF - zero filled
  kDataEvenZF(:,2:2:end,:) = reshape( kDataEvenCols, size( kDataEvenZF(:,2:2:end,:) ) );
  reconEvenCols = fftshift2( ifft2( ifftshift2( kDataEvenZF ) ) );

  % 2) Reconstruct outer image
  outerSupport = bsxfun( @times, support, outerRows );
  reconsOuter = 2 * bsxfun( @times, reconEvenCols, outerSupport );

  % 3) Interpolate aliased images onto odd column trajectory with nInnerRows in each column
  kOuter = iGrid_2D( reconsOuter, [ trajEven; trajOdd; ] );
  kOuterEven = reshape( kOuter( 1 : nEven, : ), size( kDataEvenCols ) );
  kOuterOdd = reshape( kOuter( nEven + 1 : end, : ), size( kDataOddCols ) );

  % 4) Subtract outer Fourier values from full Fourier values
  kRemainingEven = reshape( kDataEvenCols - kOuterEven, [], nCoils );
  kRemainingOdd = reshape( kDataOddCols - kOuterOdd, [], nCoils );

  % 5) Reconstruct the inner region
  reconsInner = grid_2D( [ kRemainingEven; kRemainingOdd; ], [ trajEven; trajOdd; ], [ nRows nCols ] );
  reconsInner = bsxfun( @times, reconsInner, 1 - outerRows );

  % 6) Sum to create the whole image
  recons = reconsInner + reconsOuter;

  if nCoils > 1
    recon = mri_reconRoemer( recons ) .* support;
  else
    recon = recons .* support;
  end
end

