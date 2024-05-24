
function recon = reconNonRectSupport( support, kDataEvenCols, kDataOddCols )
  % recon = reconNonRectSupport( support, kDataEvenCols, kDataOddCols )
  %
  % Inputs:
  % support - a two-dimensional binary mask with 1 indicating which pixels are in the support
  % kDataEvenCols - 
  % kDataOddCols - 
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

  sFull = size( kDataEvenCols );
  sFull(2) = sFull(2) * 2;
  outerRows = outerRowsFromSupport( support );
  innerRows = 1 - outerRows;
  nInnerRows = sum( innerRows );

  kDataEvenZF = zeros( sFull );  % ZF - zero filled
  kDataEvenZF(:,2:2:end,:) = kDataEvenCols;
  coilReconsEvenCols = fftshift2( ifft2( ifftshift2( kDataEvenZF ) ) );

  outerSupport = bsxfun( @times, support, outerRows );
  coilReconsOuter = 2 * bsxfun( @times, coilReconsEvenCols, outerSupport );

  kxFull = size2fftCoordinates( nCols );
  kxInnerMissing = kxFull(2:2:end);
  kyInner = size2fftCoordinates( nInnerRows );
  [ kxInnerMissingMesh, kyInnerMissingMesh ] = meshgrid( kxInnerMissing, kyInner );
  trajMissing = [ kyInnerMissingMesh(:) kxInnerMissingMesh(:) ];
  nxInnerMissing = numel( kxInnerMissing );
  kMissing = iGrid_2D( coilReconsEvenCols, trajMissing );
  kMissing = reshape( kMissing, [ nInnerRows nxInnerMissing nCoils ] );

  kInner = zeros( nInnerRows, nCols, nCoils );
  kInner(:,1:2:end,:) = kDataOddCols;
  kInner(:,2:2:end,:) = kMissing;

  % Now that we have kInner, we can subtract away the outer portion

  [ kxInnerMesh, kyInnerMesh ] = meshgrid( kxFull, kyInner );
  trajInner = [ kyInnerMesh(:) kxInnerMesh(:) ];

  tmp = iGrid_2D( coilReconsOuter, trajInner );
  kRemaining = kInner - reshape( tmp, [ nInnerRows nCols nCoils ] );

  coilReconsInner = grid_2D( reshape( kRemaining, [], nCoils ), trajInner, [ nRows nCols ] );
  coilReconsInner = bsxfun( @times, coilReconsInner, innerRows );

  coilRecons = coilReconsInner + coilReconsOuter;

  if nCoils > 1
    recon = mri_reconRoemer( coilRecons );
  else
    recon = coilRecons;
  end
end
