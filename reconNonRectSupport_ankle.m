
function [ recon, sMaps ] = reconNonReconSupport_ankle( support, kData )
  % [ recon, sMaps ] = reconNonReconSupport_ankle( support, kData )
  %
  % Inputs:
  % support - a 2D mask indicating the support
  % kData - an array of size M x N x nCoils where the non-zero elements correspond
  %   to the data collected.
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


  [ M, ~ ] = size( support );

  kDataSubHoriz = kData;
  kDataSubHoriz(:,1:2:end,:) = 0;
  coilReconsHoriz = fftshift2( ifft2( ifftshift2( kDataSubHoriz ) ) );
figure;  showImageCube( abs( coilReconsHoriz ) );  titlenice( 'coilReconsHoriz' );
  coilReconsHoriz = bsxfun( @times, coilReconsHoriz, support );


  coilReconsUpper = coilReconsHoriz;
  coilReconsUpper( ceil((M-1)/2)+1:end, :, : ) = 0;
figure;  showImageCube( abs( coilReconsUpper ) );  titlenice( 'coilReconsUpper' );
  kUpper = fftshift2( fft2( ifftshift2( coilReconsUpper ) ) );
  kRemaining = kData - 2 * kUpper;

coilReconsRemaining = fftshift2( ifft2( ifftshift2( kRemaining ) ) );
figure;  showImageCube( abs( coilReconsRemaining ) );  titlenice( 'coilReconsRemaining' );

  % Now we know that the remaining portion of the image has a field of view of half the
  % size in the vertical direction.
  kDataSubVert = kRemaining;
  kDataSubVert(1:2:end,:,:) = 0;
  coilReconsLower = fftshift2( ifft2( ifftshift2( kDataSubVert ) ) );
figure;  showImageCube( abs( coilReconsLower ) );  titlenice( 'coilReconsLower' );
  coilReconsLower( 1 : ceil((M-1)/2), :, : ) = 0;

  coilRecons = coilReconsLower + coilReconsUpper;
figure;  showImageCube( abs( coilRecons ) );  titlenice( 'coilRecons' );

  recon = mri_reconRoemer( coilRecons );
figure;  imshowscale( abs( recon ), 3 );

end
