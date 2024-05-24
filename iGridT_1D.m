
function out = iGridT_1D( F, traj, N, varargin )
  % out = iGridT_1D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  % Based on EE369C notes written by John Pauly
  % Detailed in the following document http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %   traj is a Mx1 array specifying the k-space locations normalized on [-0.5,0.5).
  %   N is the size of the output vector
  %
  % Optional Inputs:
  %   alpha - a float parameter specifying the oversampling factor
  %   W - an integer specifying the kernel's width
  %   nC - specifies the number of samples in the kernel
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  if nargin < 1
    disp( 'Usage:  out = iGridT_2D( F, traj, N, [ ''alpha'', alpha, ''W'', W, ''nC'', nC ] )' );
    if nargout > 0, out = []; end
    return;
  end

  checknum = @(x) numel(x) == 0 || min( isnumeric(x) & isscalar(x) & (x >= 1) ) == 1;
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  if numel( alpha ) == 0, alpha = 1.5; end

  % Make the Kaiser Bessel convolution kernel
  Gy = N;
  [kC,C,c] = makeKbKernel( Gy, N, 'alpha', alpha, 'W', W, 'nC', nC );
  Nc = numel( c );

  % Perform a circular convolution
  fftGridded = applyCT_1D( F, Nc, traj, kC, C );

  % Perform an ifft
  data = fftshift( ifft( ifftshift( fftGridded, 1 ) ), 1 );

  % Perform deapodization
  out = bsxfun( @rdivide, data, c );

  % Crop out center region if oversampling was used
  if alpha ~= 1
    sOut = size( out );
    sOut( 1 ) = N;
    out = cropData( out, sOut );
  end
end
