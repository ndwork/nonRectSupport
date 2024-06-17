
function [ trajOdd, trajEven, outerRows ] = trajFromSupports( supports, varargin )
  % [ trajOdd, trajEven, outerRows ] = trajFromSupports( supports [, downsample, downsampleEven ] )
  %
  % Inputs:
  % supports - and MxNxC array of 1s and 0s that indicates support.
  %   Note: C is the number of coils.
  %
  % Optional Inputs:
  % downsample - the amount to downsample the columns of data
  % downsampleEven - by default, assumes the same downsampling for even and odd columns
  %   If downsampleEven is supplied, then this number is used for the even columns only.
  %
  % Outputs:
  % trajOdd - an array with each row corresponding to (ky/kx) for the odd columns of data
  % trajEven - an array with each row corresponding to (ky/kx) for the even columns of data
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It is offered without any warranty expressed
  % or implied, including the implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage: [ trajOdd, trajEven, outerRows ] = trajFromSupports( supports [, downsample ] )' );
    if nargout > 0
      trajOdd = []; trajEven = []; outerRows = [];
    end
    return
  end

  p = inputParser;
  p.addOptional( 'downsample', 1, @ispositive );
  p.addOptional( 'downsampleEven', [], @ispositive );
  p.parse( varargin{:} );
  dsO = p.Results.downsample;  % downsample for odd columns
  dsE = p.Results.downsampleEven;         % downsample for even columns

  if numel( dsE ) == 0, dsE = dsO; end

  outerRows = outerRowsFromSupports( supports );
  innerRows = 1 - outerRows;
  nInnerRows = sum( innerRows );
  if numel( nInnerRows ) > 1, nInnerRows = max( nInnerRows ); end

  kyOdd = size2fftCoordinates( nInnerRows );

  ks = size2fftCoordinates( size( supports ) );
  ky = ks{1};  kx = ks{2};

  [  kxOdd,  kyOdd ] = meshgrid( kx(1:2:end), kyOdd( ceil(dsO/2) : dsO : end, : ) );
  [ kxEven, kyEven ] = meshgrid( kx(2:2:end), ky( ceil(dsE/2) : dsE : end, : ) );

  trajEven = [ kyEven(:), kxEven(:) ];
  trajOdd  = [  kyOdd(:),  kxOdd(:) ];

end
