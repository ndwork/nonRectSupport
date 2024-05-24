
function run_nonRectSupport
  close all; clear; rng(1);

  datacase = 5;
  showScale = 3;

  [ supports, kDataEvenCols, kDataOddCols, reconFull ] = loadDatacase( datacase );
%load( 'junk.mat', 'supports', 'kDataEvenCols', 'kDataOddCols', 'reconFull' );

  tic;
  if datacase == 0
    recon = reconNonRectSupport_ankle( kData, supports );
  else
    recon = reconNonRectSupports( supports, kDataEvenCols, kDataOddCols );
  end
  timeTaken = toc;

  if ~ismatrix( supports)
    support = max( supports, [], 3 );
  else
    support = supports;
  end
  disp([ 'Time taken: ', num2str( timeTaken ) ]);
  figure;  imshowscale( abs( recon ) .* support, showScale );  titlenice( 'recon' );

  if numel( reconFull ) > 0
    figure;  imshowscale( abs( reconFull ) .* support, showScale );  titlenice( 'reconFull' );
    diff = reconFull - recon;
    figure;  imshowscale( abs( diff ) .* support, showScale );
    titlenice( 'diff' );  colorbarnice;
  end

  samplingBurden = ( prod( size( kDataEvenCols, [1 2] ) ) + prod( size( kDataOddCols, [1 2] ) ) ) / ...
    prod( size( support ) );   %#ok<PSIZE>
  disp([ 'Sampling burden: ', num2str(samplingBurden) ]);
end
