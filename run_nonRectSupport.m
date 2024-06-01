
function run_nonRectSupport
  close all; clear; rng(1);

  datacase = 3;
  showScale = 3;
  mbr = true;

  %[ supports, kDataEvenCols, kDataOddCols, reconFull, sMaps ] = loadDatacase( datacase );
%save( 'junk.mat', 'supports', 'kDataEvenCols', 'kDataOddCols', 'reconFull', 'sMaps' );
load( 'junk.mat', 'supports', 'kDataEvenCols', 'kDataOddCols', 'reconFull', 'sMaps' );

  tic;
  if mbr == true
    support = max( supports, [], 3 );
    [ trajOdd, trajEven ] = trajFromSupports( supports );
    kTraj = [ trajOdd; trajEven; ];
    nCoils = size( kDataOddCols, 3 );
    kData = [ reshape( kDataOddCols, [], nCoils ); reshape( kDataEvenCols, [], nCoils ); ];
    if numel( sMaps ) == 0
      sImg = size( supports, [1 2] );
      sMaps = mri_makeSensitivityMaps( kData, 'sImg', sImg, 'traj', kTraj );
    end
    recon = mri_reconModelBased( kData, sMaps, 'traj', kTraj, 'support', support );
  else
    if datacase == 0
      recon = reconNonRectSupport_ankle( kData, supports );
    else
      recon = reconNonRectSupports( supports, kDataEvenCols, kDataOddCols );
    end
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
