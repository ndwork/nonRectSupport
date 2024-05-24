

function run_sampleSpacings
  close all; clear; rng(1);

  dataDir = '/Users/nicholasdwork/Documents/Data/';
  outDir = './sampleSpacings';
  if ~exist( outDir, 'dir' ), mkdir( outDir ); end

  load( [ dataDir, './espiritData/brain_8ch.mat' ], 'DATA' );
  kDataFull = squeeze( permute( DATA, [1 2 4 3] ) );


  reconFullCoils = mri_ifftRecon( kDataFull );
  reconFull = abs( mri_reconRoemer( reconFullCoils ) );
  reconFull = reconFull / max( abs( reconFull(:) ) );

  ky = size2fftCoordinates( size(reconFull,1) );

  spacings = 1 : 0.2 : 2;
  parfor spacingIndx = 1 : numel( spacings )
    spacing = spacings( spacingIndx );

    kx = size2fftCoordinates( round( size(reconFull,2) / spacing ) );

    [kys,kxs] = ndgrid( ky, kx );
    traj = [ kys(:) kxs(:) ];

    kValues = iGrid_2D( reconFull, traj );

    recon = grid_2D( kValues, traj, size(reconFull) );
    imwrite( abs(recon), [outDir, filesep, 'spacing_', num2str(spacing, '%3.2f'), '.png' ] );

  end


end

