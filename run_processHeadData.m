
function run_processHeadData
  close all; clear; rng(1);

  dataDir = '/Users/nicholasdwork/Documents/Data/VP2.001.2_HEAD/TWIX';
  %slice2Recon = 155;  %91, 95
  slice2Recon = 153;  %92
  %slice2Recon = 150;  %93
  %slice2Recon = 145;  %99

  slice2Recon = 153;
% slice2Recon = 145;
% slice2Recon = 135;
% slice2Recon = 125;
% slice2Recon = 115;

  % Fully sampled image - meas_MID00085_FID11343_MPRAGE
  twix = mapVBVD( [ dataDir, '/meas_MID00085_FID11343_MPRAGE' ] );
  twix = twix{2};
  kDataFull = squeeze( twix.image() );
  kDataFull = permute( kDataFull, [3 4 1 2 5] );
  kDataFull = mri_reconIFFT( kDataFull, 'dims', 3 );
  kDataFull = squeeze( kDataFull(:,:,slice2Recon,:,:) );

  kDataFull = kDataFull ./ kDataFull( 65,41,1,1 );
  fullRecon = mri_reconRoemer( mri_reconIFFT( kDataFull(:,:,:,1) ) );
  figure;  imshowscale( abs( fullRecon ), 3 );

  % Under sampled image - meas_MID00087_FID11345_MPRAGE_WITH_WRAP
  twix = mapVBVD( [ dataDir, '/meas_MID00087_FID11345_MPRAGE_WITH_WRAP' ] );
  twix = twix{2};
  kDataWrap = squeeze( twix.image() );
  kDataWrap = permute( kDataWrap, [3 4 1 2 5] );
  kDataWrap = mri_reconIFFT( kDataWrap, 'dims', 3 );
  kDataWrap = squeeze( kDataWrap(:,:,slice2Recon,:,:) );

%load( 'junk.mat', 'kDataFull', 'kDataWrap' );


kDataFull = kDataFull ./ kDataFull( 65,41,1,2 );
kDataWrap = kDataWrap ./ kDataWrap( 47,41,1,2 );


  reconFullCoils = mri_reconIFFT( kDataFull(:,:,:,1), 'dims', [1 2 3] );
  reconFull = abs( mri_reconRoemer( reconFullCoils ) );
  noise = reconFull( 110:end, 1:20, : );
  support = reconFull > ( mean( noise(:) ) + 5 * std( noise(:) ) );

  support = imerode( support, strel( "disk", 3 ) );
  support = imdilate( support, strel( "disk", 8 ) );
  support = imerode( support, strel( "disk", 5 ) );

  kDataEvenCols = squeeze( kDataFull( :, 2:2:end, :, 2 ) );
  kDataOddCols = squeeze( kDataWrap( :, 1:2:end, :, 2 ) );

kDataEvenColsPro = squeeze( kDataFull( :, 2:2:end, :, 2 ) );
kDataOddColsPro = squeeze( kDataWrap( :, 1:2:end, :, 2 ) );


      kDataFull = kDataFull(:,:,:,2);
      M = size( kDataFull, 1 );
      N = size( kDataFull, 2 );
      nCoils = size( kDataFull, 3 );

      reconFullCoils = mri_reconIFFT( kDataFull, 'multiSlice', true );
      reconFull = abs( mri_reconRoemer( reconFullCoils ) );
      noise = reconFull(1:5,1:5);
      support = reconFull > ( mean( noise(:) ) + 5 * std( noise(:) ) );

      support = imerode( support, strel( "disk", 2 ) );
      support = imdilate( support, strel( "disk", 8 ) );
      support = imerode( support, strel( "disk", 6 ) );

      kDataEvenCols = kDataFull(:,2:2:end,:);
      kxFull = size2fftCoordinates( N );

      outerRows = outerRowsFromSupport( support );
      innerRows = 1 - outerRows;
      nInnerRows = sum( innerRows );

      kyOdd = size2fftCoordinates( nInnerRows );
      [ kxOddMesh, kyOddMesh ] = meshgrid( kxFull, kyOdd );
      trajOdd = [ kyOddMesh(:) kxOddMesh(:) ];
      kDataOddCols = cell( 1, 1, nCoils );
      parfor coilIndx = 1 : nCoils
        tmp = iGrid_2D( reconFullCoils(:,:,coilIndx), trajOdd );
        kDataOddCols{ 1, 1, coilIndx } = reshape( tmp, [ nInnerRows, N ] );
      end
      kDataOddCols = cell2mat( kDataOddCols );
      kDataOddCols = kDataOddCols(:,1:2:end,:);

  %recon = reconNonRectSupport_new( support, kDataEvenColsPro, kDataOddColsProFixed );
  recon = reconNonRectSupport( support, kDataEvenCols, kDataOddCols );

  figure;  imshowscale( abs( recon ), 3 );

end

