
function [ support, kDataEvenCols, kDataOddCols, reconFull ] = loadDatacase( datacase )

  dataDir = '/Volumes/NDWORK128GB/';
  if ~exist( dataDir, 'dir' )
    dataDir = '/Volumes/Data Drive/';
  end
  if ~exist( dataDir, 'dir' )
    dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/';
  end
  if ~exist( dataDir, 'dir' )
    dataDir = '/Users/ndwork/Documents/Data/';
  end
  if ~exist( dataDir, 'dir' )
    dataDir = '/Users/nicholasdwork/Documents/Data/';
  end
  if ~exist( dataDir, 'dir' ), error( 'dataDir not found' ); end

  reconFull = [];

  switch datacase

    case 0
      % Ankle
      [kDataFull,header] = read_MR_rawdata( [ dataDir, '/fullySampled3D/body3d/P27648.7' ] );   %#ok<ASGLU>
      kDataFull = squeeze( kDataFull );
      kDataFull = ifft( ifftshift( kDataFull, 3 ), [], 3 );
      kDataFull = squeeze( kDataFull(:,:,50,:) );

      for coilIndx = 1 : size( kDataFull, 3 )
        kDataFull(:,:,coilIndx) = rot90( kDataFull(:,:,coilIndx), -1 );
      end

      M = size( kDataFull, 1 );
      N = size( kDataFull, 2 );
      support = ones( M, N );
      support( 1 : ceil((M-1)/2), ceil((N-1)/2) : end ) = 0;

      %reconFull = mri_reconRoemer( mri_reconIFFT( kDataFull ) );

      kData = zeros( size( kDataFull ) );
      kData(2:2:end,:,:) = kDataFull(2:2:end,:,:);
      kData(:,2:2:end,:) = kDataFull(:,2:2:end,:);
      recon = reconNonRectSupport_ankle( support, kData );
      error( 'This is the special algorithm for the ankle.' );

    case 1
      % Ankle
      [kDataFull,header] = read_MR_rawdata( [ dataDir, '/fullySampled3D/body3d/P27648.7' ] );   %#ok<ASGLU>
      kDataFull = squeeze( kDataFull );
      kDataFull = ifft( ifftshift( kDataFull, 3 ), [], 3 );
      kDataFull = squeeze( kDataFull(:,:,50,:) );
      kDataFull = rot90( kDataFull, -1 );

      %reconFull = mri_reconRoemer( mri_reconIFFT( kDataFull, 'multiSlice', true ) );

      M = size( kDataFull, 1 );
      N = size( kDataFull, 2 );
      nCoils = size( kDataFull, 3 );
      support = ones( M, N );
      support( 1 : ceil((M-1)/2), ceil((N-1)/2) : end ) = 0;

      kDataEvenCols = kDataFull(:,2:2:end,:);
      kxFull = size2fftCoordinates( N );

      outerRows = outerRowsFromSupport( support );
      innerRows = 1 - outerRows;
      nInnerRows = sum( innerRows );

      coilReconsFull = fftshift2( ifft2( ifftshift2( kDataFull ) ) );

      kyOdd = size2fftCoordinates( nInnerRows );
      [ kxOddMesh, kyOddMesh ] = meshgrid( kxFull, kyOdd );
      trajOdd = [ kyOddMesh(:) kxOddMesh(:) ];
      kDataOddCols = cell( 1, 1, nCoils );
      parfor coilIndx = 1 : nCoils
        tmp = iGrid_2D( coilReconsFull(:,:,coilIndx), trajOdd );
        kDataOddCols{ 1, 1, coilIndx } = reshape( tmp, [ nInnerRows, N ] );
      end
      kDataOddCols = cell2mat( kDataOddCols );
      kDataOddCols = kDataOddCols(:,1:2:end,:);
      disp( 'Data loaded for case 1' );

    case 2
      % ESPIRIT Brain
      load( [ dataDir, './espiritData/brain_8ch.mat' ], 'DATA' );
      kDataFull = squeeze( permute( DATA, [1 2 4 3] ) );

      M = size( kDataFull, 1 );
      N = size( kDataFull, 2 );
      nCoils = size( kDataFull, 3 );

      reconFullCoils = mri_reconIFFT( kDataFull, 'multiSlice', true );
      reconFull = abs( mri_reconRoemer( reconFullCoils ) );
      noise = reconFull(1:25,1:25);
      support = reconFull > ( mean( noise(:) ) + 2.5 * std( noise(:) ) );

      support = imerode( support, strel( "disk", 2 ) );
      support = imdilate( support, strel( "disk", 6 ) );
      support = imerode( support, strel( "disk", 4 ) );

      kDataEvenCols = kDataFull(:,2:2:end,:);
      kxFull = size2fftCoordinates( N );

      outerRows = outerRowsFromSupport( support );
      innerRows = 1 - outerRows;
      nInnerRows = sum( innerRows );

      coilReconsFull = fftshift2( ifft2( ifftshift2( kDataFull ) ) );

      kyOdd = size2fftCoordinates( nInnerRows );
      [ kxOddMesh, kyOddMesh ] = meshgrid( kxFull, kyOdd );
      trajOdd = [ kyOddMesh(:) kxOddMesh(:) ];
      kDataOddCols = cell( 1, 1, nCoils );
      parfor coilIndx = 1 : nCoils
        tmp = iGrid_2D( coilReconsFull(:,:,coilIndx), trajOdd );
        kDataOddCols{ 1, 1, coilIndx } = reshape( tmp, [ nInnerRows, N ] );
      end
      kDataOddCols = cell2mat( kDataOddCols );
      kDataOddCols = kDataOddCols(:,1:2:end,:);

    case 3
      % UCSF Brain (my brain -Nicholas Dwork)
      [kDataFull,header] = read_MR_rawdata( [ dataDir, '/fullySampled3D/brain3d/P73216.7' ] );   %#ok<ASGLU>
      kDataFull = squeeze( kDataFull );
      kDataFull = ifft( ifftshift( kDataFull, 3 ), [], 3 );
      kDataFull = rot90( squeeze( kDataFull(:,:,80,:) ), -1 );

      reconFullCoils = mri_reconIFFT( kDataFull, 'multiSlice', true );
      reconFullCoils = reconFullCoils( 34:222, 53:210, : );
      kDataFull = fftshift2( fft2( ifftshift2( reconFullCoils ) ) );
      M = size( kDataFull, 1 );
      N = size( kDataFull, 2 );
      nCoils = size( kDataFull, 3 );

      reconFull = abs( mri_reconRoemer( reconFullCoils ) );
      noise = reconFull(1:5,1:5);
      support = reconFull > ( mean( noise(:) ) + 5 * std( noise(:) ) );

      support = imerode( support, strel( "disk", 2 ) );
      support = imdilate( support, strel( "disk", 8 ) );
      support = imerode( support, strel( "disk", 4 ) );

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

    case 4
      %pineapple
      sliceIndx = 252;  % 96!!!
      dataFile = [ dataDir, '/pineapple3T/raw_004.list' ];
      kDataFull = loadPhilipsKSpace( dataFile );

      kSlices = fftshift( ifft( ifftshift( kDataFull, 1 ), [], 1 ), 1 );
      kDataFull = squeeze( kSlices( sliceIndx, :, :, : ) );
      clear kSlices;

      imgsFull = ifft2( fftshift2( kDataFull ) );
      kDataFull = fftshift2( fft2( ifftshift2( imgsFull ) ) );

      dataFile = [ dataDir, '/pineapple3T/raw_007.list' ];
      kDataUnder = loadPhilipsKSpace( dataFile );

      kSlices = fftshift( ifft( ifftshift( kDataUnder, 1 ), [], 1 ), 1 );
      kDataUnder = squeeze( kSlices( sliceIndx, :, :, : ) );
      clear kSlices;

      imgsUnder = ifft2( ifftshift2( kDataUnder ) );
      kDataUnder = fftshift2( fft2( ifftshift2( imgsUnder ) ) );

      kDataEvenCols = kDataFull(:,2:2:end,:);
      kDataOddCols = kDataUnder(:,1:2:end,:);

      reconFullCoils = mri_reconIFFT( kDataFull, 'multiSlice', true );
      reconFull = abs( mri_reconRoemer( reconFullCoils ) );
      noise = reconFull(1:5,1:5);

      support = reconFull > ( mean( noise(:) ) + 5 * std( noise(:) ) );
      labels = bwlabel( 1 - support );
      support(:) = 1;
      support( labels==labels(1,1) ) = 0;
      support( labels==labels(end,1) ) = 0;
      support( labels==labels(end,end) ) = 0;
      support( labels==labels(1,end) ) = 0;

      support = imerode( support, strel( "disk", 2 ) );
      support = imdilate( support, strel( "disk", 6 ) );
      support = imerode( support, strel( "disk", 3 ) );
      disp( 'Pineapple data loaded.' );

    case 5
      % Ankle - parallel imaging, multiple supports
      [kDataFull,header] = read_MR_rawdata( [ dataDir, '/fullySampled3D/body3d/P27648.7' ] );   %#ok<ASGLU>
      kDataFull = squeeze( kDataFull );
      kDataFull = ifft( ifftshift( kDataFull, 3 ), [], 3 );
      kDataFull = squeeze( kDataFull(:,:,50,:) );
      kDataFull = rot90( squeeze( kDataFull ), -1 );
      kDataFull = kDataFull( :, :, [ 1 3 5 6 7 8] );
      coilReconsFull = mri_reconIFFT( kDataFull, 'multiSlice', true );
      reconFull = mri_reconRoemer( coilReconsFull );

      noise = coilReconsFull( 1 : 50, end-50 : end, : );

      nCoils = size( kDataFull, 3 );

      % Come up with independent supports for each coil
      supports = zeros( size( coilReconsFull ) );

      coilSupport = 20*log10( abs( coilReconsFull(:,:,1) ) ) > -47;
      coilSupport = keepLargestRegion( coilSupport );
      coilSupport(1:129,:) = 0;
      coilSupport(141:155,12:16) = 1;
      coilSupport(155:165,10:20) = 1;
      coilSupport( bwlabel( 1 - coilSupport ) == 46 ) = 1;
      coilSupport( bwlabel( 1 - coilSupport ) == 59 ) = 1;
      coilSupport( bwlabel( 1 - coilSupport ) ==  2 ) = 1;
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      coilSupport = imerode( coilSupport, strel( "disk", 3 ) );
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      supports(:,:,1) = coilSupport;
      supports(:,:,2) = coilSupport;

      coilSupport = 20*log10( abs( coilReconsFull(:,:,3) ) ) > -42;
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      coilSupport = imerode( coilSupport, strel( "disk", 3 ) );
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      coilSupport = keepLargestRegion( coilSupport );
      labels = bwlabel( 1 - coilSupport );
      coilSupport( labels == 2 ) = 1;
      coilSupport( labels == 3 ) = 1;
      coilSupport( labels == 4 ) = 1;
      coilSupport( labels == 5 ) = 1;
      supports(:,:,3) = coilSupport;

      coilSupport = abs( coilReconsFull(:,:,4) ) > 0.012;
      coilSupport(91,15:25) = 1;
      coilSupport = keepLargestRegion( coilSupport );
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      coilSupport = imerode( coilSupport, strel( "disk", 3 ) );
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      supports(:,:,4) = coilSupport;

      coilSupport = abs( coilReconsFull(:,:,5) ) > 0.013;
      coilSupport = imerode( coilSupport, strel( "disk", 2 ) );
      coilSupport = imdilate( coilSupport, strel( "disk", 6 ) );
      coilSupport = keepLargestRegion( coilSupport );
      supports(:,:,5) = coilSupport;

      coilSupport = abs( coilReconsFull(:,:,6) ) > 0.008;
      coilSupport = keepLargestRegion( coilSupport );
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      coilSupport = imerode( coilSupport, strel( "disk", 4 ) );
      coilSupport = imdilate( coilSupport, strel( "disk", 2 ) );
      labels = bwlabel( 1 - coilSupport );
      coilSupport( labels == 2 ) = 1;
      coilSupport( labels == 3 ) = 1;
      coilSupport( labels == 4 ) = 1;
      supports(:,:,6) = coilSupport;

      nCoils = size( kDataFull, 3 );
      kDataEvenCols = kDataFull(:,2:2:end,:);

      N = size( kDataFull, 2 );
      kxFull = size2fftCoordinates( N );

      outerRows = outerRowsFromSupports( supports );
      innerRows = 1 - outerRows;
      nInnerRows = max( sum( innerRows, 1 ) );

      kyOdd = size2fftCoordinates( nInnerRows );
      [ kxOddMesh, kyOddMesh ] = meshgrid( kxFull, kyOdd );
      trajOdd = [ kyOddMesh(:) kxOddMesh(:) ];
      kDataOddCols = cell( 1, 1, nCoils );
      parfor coilIndx = 1 : nCoils
        tmp = iGrid_2D( coilReconsFull(:,:,coilIndx), trajOdd );
        kDataOddCols{ 1, 1, coilIndx } = reshape( tmp, [ nInnerRows, N ] );
      end
      kDataOddCols = cell2mat( kDataOddCols );
      kDataOddCols = kDataOddCols(:,1:2:end,:);

      support = supports;

    otherwise
      error( 'This datacase doesn''t exist' );
  end

end

