

function makeMasks
  close all; clear; clc;

  showScale = 3;

  % support mask 
  supportMask = zeros(256);  supportMask (2:2:end,:)=1;  supportMask (:,2:2:end)=1;
  figure;  imshowscale( supportMask, showScale );  titlenice( 'supportMask' );

  % auto-calibration region
  acr = zeros(256);  acr(128-14:128+15,128-14:128+15)=1;

  % partial Fourier mask
  pfMask = ones(256);  pfMask(:,129:end)=0;
  figure;  imshowscale( pfMask | acr, showScale );  titlenice( 'pfMask' );

  % grappa mask with acceleration factor of 3
  grappaMask = zeros(256);  grappaMask(2:6:end,:) = 1;
  figure;  imshowscale( grappaMask | acr, showScale );  titlenice( 'grappaMask' );

  % compressed sensing mask
  vdMask = vdSampleMask( [ 256 256 ], 75, 16384 );
  figure;  imshowscale( vdMask | acr, showScale );  titlenice( 'vdMask' );


  % support and pf
  supportPfMask = ( supportMask & pfMask ) | acr;
  figure;  imshowscale( supportPfMask, showScale );  titlenice( 'supportPfMask' );

  % support and grappa
  grappaComplementRows = zeros(256);
  grappaComplementRows(4:6:end,:) = 1;
  grappaComplementRows(6:6:end,:) = 1;
  supportGrappaMask = ( supportMask - grappaComplementRows ) | acr;
  figure;  imshowscale( supportGrappaMask, showScale );  titlenice( 'supportGrappaMask' );

  % support cs mask
  supportCsMask = ( supportMask & vdMask ) | acr;
  figure;  imshowscale( supportCsMask, showScale );  titlenice( 'supportCsMask' );

  % all mask
  allMask = supportPfMask & supportGrappaMask & supportCsMask;
  figure;  imshowscale( allMask, showScale );  titlenice( 'allMask' );

  disp([ 'All sampling percentage: ', num2str( sum( allMask(:) ) / numel( allMask ) ) ]);

end

