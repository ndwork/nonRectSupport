
function out = keepLargestRegion( in )

  [ labels, nLabels ] = bwlabel( in );
  nPerLabel = zeros( nLabels, 1 );
  for labelIndx = 1 : nLabels
    nPerLabel( labelIndx ) = ( sum( find( labels(:) == labelIndx ) ) );
  end

  [~,largestLabel] = max( nPerLabel );
  out = zeros( size( in ) );
  out( labels == largestLabel ) = 1;
end
