
function outerRows = outerRowsFromSupport( support )

  nRows = size( support, 1 );
  nCols = size( support, 2 );
  nCoils = size( support, 3 );
  shiftedSupport = circshift( support, floor(nCols/2), 2 );
  summedSupport = support + shiftedSupport;
  
  outerRows = zeros( nRows, nCoils );
  for coilIndx = 1 : nCoils
    maxRowsSummedSupport = max( summedSupport(:,:,coilIndx), [], 2 );
    firstOuterRowIndx = find( maxRowsSummedSupport > 1, 1, 'first' ) - 1;
    lastOuterRowIndx = find( maxRowsSummedSupport > 1, 1, 'last' ) + 1;
    
    outerRows( 1 : firstOuterRowIndx, coilIndx ) = 1;
    outerRows( lastOuterRowIndx : end, coilIndx ) = 1;
  end
end
