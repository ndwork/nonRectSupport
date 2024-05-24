
function outerRows = outerRowsFromSupports( supports )

  nRows = size( supports, 1 );
  nCols = size( supports, 2 );
  shiftedSupports = circshift( supports, floor(nCols/2), 2 );
  summedSupports = supports + shiftedSupports;

  nCoils = size( supports, 3 );
  outerRows = zeros( nRows, nCoils );

  for coil = 1 : nCoils
    summedSupport = summedSupports(:,:,coil);
    if max( summedSupport(:) ) == 1
      outerRows(:,coil) = 1;
      continue;
    end

    maxRowsSummedSupport = max( summedSupports(:,:,coil), [], 2 );
    firstOuterRowIndx = find( maxRowsSummedSupport > 1, 1, 'first' ) - 1;
    lastOuterRowIndx = find( maxRowsSummedSupport > 1, 1, 'last' ) + 1;
    
    outerRows( 1 : firstOuterRowIndx, coil ) = 1;
    outerRows( lastOuterRowIndx : end, coil ) = 1;
  end
end
