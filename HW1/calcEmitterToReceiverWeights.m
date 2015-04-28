function [rowVec] = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols)

%assumes y axis points down
delta   = round([cos(direction), sin(direction)]);

rowVec = zeros(1, A_cols);
curPos = startXY;
while (1),
    idx = sub2ind([X_rows, X_cols], curPos(2), curPos(1));
    rowVec(idx) = norm(delta);
    disp(['Passing through cell: (' num2str(curPos(2)) ', ' num2str(curPos(1)) ')']);
    next = curPos + delta;
    curPos = next;
    if ( ~(BoundingBox(1) <= next(1) && next(1) <= BoundingBox(2) && ...
          BoundingBox(3) <= next(2) && next(2) <= BoundingBox(4)))
     break;
    end
    
end
