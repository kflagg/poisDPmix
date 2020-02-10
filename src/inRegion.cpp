// Tests if a point (x0, y0) is in the interior of the polygonal region with
// vertices (x, y) using the parity algorithm (Edelsbrunner and Harer, 2010).
// Counter needs to be initialized to zero outside of this function.
bool inRegion(const double& x0, const double& y0, double* x, double* y, const int& nV, double& x1, double& y1, double& x2, double& y2, int& counter, int& i, int& j){
  j = 0; // First point index.

  // Loop through all pairs of adjacent vertices.
  for(i = 1; i < nV; i++){
    // Check if point j is above point i and make the lower point first.
    if(y[j] > y[i]){
      x1 = x[i];
      y1 = y[i];
      x2 = x[j];
      y2 = y[j];
    } else {
      x1 = x[j];
      y1 = y[j];
      x2 = x[i];
      y2 = y[i];
    }

    // If the horizontal line through (x0, y0) passes between (x1, y1) and
    // (x2, y2) then check the determinant.
    if(!(y0 < y1 || y0 > y2)){
      // If this determinant is positive, the path (x0, y0) -> (x1, y1) ->
      // (x2, y2) turns left so (x0, y0) is to the left of the edge.
      if(x1 * y2 - x2 * y1 - x0 * y2 + x2 * y0 + x0 * y1 - x1 * y0 > 0){
        counter++;
      }
    }

    // Second point becomes first point.
    j = i;
  }
  return counter % 2;
}

// Tests if a point (x0, y0) is in the interior of a union of polygons.
bool inRegions(const double& x0, const double& y0, double** xpt, double** ypt, const int& nP, const int* nVext, double** xint, double** yint, const int* nVint, const int* offs, const int* ends, double& x1, double& y1, double& x2, double& y2, int& counter, int& i, int& j, int& polyidx, int& holeidx){
  for(polyidx = 0; polyidx < nP; polyidx++){
    // Count exterior intersections.
    if(inRegion(x0, y0, xpt[polyidx], ypt[polyidx], nVext[polyidx],
                x1, y1, x2, y2, counter, i, j)){
      // Count interior intersections.
      for(holeidx = offs[polyidx]; holeidx < ends[polyidx]; holeidx++){
        inRegion(x0, y0, xint[holeidx], yint[holeidx], nVint[holeidx],
                 x1, y1, x2, y2, counter, i, j);
      }
      polyidx = nP; // End the loop if we've found a polygon the point is in.
    }
  }
  return counter % 2;
}
