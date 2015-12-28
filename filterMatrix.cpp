/*
 * Function name  : filterMatrix.cpp
 * Author name    : Vijay Rengarajan
 * Creation date  : June 9, 2015
 *
 * Syntax:
 * targetMatrix = filterMatrix(sourceSize, filter)
 *
 * sourceSize : [nrows, ncols]
 *   Size of the source image.
 *
 * filter : nxn matrix
 * 
 * Updates:
 * June 9, 2015 : Created (for 3x3 filter)
 * June 10, 2015 : Updated for any nxn filter
 */

#include "mex.h"
#include<math.h>


/* Create matrix for the given filter.
 */
void filterMatrix(double *target, double *h, int nrows, int ncols, int nfrows, int nfcols)
{
    // Number of pixels
    int npix;
    
    // 2D point vars
    int src_row, src_col;
    
    // 1D point vars
    int tgt_idx, src_idx;
    
    // Filter vars
    int nfrowsby2, nfcolsby2;
    
    // Loop vars
    int i, j, ii, jj;
    double okFlag;
    
    okFlag = 0;
    npix = nrows * ncols;
    nfrowsby2 = (nfrows - 1 ) / 2;
    nfcolsby2 = (nfcols - 1 ) / 2;
    
    /* ------------------------------------
     * For every pixel in target
     * ------------------------------------
     */
    for(i=0; i<ncols; i++) {
        for(j=0; j<nrows; j++) {
            
            tgt_idx = i*nrows + j; 
            
            if (i - nfcolsby2 + 0 >= 0 &&  j - nfrowsby2 + 0 >= 0 && i - nfcolsby2 + nfcols < ncols && j - nfrowsby2 + nfrows < nrows) {
                okFlag = 1;
            }
            else {
                okFlag = 0;
            }
            
            for(ii=0; ii<nfcols; ii++) {
                for(jj=0;jj<nfrows; jj++) {
                    
                    /* Calculate the source index position for this element
                     * of the filter centred around ii,jj
                     */
                    src_col = i - nfcolsby2 + ii;
                    src_row = j - nfrowsby2 + jj;
                    src_idx = src_col*nrows + src_row;
                    
                    /* Add to target intensity the weighted source intensity*/
                    if (okFlag == 1 && src_row >=0 && src_row < nrows && src_col >=0 && src_col < ncols) {
                        target[(ii*nfrows + jj) * npix + tgt_idx ] = src_idx; 
                        target[(ii*nfrows + jj + nfrows*nfcols) * npix + tgt_idx] = h[ii*nfrows + jj];
                    }
                    else {
                        target[(ii*nfrows + jj) * npix + tgt_idx ] = 0; 
                        target[(ii*nfrows + jj + nfrows*nfcols) * npix + tgt_idx] = 0;
                    }
                }
            }
        }
    }
}
            

/* Matlab calls this function */
void mexFunction(int nhls, mxArray *plhs[], int nrls, const mxArray *prhs[])
{
    int nrows, ncols, nfrows, nfcols;
    
    double *siz, *target, *h;
    const int *fsiz;
    
    /* Create pointer to inputs */
    siz = mxGetPr(prhs[0]); // source image size
    h = mxGetPr(prhs[1]); // filter array
    
    nrows = (int) siz[0];
    ncols = (int) siz[1];
    
    /* Size of the filter */
    fsiz = mxGetDimensions(prhs[1]); // Each element is the size of each dim
    
    nfrows = (int) fsiz[0];
    nfcols = (int) fsiz[1];
    
    /* Create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(nrows*ncols,2*nfrows*nfcols,mxREAL);
    
    /* Get a pointer to the output matrix */
    target = mxGetPr(plhs[0]);
    
    /* Call the function */
    filterMatrix(target, h, nrows, ncols, nfrows, nfcols);
}