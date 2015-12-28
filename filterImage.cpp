/*
 * Function name  : filterImage.cpp
 * Author name    : Vijay Rengarajan
 * Creation date  : June 10, 2015
 *
 * Syntax:
 * targetImage = filterMatrix(sourceImage, filter)
 *
 * sourceImage : nrowsxncols or nrowsxncolsx3 image
 *
 * filter : nxn matrix
 *
 * Updates:
 * June 10, 2015 : Created
 */

#include "mex.h"
#include<math.h>


/* Create matrix for the given filter.
 */
void filterImage(double *source, double *target, double *h, int nrows, 
        int ncols, int nchan, int nfrows, int nfcols)
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
    int i, j, ii, jj, k;
    
    npix = nrows * ncols;
    nfrowsby2 = (nfrows - 1) / 2;
    nfcolsby2 = (nfcols - 1) / 2;
    
    /* ------------------------------------
     * For every pixel in target
     * ------------------------------------
     */
    for(i=0; i<ncols; i++) {
        for(j=0; j<nrows; j++) { 
            
            tgt_idx = i*nrows + j;
            
            /* Default intensity of target is zero */
            for (k=0; k<nchan; k++) {
                target[k*npix + tgt_idx] = 0;
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
                    if (src_row >=0 && src_row < nrows && src_col >=0 && src_col < ncols) {
                        for (k=0; k<nchan; k++) {
                            target[k*npix + tgt_idx] = target[k*npix + tgt_idx] + h[ii*nfrows + jj] * source[k*npix + src_idx];
                        }
                    }
                    
                }
            }
        }
    }
}


/* Matlab calls this function */
void mexFunction(int nhls, mxArray *plhs[], int nrls, const mxArray *prhs[])
{
    int nrows, ncols, nfrows, nfcols, ndim, nchan;
    
    double *source, *target, *h;
    const int *siz, *fsiz;
    
    /* Create pointer to inputs */
    source = mxGetPr(prhs[0]); // source image
    h = mxGetPr(prhs[1]); // filter matrix
    
    /* Size of the source image */
    ndim = mxGetNumberOfDimensions(prhs[0]); // Number of dimensions
    siz = mxGetDimensions(prhs[0]); // Each element is the size of each dim
    
    /* If the image is neither grayscale nor colour, return */
    if (ndim < 2 || ndim > 3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notImage","Not a valid image\n");
    }
    
    nrows = (int) siz[0];
    ncols = (int) siz[1];
    if (ndim == 2) { nchan = 1; }
    else { nchan = 3; }
    
    /* Size of the filter */
    fsiz = mxGetDimensions(prhs[1]); // Each element is the size of each dim
    
    nrows = (int) siz[0];
    ncols = (int) siz[1];
    
    nfrows = (int) fsiz[0];
    nfcols = (int) fsiz[1];
    
    /* Create the output matrix */
    plhs[0] = mxCreateNumericArray(ndim,siz,mxDOUBLE_CLASS,mxREAL);
    
    /* Get a pointer to the output matrix */
    target = mxGetPr(plhs[0]);
    
    /* Call the function */
    filterImage(source, target, h, nrows, ncols, nchan, nfrows, nfcols);
}