/*
 * Function name  : warpImage.cpp
 * Author name    : Vijay Rengarajan
 * Creation date  : June 5, 2015
 *
 * Syntax:
 * targetImage = warpImage(sourceImage, homography, offset, targetSize, sourceOrigin)
 *
 * sourceImage : grayscale (nrows x ncols) or colour (nrows x ncols x 3)
 *
 * homography : 3x3 matrix
 *   This has to invertible. When homography is applied on the source, we
 *   get the target image. Inside the function, we find the inverse to map
 *   from the target grid to souce grid.
 *
 * offset : offset of target wrt source.
 *   For example, if offset = [-nrows/2,-ncols/2], then there will be zeros
 *   appended to the top-half and left-half of the image. That is, the
 *   origin of the targetImage is shifted by [-nrows/2,-ncols/2].
 *
 * targetSize : size of targetImage
 *   offset and targetSize are useful in mosaicing to increase the size of
 *   the canvas. The first element is row number and second, col number.
 *
 * sourceOrigin: The (0,0) location in sourceImage
 *   This is the point around which the image will be rotated.
 *
 * Updates:
 * Created : June 5, 2015
 */

#include "mex.h"
#include<math.h>
#include "myMatrix.h"

/* This is the warping function. I use target-to-source mapping. For every
 * pixel in the target grid, I find the correponding pixel in the source
 * grid using the inverse homography. I use bilinear interpolation with
 * zero padding outside the source image region. The target grid is offset
 * by the user input, and the origin is provided by the user.
 */
void warpImage(double *source, double *target, double *origH,
        int tgt_nrows, int tgt_ncols, int src_nrows, int src_ncols,
        int nchan, double src_orig_row, double src_orig_col,
        int st_row, int st_col)
{
    // Number of pixels
    int tgt_npix, src_npix;
    
    // 2D point vars
    int src_row, src_col, f_src_row, f_src_col;
    double tgt_row, tgt_col, src_row_ni, src_col_ni;
    
    // 1D point vars
    int tgt_idx, src_idx;
    
    // Loop vars
    int i, j, k, ii, jj;
    
    // Bilinear interp vars
    int *src_col_list = new int[2];
    int *src_row_list = new int[2];
    double *wt_list = new double[4];
    double *H = new double[9];
    double wt;
    
    tgt_npix = tgt_nrows * tgt_ncols;
    src_npix = src_nrows * src_ncols;
    
    /* Calculate homography inverse for target to source mapping */
    matInverse(origH, H);
    
    /* ------------------------------------
     * For every pixel in target
     * ------------------------------------
     */
    for(i=0; i<tgt_ncols; i++) {
        for(j=0; j<tgt_nrows; j++) {
            
            /* --------------------------------------------------------
             * Homography mapping from target to source
             * --------------------------------------------------------
             */
            
            /* Get row number and col number of the target grid. Offset
             * and move the origin.
             */
            tgt_row = j + st_row - src_orig_row;
            tgt_col = i + st_col - src_orig_col;
            
            /* Apply homography to get non-int row,col of source grid */
            src_col_ni = (H[0]*tgt_col + H[3]*tgt_row + H[6]) / (H[2]*tgt_col + H[5]*tgt_row + H[8]);
            src_row_ni = (H[1]*tgt_col + H[4]*tgt_row + H[7]) / (H[2]*tgt_col + H[5]*tgt_row + H[8]);
            
            /* Uncompensate the origin to get back to the source row and
             * col numbers.
             */
            src_col_ni = src_col_ni + src_orig_col;
            src_row_ni = src_row_ni + src_orig_row;
            
            /* ---------------------------
             * Bilinear interpolation
             * ---------------------------
             */
            /* For bilinear interpolation */
            f_src_col = (int)floor(src_col_ni);
            f_src_row = (int)floor(src_row_ni);
            
            src_col_list[0] = f_src_col;
            src_col_list[1] = f_src_col+1;
            src_row_list[0] = f_src_row;
            src_row_list[1] = f_src_row+1;
            
            /* Weights for four corners.
             * The order is top-left, bottom-left, top-right, bottom-right.
             */
            wt_list[0] = (1 - src_col_ni + f_src_col) * (1 - src_row_ni + f_src_row);
            wt_list[1] = (1 - src_col_ni + f_src_col) * (src_row_ni - f_src_row);
            wt_list[2] = (src_col_ni - f_src_col) * (1 - src_row_ni + f_src_row);
            wt_list[3] = (src_col_ni - f_src_col) * (src_row_ni - f_src_row);
            
            tgt_idx = i*tgt_nrows + j;
            /* Default intensity of target is zero */
            for (k=0; k<nchan; k++) {
                target[k*tgt_npix + tgt_idx] = 0;
            }
            
            /* For every pair of row,col in the four-point set */
            for (ii=0; ii<2; ii++) {
                for (jj=0; jj<2; jj++) {
                    /* Get this row and this col number, and weight */
                    src_col = src_col_list[ii];
                    src_row = src_row_list[jj];
                    wt = wt_list[ii*2 + jj];
                    
                    /* Get 1D index of target and source points */
                    src_idx = (int)(src_col*src_nrows + src_row);
                    
                    /* Map intensity to target only if the source point is
                     * within the source image. Else, assign 0. Do for all
                     * colour channels
                     */
                    if (src_row >=0 && src_row < src_nrows && src_col >=0 && src_col < src_ncols) {
                        /* Map intensity */
                        for (k=0; k<nchan; k++) {
                            target[k*tgt_npix + tgt_idx] = target[k*tgt_npix + tgt_idx] + wt * source[k*src_npix + src_idx];
                        }
                    }
                } // end for jj
            } // end for ii
            
        }
    }
    delete[] src_col_list;
    delete[] src_row_list;
    delete[] wt_list;
    delete[] H;
}

/* Matlab calls this function */
void mexFunction(int nhls, mxArray *plhs[], int nrls, const mxArray *prhs[])
{
    int nrows, ncols, tgt_nrows, tgt_ncols, ndim, nchan, i, j;
    
    double *source, *target, *H, *off, *tgt_siz, *orig;
    const int *siz;
    int *tgt_siz_1;
    
    /* Create pointer to inputs */
    source = mxGetPr(prhs[0]); // source image
    H = mxGetPr(prhs[1]); // homography
    off = mxGetPr(prhs[2]); // offset of the target image
    tgt_siz = mxGetPr(prhs[3]); // size of the target image
    orig = mxGetPr(prhs[4]); // origin of the source grid
    
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
    
    /* Size of the target image */
    tgt_nrows = (int) tgt_siz[0];
    tgt_ncols = (int) tgt_siz[1];
    
    /* tgt_siz_1 is an int pointer, which has to be given as input to
     * mxCreateNumericArray function. The original pointer tgt_siz raises
     * an error, since it's a double pointer. Therefore, I create this
     * new int pointer which contains the size of the target image
     * (including the number of colour channels, 3) if it is colour.
     */
    if (nchan == 1) {
        tgt_siz_1 = new int[2];
        tgt_siz_1[0] = tgt_nrows;
        tgt_siz_1[1] = tgt_ncols;
    }
    else {
        tgt_siz_1 = new int[3];
        tgt_siz_1[0] = tgt_nrows;
        tgt_siz_1[1] = tgt_ncols;
        tgt_siz_1[2] = 3;
    }
    
    /* Create the output matrix */
    plhs[0] = mxCreateNumericArray(ndim,tgt_siz_1,mxDOUBLE_CLASS,mxREAL);
    
    /* Get a pointer to the output matrix */
    target = mxGetPr(plhs[0]);
    
    /* Call the function */
    warpImage(source, target, H, tgt_nrows, tgt_ncols, nrows, ncols,
            nchan, orig[0], orig[1], (int)off[0], (int)off[1]);
    
    delete[] tgt_siz_1;
}