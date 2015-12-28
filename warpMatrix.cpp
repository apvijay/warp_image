/*
 * Function name  : mbMatrix.cpp
 * Author name    : Vijay Rengarajan
 * Creation date  : June 9, 2015
 *
 * Syntax:
 * targetImage = mbMatrix(sourceSize, homography, offset, targetSize, sourceOrigin)
 *
 * sourceSize : [nrows, ncols]
 *   Size of the source image.
 *
 * homographies : 3x3 matrix
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
 * targetMatrix : (nrows*ncols)x8 matrix
 *   Returns a matrix which contains source indices and weights which have 
 *   to be mapped for every target index. Every row is for a target index. 
 *   The first four elements of a column correspond to the source 
 *   indices and the next four correspond to their weights (according to 
 *   bilinear interpolation). By this way, the columns in each row are 
 *   created. If the source indices map outside the size of the image, we 
 *   map zero to both indices and weights.
 * 
 * Updates:
 * June 9, 2015 : Created
 */

#include "mex.h"
#include<math.h>
#include "myMatrix.h"

/* This is the motion blur function. I use target-to-source mapping. For
 * every pixel in the target grid, for every homography and weight, I find
 * the correponding pixel in the source grid using the inverse homography.
 * I use bilinear interpolation with zero padding outside the source image
 * region. The target grid is offset by the user input, and the origin is
 * provided by the user.
 */
void warpMatrix(double *target, double *origH, int tgt_nrows, int tgt_ncols, 
        int src_nrows, int src_ncols, double src_orig_row, 
        double src_orig_col, int st_row, int st_col)
{
    // Number of pixels
    int tgt_npix, src_npix;
    
    // 2D point vars
    int src_row, src_col, f_src_row, f_src_col;
    double tgt_row, tgt_col, src_row_ni, src_col_ni;
    
    // 1D point vars
    int tgt_idx, src_idx;
    
    // Loop vars
    int i, j, k, ii, jj, m;
    
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
            
            tgt_idx = i*tgt_nrows + j; 
                
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
                
                /* For every pair of row,col in the four-point set */
                for (ii=0; ii<2; ii++) {
                    for (jj=0; jj<2; jj++) {
                        /* Get this row and this col number, and weight */
                        src_col = src_col_list[ii];
                        src_row = src_row_list[jj];
                        wt = wt_list[ii*2 + jj];
                        
                        /* Get 1D index of target and source points */
                        src_idx = (int)(src_col*src_nrows + src_row);
                        
                        if (src_row >=0 && src_row < src_nrows && src_col >=0 && src_col < src_ncols) {    
                            target[(ii*2 + jj)*tgt_npix + tgt_idx] = src_idx;
                            target[(4 + ii*2 + jj)*tgt_npix + tgt_idx] = wt;
                        }
                        else {
                            target[(ii*2 + jj)*tgt_npix + tgt_idx] = 0;
                            target[(4 + ii*2 + jj)*tgt_npix + tgt_idx] = 0;
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
    int nrows, ncols, tgt_nrows, tgt_ncols, ndim, nchan, i, j, nhomo;
    
    double *siz, *target, *H, *off, *tgt_siz, *orig;
    
    /* Create pointer to inputs */
    siz = mxGetPr(prhs[0]); // source image size
    H = mxGetPr(prhs[1]); // homography
    off = mxGetPr(prhs[2]); // offset of the target image
    tgt_siz = mxGetPr(prhs[3]); // size of the target image
    orig = mxGetPr(prhs[4]); // origin of the source grid
    
    nrows = (int) siz[0];
    ncols = (int) siz[1];
    
    /* Size of the target image */
    tgt_nrows = (int) tgt_siz[0];
    tgt_ncols = (int) tgt_siz[1];
   
    /* Create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(tgt_nrows*tgt_ncols,8,mxREAL);
    
    /* Get a pointer to the output matrix */
    target = mxGetPr(plhs[0]);
    
    /* Call the function */
    warpMatrix(target, H, tgt_nrows, tgt_ncols, nrows, ncols,
            orig[0], orig[1], (int)off[0], (int)off[1]);
}