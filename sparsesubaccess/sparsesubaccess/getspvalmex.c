/**************************************************************************
 * mex function getspvalmex.c
 * Get value of sparse matrix from subindices
 * Dichotomy search of indices
 * Calling:
 *   val = getspvalmex(S, i, j);
 *      S is the sparse matrix
 *      where i, j are subindex vectors a same length (double or int32)
 *      The output val is returned in column vector, same class as S
 *   [val nz] = getspvalmex(...) returns the number of zero in val
 *              (assuming S does not contain explicit zero according
 *               to Matlab convention)
 * Compile on 32-bit platform
 *  mex -O -v getspvalmex.c
 * On 64-bit platform
 *  mex -v -O -largeArrayDims getspvalmex.c
 * Author: Bruno Luong
 * Original: 11/April/2009
 **************************************************************************/

#include "mex.h"
#include "matrix.h"
#include "string.h"

/* Flag used only for debug purpose */
/* #define FINDSUBIDX_DEBUG */

/* Gateway routine getspvalmex */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    mwSize rows, columns;
    mwIndex nnz; /* , nzmax; */
    mwIndex i1, i9, imid;
    mwIndex *irs, *jcs;
    mwIndex jcsleft, jcsright;
    mwIndex k, ni, nj;
    double *idouble, *jdouble;
    int *iint32, *jint32;
    double *doublePr, *doublePi;
    mxLogical *logicalPr;
    double *si, *sr;
    mxLogical *sl;
    mwIndex i, j;
    mwSize dims[2];
    mxArray *val, *nzvalout;
    int cmplx;
    mxComplexity ComplexFlag;
    mxClassID IdxClass, DataClass;
    mwIndex nzval;
    
    if (nrhs<3)
        mexErrMsgTxt("Three input arguments required.");
    
    /* Check data type of input argument  */
    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("First input argument must be a sparse matrix.");
    
    if (!mxIsNumeric(prhs[1]))
        mexErrMsgTxt("Second input argument must be numeric.");
    
    if (!mxIsNumeric(prhs[2]))
        mexErrMsgTxt("Third input argument must be numeric.");
    
    /* Get the number of elements of the list of subindexes */
    ni = mxGetM(prhs[1])*mxGetN(prhs[1]);
    nj = mxGetM(prhs[2])*mxGetN(prhs[2]);
    
    /* They must be equal */
    if (ni!=nj)
        mexErrMsgTxt("i and j must have the same number of elements.");
    
    IdxClass = mxGetClassID(prhs[1]);
    /* They must be equal */
    if ((IdxClass = mxGetClassID(prhs[1])) != mxGetClassID(prhs[2]))
        mexErrMsgTxt("i and j must have the same class.");
        
    /* Get the pointer ob subindex arrays */
    irs = mxGetIr(prhs[0]);
    jcs = mxGetJc(prhs[0]);
    
    /* Pointers on data */
    DataClass =  mxGetClassID(prhs[0]);
    switch (DataClass) {
        case mxLOGICAL_CLASS:
#ifdef FINDSUBIDX_DEBUG
            mexPrintf("Logical\n");
#endif            
            sl  = mxGetLogicals(prhs[0]);
            sr = NULL;
            si = NULL;
            break;
        case mxDOUBLE_CLASS:
#ifdef FINDSUBIDX_DEBUG
            mexPrintf("Double\n");
#endif               
            sl = NULL;
            sr  = mxGetPr(prhs[0]);
            si  = mxGetPi(prhs[0]);
            break;
        default:
            mexErrMsgTxt("Only support data of class double and logical.");
    }
    
    /* Check if the input matrix is complex or real */
    cmplx = (si==NULL ? 0 : 1);
    
    /* NOTE: nnz is the actual number of nonzeros and is stored as the
     * last element of the jc array where the size of the jc array is the
     * number of columns + 1 */
    rows = mxGetM(prhs[0]);
    columns = mxGetN(prhs[0]);
    nnz = *(mxGetJc(prhs[0]) + columns);
    /*nzmax = mxGetNzmax(prhs[0]);*/
    
    if (cmplx)
        ComplexFlag = mxCOMPLEX;
    else
        ComplexFlag = mxREAL;
    
    /* Create a scalar containing a number of (implicit) zeros on second output */
    if (nlhs>=2) {
        dims[0] = 1; dims[1] = 1;
        plhs[1] = nzvalout = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        if (nzvalout==NULL)
        {
            mexErrMsgTxt("Out of memory.");
            return;
        }
        *mxGetPr(nzvalout) = 0;
    }
    
    /* Create output val as a long column vector */
    dims[0] = ni; dims[1] = 1;
    plhs[0] = val = mxCreateNumericArray(2, dims, DataClass, ComplexFlag);
    if (val==NULL) /* Cannot allocate memory */ {
        /* Try anyway to teturn an empty array */
        dims[0] = 0; dims[1] = 0;
        plhs[0] = mxCreateNumericArray(2, dims, DataClass, ComplexFlag);
        mexErrMsgTxt("Out of memory.");
        return;
    }
    
#ifdef FINDSUBIDX_DEBUG
    mexPrintf("nnz=%d\n", nnz);
#endif

    /* Fill with 0 */
    switch (DataClass) {
        case mxLOGICAL_CLASS:
            logicalPr = mxGetLogicals(val);
            memset((void*)logicalPr, 0, ni*sizeof(mxLogical));
            doublePr = (doublePi = NULL);
            break;
        case mxDOUBLE_CLASS:
            /* Get output data pointers */
            logicalPr = NULL;
            doublePr = mxGetPr(val);
            doublePi = mxGetPi(val);
            memset((void*)doublePr, 0, ni*sizeof(double));
            if (cmplx) memset((void*)doublePi, 0, ni*sizeof(double));
            break;
    }

    /* This boolean is to keep track if all value are different zero */
    nzval = 0; /* Yes by default */
    
    /* Here we go, working hard from here */
    switch (IdxClass) {
        
        case mxDOUBLE_CLASS:
            /* Get the pointer to subindex arrays */
            idouble = mxGetPr(prhs[1]);
            jdouble = mxGetPr(prhs[2]); 
            /* Loop on the list of subindexes */
            for (k=ni; k--;) /* Reverse of for (k=0; k<ni; k++) {...} */ {
                /* Get subindex value, substract 1 so they start from 0 */
                i = (mwIndex)(idouble[k]-1);
                j = (mwIndex)(jdouble[k]-1);
#ifdef FINDSUBIDX_DEBUG
                mexPrintf("look whereis (%d,%d)\n", i+1, j+1);
#endif
                /* Outside */
                if (i>=rows || j>=columns) {
#ifdef FINDSUBIDX_DEBUG
                 mexPrintf("\tfall outside\n");
#endif
                    nzval++;
                    continue;
                }
                
                /* jth column is empty */
                if ((jcsleft = jcs[j]) >= (jcsright = jcs[j+1])) {
#ifdef FINDSUBIDX_DEBUG
                    mexPrintf("\tcolumn is empty (all zero)\n");
#endif
                    nzval++;
                    continue;
                }
#ifdef FINDSUBIDX_DEBUG
                mexPrintf("\tjcsleft =  %d\n", jcsleft);
                mexPrintf("\tjcsright = %d\n", jcsright);
#endif
                
                /* Bracket */
                i1=0;
                i9=nnz-1;
                while (i9>i1+1) /* Dichotomy search */ {
                    imid = (i1+i9+1)/2; /* ceil of the mean */
                    if ((jcsleft>imid) /* Current column < j */ ||
                        ((jcsright>imid) && (i>irs[imid])))
                        /* Current column == j, current row > i */
                        i1=imid;
                    else
                        i9=imid;
                } /* of while loop */
                
#ifdef FINDSUBIDX_DEBUG
                mexPrintf("\ti1=%d\n", i1);
                mexPrintf("\ti9=%d\n", i9);
                mexPrintf("\tirs[i1]=%d\n", irs[i1]);
                mexPrintf("\tirs[i9]=%d\n", irs[i9]);
#endif
                /* i9 corresponds to (i,j) */
                if ((i9<jcsright) && (i==irs[i9])) {
#ifdef FINDSUBIDX_DEBUG
                    mexPrintf("\tmatch at i9\n");
#endif
                    if (sl)
                    {
                        logicalPr[k] = sl[i9];
                    }
                    else {
                        doublePr[k] = sr[i9];
                        if (cmplx) doublePi[k] = si[i9];
                    }
                }
                else if ((jcsleft<=i1) && (i==irs[i1])) /* i1 corresponds to (i,j) */ {
#ifdef FINDSUBIDX_DEBUG
                    mexPrintf("\tmatch at i1\n");
#endif
                    if (sl)
                    {
                        logicalPr[k] = sl[i1];
                    }
                    else
                    {
                        doublePr[k] = sr[i1];
                        if (cmplx) doublePi[k] = si[i1];
                    }
                }
                else nzval++; /* Not found */
                    
                
            } /* for-loop */
            break; /* mxDOUBLE_CLASS */
            
        case mxINT32_CLASS:
            iint32 = (int*)mxGetPr(prhs[1]);
            jint32 = (int*)mxGetPr(prhs[2]);
            /* Loop on the list of subindexes */
            for (k=ni; k--;) /* Reverse of for (k=0; k<ni; k++) {...} */ {
                /* Get subindex value, substract 1 so they start from 0 */
                i = (mwIndex)(iint32[k]-1);
                j = (mwIndex)(jint32[k]-1);
#ifdef FINDSUBIDX_DEBUG
                mexPrintf("look whereis (%d,%d)\n", i+1, j+1);
#endif
                /* Outside */
                if (i>=rows || j>=columns) {
#ifdef FINDSUBIDX_DEBUG
                 mexPrintf("\tfall outside\n");
#endif
                    nzval++;
                    continue;
                }
                
                /* jth column is empty */
                if ((jcsleft = jcs[j]) >= (jcsright = jcs[j+1])) {
#ifdef FINDSUBIDX_DEBUG
                    mexPrintf("\tcolumn is empty (all zero)\n");
#endif
                    nzval++;
                    continue;
                }
#ifdef FINDSUBIDX_DEBUG
                mexPrintf("\tjcsleft =  %d\n", jcsleft);
                mexPrintf("\tjcsright = %d\n", jcsright);
#endif
                
                /* Bracket */
                i1=0;
                i9=nnz-1;
                while (i9>i1+1) /* Dichotomy search */ {
                    imid = (i1+i9+1)/2; /* ceil of the mean */
                    if ((jcsleft>imid) /* Current column < j */ ||
                        ((jcsright>imid) && (i>irs[imid])))
                        /* Current column == j, current row > i */
                        i1=imid;
                    else
                        i9=imid;
                } /* of while loop */
                
#ifdef FINDSUBIDX_DEBUG
                mexPrintf("\ti1=%d\n", i1);
                mexPrintf("\ti9=%d\n", i9);
                mexPrintf("\tirs[i1]=%d\n", irs[i1]);
                mexPrintf("\tirs[i9]=%d\n", irs[i9]);
#endif
                /* i9 corresponds to (i,j) */
                if ((i9<jcsright) && (i==irs[i9])) {
#ifdef FINDSUBIDX_DEBUG
                    mexPrintf("\tmatch at i9\n");
#endif
                    if (sl)
                    {
                        logicalPr[k] = sl[i9];
                    }
                    else {
                        doublePr[k] = sr[i9];
                        if (cmplx) doublePi[k] = si[i9];
                    }
                }
                else if ((jcsleft<=i1) && (i==irs[i1])) /* i1 corresponds to (i,j) */ {
#ifdef FINDSUBIDX_DEBUG
                    mexPrintf("\tmatch at i1\n");
#endif
                    if (sl)
                    {
                        logicalPr[k] = sl[i1];
                    }
                    else {
                        doublePr[k] = sr[i1];
                        if (cmplx) doublePi[k] = si[i1];
                    }
                }
                else nzval++;
                
            } /* for-loop */            
            break; /* mxINT32_CLASS */
        default:
            mexErrMsgTxt("Only support subindex of class double and int32.");
    } /* switch */
    
    if (nlhs>=2) *mxGetPr(nzvalout) = (double)nzval;
  
    return;
    
} /* Gateway routine */

