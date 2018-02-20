/**************************************************************************
 * MATLAB mex function setspvalmex.c
 * Set value of sparse matrix from subindices
 * Calling:
 *   S = setspvalmex(S, i, j, val);
 *      S is the sparse matrix
 *      where i, j are subindex vectors a same length (double or int32)
 *      val is the same length tha i,j, and must have the same class as S
 *      Return S with new assigned values: S(i(k),j(k)) = val(k) for all k
 *   Use last value to assign in case of conflict of position
 *   Dimension of S is automatically enlarged to fit all the indices (i,j)
 *
 * Compile on 32-bit platform
 *  mex -O -v setspvalmex.c
 * On 64-bit platform
 *  mex -v -O -largeArrayDims setspvalmex.c
 *
 * Algorithm: couple of subindices (i,j) are converted to unique 64-bit
 *            integers, quicksort/insertion-sort engine work on int64
 *            with pivot strategy using median of three fixed points
 *            The sort result are then merge with input sparse matrix
 *            Explicit zeros are then removed
 *
 * Author Bruno Luong <brunoluong@yahoo.com>
 * Last update:
 *  18-Aug-2009 tolerate setting UINT8/INT8 to Logical sparse
 *  16-Dec-2009 replace "|" operator by "+" at three places as workaround
 *              of LCC bug
 *  10-Jan-2010 use introsort (safer than quicksort)
 *              correct bug in Median of three pivoting
 *************************************************************************/

#include "mex.h"
#include "matrix.h"
#include "string.h"
#include "math.h"

/* Define correct type depending on platform */
#if defined(_MSC_VER) || defined(__BORLANDC__)
typedef unsigned __int32 uint32;
typedef signed __int32 int32;
typedef unsigned __int64 ulong64;
typedef signed __int64 long64;
#define INTMIN64 0x8000000000000000
#else
typedef unsigned int uint32;
typedef signed int int32;
typedef unsigned long long ulong64;
typedef signed long long long64;
#define INTMIN64 0x8000000000000000ULL
#endif

/* Different options for Pivot selection */
/* Different options for Pivot selection */
#define MIDPOINT      0
#define MEDIAN3       1
#define MEDIANMEDIANS 2

/* Pivot Strategy, use one of the above */
#define PIVOT MEDIAN3

#define INSERTSORT  0
#define SHELLSORT   1
#define QUICKSORT   2
#define HEAPSORT    3
#define INTROSORT   4

/* Method for sorting small arrays */
/* INSERTSORT | SHELLSORT */
#define BASESORT INSERTSORT 

/* Method dor sorting large arrays */
/* QUICKSORT | HEAPSORT | INTROSORT */
#define MAINSORT INTROSORT

/* Increment for shellshort */
#define SizeIncs 3
mwIndex incs[SizeIncs] = { 7, 3, 1 };

/* A limit below the recusion stops, we which we switch to InsertionSort */
/* This must be smaller than */
#define SWITCH_ALGO_THRESHOLD 12

/* Type of elements to be sorted */
#define ARRAYTYPE long64

/* Global variables, used to avoid stacking them during recusive call since
   they do not change */
mwIndex *pos;
ARRAYTYPE *list;

/*************************************************************************/
/* Remove explicite zeros from sparse matrix, conform to Matlab
 * requirement. Modify data under the pointer directly (in-place).
 * nnz might changes but nzmax stays. 
 * S can be logical, real-double or real-complex sparse matrix. */
void SpRemoveZ(mxArray *S)
{
    mxLogical *sl, *src_sl, *dst_sl;
    double *sr, *src_sr, *dst_sr;
    double *si, *src_si, *dst_si;
    mwIndex *src_irs, *dst_irs;
    mwIndex nz, knextj, k, nnz;
    mwIndex *irs, *jcs;
    mwSize j, columns;
    mxClassID DataClass;
    
    /* Check if the matrix is sparse */
    if (!mxIsSparse(S))
        mexErrMsgTxt("First input argument must be a sparse matrix.");
    
    /* Get class, number of colmns, and number of non-zero elements */
    columns = mxGetN(S);
    DataClass = mxGetClassID(S);
    nnz = *(mxGetJc(S) + columns);

    /* Check the type and initialize data pointers */
    switch (DataClass) {
        case mxLOGICAL_CLASS:
            sl = mxGetLogicals(S);
            sr = NULL;
            si = NULL;            
            break;
        case mxDOUBLE_CLASS:
            sl = NULL;
            sr = mxGetPr(S);
            si = mxGetPi(S);
            break;
    } /* Switch */

    /* Get the index pointers */
    irs = mxGetIr(S);
    jcs = mxGetJc(S);
    
    /* Initialize source and destination pointers (irs) */
    src_irs = irs; dst_irs = irs;
    /* Look for first non empty column */
    j=0;
    while (j+1<columns && jcs[j+1]==0) j++;
    knextj = jcs[j+1];
    /* variable to keep track the number of explicite zeros */
    nz = 0;
    if (sl) { /* Logicals */
        
        src_sl = sl; dst_sl = sl;
        
        /* Loop to remove zeros */
        for (k=0; k<nnz; k++) {
            /* Update the column index j */
            while (k>=knextj) {
                jcs[++j] -= nz;
                knextj = jcs[j+1];                     
             }
             if (*(src_sl)!=0) /* Non-zero element */
             {
                *(dst_sl++) = *(src_sl++);
                *(dst_irs++) = *(src_irs++);
             }
             else /* Ignore explicit zero */ {                      
                src_sl++;
                src_irs++;
                nz++;
             }
        } /* for-loop */
    }
    else if (!si) { /* Non-complex double only */
        
        src_sr = sr; dst_sr = sr;
        
        /* Loop to remove zeros */
        for (k=0; k<nnz; k++) {
            /* Update the column index j */
            while (k>=knextj) {
                jcs[++j] -= nz;
                knextj = jcs[j+1];                       
            }
            if (*(src_sr)!=0.0) /* Non-zero element */ {
                *(dst_sr++) = *(src_sr++);
                *(dst_irs++) = *(src_irs++);
            }
            else /* Ignore explicit zero */ {                       
                src_sr++;
                src_irs++;
                nz++;
            }
        } /* for-loop */
    }
    else { /* Complex array */
        
        src_sr = sr; dst_sr = sr;
        src_si = si; dst_si = si;
        
        /* Loop to remove zeros */
        for (k=0; k<nnz; k++) {
            /* Update the column index j */
            while (k>=knextj) {
                jcs[++j] -= nz;
                knextj = jcs[j+1];                      
            }
            if (*(src_sr)!=0.0 || *(src_si)!=0.0) /* Non-zero element */ {
                *(dst_sr++) = *(src_sr++);
                *(dst_si++) = *(src_si++);
                *(dst_irs++) = *(src_irs++);
            }
            else /* Ignore explicit zero */ {                        
                src_sr++;
                src_si++;
                src_irs++;
                nz++;
            }
        } /* for-loop */
    } /* All cases (Logical, real, imag) are treated */
                
    /* Adjust the number of non-zero for the rest of jcs */
    nnz = nnz - nz;
    for (++j; j<=columns; j++) jcs[j] = nnz;
        
    return;
} /* SpRemoveZ */

/*************************************************************************/
/* Insertion sort, this is stable */
void InsertionSort(mwIndex left, mwIndex right)
{
    mwIndex pival;
    mwIndex *pi, *pj, *pleft, *pright;
    ARRAYTYPE Value;

    pleft = pos+left;
    pright = pos+right;
    for (pi=pleft; pi<pright; pi++) {
        Value = list[pival=(*(pi+1))];
        pj = pi;
        
        while ((pj>=pleft) && (list[*pj]>Value)) { /* Comparison */
            *(pj+1) = *pj; /* Permute */
            pj--;
        }
        *(pj+1) = pival; /* Insert */
    } /* for */
    
    return;
} /* InsertionSort */

/*************************************************************************/

/* Shell sort */            
void ShellSort(mwIndex left, mwIndex right)
{
    mwIndex pival, pjval, h, k;
    mwIndex *pi, *pj, *ph;
    mwIndex *pleft, *pright;
    ARRAYTYPE Value;
    
    if (right > left) {
        
        pleft = pos+left;
        pright = pos+right;
        
        /* Loop on gap size */
        for ( k = 0; k < SizeIncs; k++) {
            h = incs[k];
            ph = pos + h;
            /* Insertion loop */
            for (pi=pleft+h; pi<=pright; pi++) 
            {
                Value = list[pival=*pi];
                pj = pi;
                while ((pj>=ph) && (list[pjval=*(pj-h)]>Value)) { /* Comparison */
                    *pj = pjval; /* Permute */
                    pj -= h;
                }
                *pj = pival; /* Insert */
            }
        }
    }
    
    return;
} /* ShellSort */

/*************************************************************************/
/*Find the index of the Median of the elements
of array that occur at every "shift" positions.*/
mwIndex findMedianIndex(mwIndex left, mwIndex right, mwIndex shift)
{
    mwIndex tmp, groups, n;
    ARRAYTYPE minValue;
    mwIndex *pi, *pj, *pk, *pright, *pminIndex;
    
    groups = (right-left)/shift + 1;
    pk = pos + (n = left + (groups/2)*shift);
    pright = pos + right;
    for (pi=pos+left; pi<=pk; pi+= shift)
    {
        pminIndex = pi;
        minValue = list[*pminIndex];
        
        for (pj=pi; pj<=pright; pj+=shift)
            if (list[*pj]<minValue) /* Comparison */
                minValue = list[*(pminIndex=pj)];
        /* Swap pos[i] with pos[minIndex] */
        tmp = *pi;
        *pi = *pminIndex;
        *pminIndex = tmp;
    }
    
    return n;
    
} /* findMedianIndex */

/*Computes the median of each group of 5 elements and stores
  it as the first element of the group (left). Recursively does this
  till there is only one group and hence only one Median */
mwIndex findMedianOfMedians(mwIndex left, mwIndex right)
{
    mwIndex i, shift, step, tmp;
    mwIndex endIndex, medianIndex;
           
    if (left==right) return left;
   
    shift = 1;
    while (shift <= (right-left))
    {
        step=shift*5;
        for (i=left; i<=right; i+=step)
        {
            if ((endIndex=i+step-1)>=right)
                endIndex=right;
            medianIndex = findMedianIndex(i, endIndex, shift);
            /* Swap pos[i] with pos[medianIndex] */
            tmp = pos[i];
            pos[i] = pos[medianIndex];
            pos[medianIndex] = tmp;
        }
        shift = step;
    }
    return left;
} /* findMedianOfMedians */

/*************************************************************************/
/*Computes the median of three points (left,right,and mid) */
mwIndex findMedianThree(mwIndex left, mwIndex right)
{
    ARRAYTYPE vleft, vright, vmid;
    mwIndex mid;
    
    if (left==right) return left;
    
    vleft = list[pos[left]];
    vright = list[pos[right]];
    vmid = list[pos[mid = (left+right+1)/2]];
    
    if (vleft<vright)
    {
        if (vmid>vright)
            return right;
        else if (vmid<vleft)
            return left;
        else
            return mid;
        
    } else { /* (vleft>=vright) */
        
        if (vmid>vleft)
            return left;
        else if (vmid<vright)
            return right;
        else
            return mid;
        
    }       
} /* findMedianThree */

/*************************************************************************/
/* Partitioning the list around pivot pivotValue := l[pivotIndex];
 * After runing, at exit we obtain: 
   l[left]...l[index-1] < pivotValue <= l[index] ... l[right]
   where l[i] := list[pos[i]] for all i */
mwIndex PartitionV1(mwIndex left, mwIndex right, mwIndex pivotIndex) {
    
    ARRAYTYPE pivotValue;
    mwIndex *pindex, *pi, *pright;
    mwIndex tmp;
    
    pright=pos+right;
    pindex=pos+pivotIndex;
    pivotValue = list[tmp = *pindex];
    /* Swap pos[pivotIndex] with pos[right] */
    *pindex = *pright;
    *pright = tmp;
    
    pindex=pos+left;
    for (pi=pindex; pi<pright; pi++)
        /* Comparison with pivotValue */
        if (list[*pi] < pivotValue) {
             /* if smaller; Swap pos[index] with pos[i] */
            tmp = *pindex;
            *pindex = *pi;
            *pi = tmp;           
            pindex++;
        }

     /* Swap pos[index] with pos[right] */
    tmp = *pindex;
    *pindex = *pright;
    *pright = tmp;  
    
    return (mwIndex)(pindex-pos); /* Pointer arithmetic */
} /* PartitionV1 */

mwIndex PartitionV2(mwIndex left, mwIndex right, mwIndex pivotIndex) {
    
    ARRAYTYPE pivotValue;
    mwIndex *pindex, *pright, *pleft;
    mwIndex *pfirst, *plast;
    mwIndex tmp;
    
    plast=pos+right;
    pindex=pos+pivotIndex;
    pivotValue = list[tmp = *pindex];
    /* Swap pos[pivotIndex] with pos[right] */
    *pindex = *plast;
    *plast = tmp;
    
    pfirst = pleft = pos+left;
    pright = plast-1;
    for (;;) {
        while (list[*pleft]<pivotValue)
            pleft++;
        while ((pright>pfirst) && (list[*pright]>pivotValue))
            pright--;
        if (pleft<pright) {
            tmp = *pleft;
            *pleft = *pright;
            *pright = tmp;
            pleft++, pright--;
        }
        else {
            tmp = *plast;
            *plast = *pleft;
            *pleft = tmp;
            
            return (mwIndex)(pleft-pos); /* Pointer arithmetic */
        }
    }
} /* PartitionV2 */

int CheckPartition(mwIndex left, mwIndex right, mwIndex pivotIndex)
{
    mwIndex i;
    ARRAYTYPE pivotValue;
    
    pivotValue = list[pos[pivotIndex]];
    for (i=left;i<=pivotIndex-1;i++)
        if (list[pos[i]]>pivotValue)
            return 0;
    for (i=right;i>=pivotIndex+1;i--)
        if (list[pos[i]]<pivotValue)
            return 0;
    return 1;
}

/*************************************************************************/
void DownHeap(mwIndex i, mwIndex n, mwIndex lo) {
    ARRAYTYPE d;
    mwIndex child, posd, nhalf;
    mwIndex *plo, *pchild; 
    
    plo = pos + lo;
    posd = *(plo+i-1);
    d = list[posd];
    nhalf = n/2;
    while (i<=nhalf) {
        child = 2*i;
        pchild = plo + child;
        if ((child < n) && 
            (list[*(pchild-1)] < list[*pchild]))
        {
            child++;
            pchild++;
        }
        if (d >= list[*(pchild-1)]) break;
        *(plo+i-1) = *(pchild-1);
        i = child;
    }
    *(plo+i-1) = posd;
    
    return;
} /* DownHeap */

/*http://users.encs.concordia.ca/~chvatal/notes/hsort.html*/
void HeapSort(mwIndex lo, mwIndex hi) { /* (lo,hi) included */
    mwIndex n, i, tmp;
    mwIndex *plo;
    
    hi++;
    n = hi-lo;
            
    for (i=n/2; i>=1; i--)
        DownHeap(i, n, lo);

    plo = pos + lo;
    for (i=n; i>1; i--) {
        
        tmp = *(plo);
        *(plo) = *(plo+i-1);
        *(plo+i-1) = tmp;
        
        DownHeap(1, i-1, lo);
    }
    return;
    
} /* HeapSort */

/*************************************************************************/

void introsort_loop(mwIndex lo, mwIndex hi, int depth_limit) 
 {
    mwIndex pivotIndex;
    while (hi+1 >= lo + SWITCH_ALGO_THRESHOLD) {
        if (depth_limit == 0) {
            HeapSort(lo, hi);
            return;
        }
        depth_limit--;
        
#if (PIVOT==MEDIANMEDIANS)
        pivotIndex = findMedianOfMedians(lo, hi);
#elif (PIVOT==MEDIAN3)
        pivotIndex = findMedianThree(lo, hi);
#else /* MIDPOINT */
        pivotIndex = (lo+hi+1)/2;
#endif       
        pivotIndex = PartitionV2(lo, hi, pivotIndex);        
        introsort_loop(pivotIndex+1, hi, depth_limit);
        hi = pivotIndex-1;
    }
    return;
} /* introsort_loop */

int floor_lg(mwIndex a) {
    return (int)(floor(log((double)a)/log(2.0)));
} /* floor_lg */


void IntroSort(mwIndex left, mwIndex right)
{
    int depthlimit;
    depthlimit=2*floor_lg(right-left+1);
    
    introsort_loop(left, right, depthlimit);
#if (BASESORT==INSERTSORT)
    InsertionSort(left, right);
#else /* #if (BASESORT==SHELLSORT)*/
    ShellSort(left, right);
#endif  
    
    return;
} /* IntroSort */

/*************************************************************************/
/* Recursive engine quicksort */
void quicksort(mwIndex left, mwIndex right) {
    
    mwIndex pivotIndex;

    /* Switch to Insertion sort for small size array */
    if (right+1 < left + SWITCH_ALGO_THRESHOLD)
    {
#if (BASESORT==INSERTSORT)
        InsertionSort(left, right);
#else /* #if (BASESORT==SHELLSORT)*/
        ShellSort(left, right);
#endif        
    }
    else
    {        
#if (PIVOT==MEDIANMEDIANS)
        pivotIndex = findMedianOfMedians(left, right);
#elif (PIVOT==MEDIAN3)
        pivotIndex = findMedianThree(left, right);
#else /* MIDPOINT */
        pivotIndex = (left+right+1)/2;
#endif
        
        pivotIndex = PartitionV2(left, right, pivotIndex);
        quicksort(left, pivotIndex-1);
        quicksort(pivotIndex+1, right);
    }
    
    return;
} /* quicksort */

/* Alloacate and Initialize the table of position before sorting */
ARRAYTYPE* AllocatePos(mwIndex l, mxClassID IdxClass, void *i, void *j)
{
    mwIndex n;
    int32 *i32, *j32;
    double *idbl, *jdlb;
    ARRAYTYPE *l64;
    
    /* clean programming */    
    pos=NULL;
    list=NULL;
   
    if (l>0) {
        pos = mxMalloc(sizeof(mwIndex)*l);
        if (pos==NULL)
            mexErrMsgTxt("Out of memory (pos).");
        /* Initialize the array of position (zero-based index) */
        for (n=0; n<l; n++) pos[n]=n;
        
        list = mxMalloc(sizeof(ARRAYTYPE)*l);
        if (list==NULL)
            mexErrMsgTxt("Out of memory (list).");
    }
    
    /* Point at the begining of the list */
    l64 = list;
    
    /* Convert couple of indexes (zero-based) into a single int64 number
     * to facilite sorting, column index J has greater priority */
    switch (IdxClass) {        
        case mxDOUBLE_CLASS:
            idbl = (double*)i; jdlb = (double*)j;
            for (n=0; n<l; n++)
                *(l64++) = ((ARRAYTYPE)(*(jdlb++)-1) << 32) + /* LCC BUG */
                            (ARRAYTYPE)(*(idbl++)-1);
            break;
        case mxINT32_CLASS:
            i32 = (int32*)i; j32 = (int32*)j;
            /* Could do faster here with interleave copy
               and indian characteristic of the platform */
            for (n=0; n<l; n++) 
                *(l64++) = ((ARRAYTYPE)(*(j32++)-1) << 32) + /* LCC BUG */
                            (ARRAYTYPE)(*(i32++)-1);
            break;
        default:
            mexErrMsgTxt("Only support subindex of class double and int32.");
    }
    
    return list;
} /* AllocatePos */

/* Free the position */
void FreePos(void)
{
    /* Free the array of position */
    if (pos) mxFree(pos);
    pos = NULL; /* clean programming */
    
    /* Free the array of values int64 */
    if (list) mxFree(list);
    list = NULL; /* clean programming */    
    
    return;
} /* FreePos */

         
/* Return the platform dependent maximum size of for mwIndex */
mwIndex MaxSize(void) {
    mwIndex mxsz;
    /*mxsz = 0x7fffffff;*/
    switch (sizeof(mxsz)) {
        case 4:
            mxsz = (mwIndex)0x7fffffff; /* 2^31-1 */
            break;
        case 8:
            mxsz = (mwIndex)0xffffffffffff; /* 2^48-1 */
            break;
    }
    return mxsz;
} /* MaxSize */

/* Find index in sparse matrix S, allocate a buffer (with extra space
 * in the tail for later mergin purpose), return the result in array of 
 * type int64 for combined sub-indexes */
ARRAYTYPE* spfind(const mxArray *S)
{
    mwIndex n, j, nnz;
    mwIndex *irs, *jcs;
    mwSize columns;
    ARRAYTYPE *buffer, *l64;
    
    /* Get class, number of colmns, and number of non-zero elements */
    columns = mxGetN(S);
    nnz = *(mxGetJc(S) + columns);

    /* Get the index pointers */
    irs = mxGetIr(S);
    jcs = mxGetJc(S);
    
    /* Allocate arrays of (nnz + extralength) */
    if (nnz==0)
        return NULL;        
    
    buffer = mxMalloc(sizeof(ARRAYTYPE)*nnz);
    if (buffer==NULL)
        mexErrMsgTxt("Out of memory (subidx buffer).");
    
    /* Look for first non empty column */
    j=0;
    while (j+1<columns && jcs[j+1]==0) j++;
    /* Move to the tail, leave the head free */ 
    l64 = buffer;
    for (n=0; n<nnz; n++) {
        /* Update the column */
        while (n>=jcs[j+1]) j++;
        /* Combined index */
        *(l64++) = ((ARRAYTYPE)j << 32) + (ARRAYTYPE)(*(irs++)); /* LCC BUG */
    } /* for-loop */
   
    return buffer;
} /* spfind */


/* Mege engine for complex double arrays */
void mergecmplx(ARRAYTYPE* list, mwIndex* pos, ARRAYTYPE* buffer, 
                mwIndex i, mwIndex upper, mwIndex *irs, mwIndex *jcs,
                double *headdbl, double *taildbl, double*sr,
                double *headimg, double *tailimg, double *si)
{    
    mwIndex *p, pprev;
    ARRAYTYPE idx, iprev;
    ARRAYTYPE* source;
    mwIndex nnz;
    mwIndex j, k; 
 
    /* Fill imaginary by zero */
    memset((void*)si, 0, sizeof(double)*upper);
    
    /* p points on the pos */
    p = pos;    
    /* copy back next-greatest element at each time */
    iprev = INTMIN64; /* Smallest int64 */
    pprev = -1; /* Something smaller than 0 */
    /* k is the number of elements */
    k = 0;
    nnz = 0;
    source = buffer; 
    /* Main loop of merge */
    while (k<i && i<upper)
    {
        if (list[*p] <= (*source)) /* Copy new element in case of draw */
        {
            /* Add a new element into a matrix */
            if ((idx = list[*p]) > iprev) {            
                iprev = idx; /* Keep track old indices */
                *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
                j = (mwIndex)(idx >> 32); /* Last 32-bits */
                jcs[j+1] = ++nnz;
                *(sr++) = headdbl[pprev=(*p)]; /* Copy new data */
                if (headimg)
                    *(si++) = headimg[*p];
                else
                    si++;
            } else if ((*p)>pprev) /* Select the last element if draw */ {
                *(sr-1) = headdbl[pprev=(*p)];
                if (headimg)
                    *(si-1) = headimg[*p];              
            }
            if (idx == (*source)) { /* Ignore data on source if conflict */
                source++;
                taildbl++;
                if (tailimg) tailimg++;
                i++;
                k++;
            }
            p++;
        }
        else /* existing data */
        {            
            idx = *(source++);
            *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
            j = (mwIndex)(idx >> 32);  /* Last 32-bits */
            jcs[j+1] = ++nnz;            
            i++;
            *(sr++) = *(taildbl++); /* Copy existing data */ 
            if (tailimg)
                *(si++) = *(tailimg++);
            else si++;
        }
        k++;
    }
   
    /* copy back remaining elements of new elements list (if any) */
    while (k<i) /* source is equal to upper if we enter here */
    {     
        /* Add a new element into a matrix */
        if ((idx = list[*p]) > iprev) {
            iprev = idx;
            *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
            j = (mwIndex)(idx >> 32); /* Last 32-bits */
            jcs[j+1] = ++nnz;
            *(sr++) = headdbl[pprev=(*p)];
            if (headimg)
                *(si++) = headimg[*p];
            else
                si++;
        } else if (*p>pprev) /* Select the last element if draw */ {
            *(sr-1) = headdbl[pprev=(*p)];
            if (headimg)
                *(si-1) = headimg[*p];
        }
        k++;
        p++;
    }

    /* Copy back remaining elements of existing elements (if any) */
    while (i<upper)
    {        
        idx = *(source++);
        *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
        j = (mwIndex)(idx >> 32); /* Last 32-bits */
        jcs[j+1] = ++nnz;
        *(sr++) = *(taildbl++);
        if (tailimg)
            *(si++) = *(tailimg++);
        i++;       
    }
    
    return;
} /* mergecmplx */

/* Mege engine for double arrays */
void mergedbl(ARRAYTYPE* list, mwIndex* pos, ARRAYTYPE* buffer, 
              mwIndex i, mwIndex upper, mwIndex *irs, mwIndex *jcs,
              double *headdbl, double *taildbl, double *sr)
{    
    mwIndex *p, pprev;
    ARRAYTYPE idx, iprev;
    ARRAYTYPE* source;
    mwIndex nnz;    
    mwIndex j, k; 
    
    /* p points on the pos */
    p = pos;    
    /* copy back next-greatest element at each time */
    iprev = INTMIN64; /* Smallest int64 */
    pprev = -1; /* Something smaller than 0 */
    /* k is the number of elements */
    k = 0;
    nnz = 0;
    source = buffer; 
    /* Main loop of merge */
    while (k<i && i<upper)
    {
        if (list[*p] <= (*source)) /* Copy new element in case of draw */
        {
            /* Add a new element into a matrix */
            if ((idx = list[*p]) > iprev) {            
                iprev = idx; /* Keep track old indices */
                *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
                j = (mwIndex)(idx >> 32); /* Last 32-bits */
                jcs[j+1] = ++nnz;
                *(sr++) = headdbl[pprev=(*p)]; /* Copy new data */                              
            } else if ((*p)>pprev) /* Select the last element if draw */
                *(sr-1) = headdbl[pprev=(*p)];              
            if (idx == (*source)) { /* Ignore data on source if conflict */
                source++;
                taildbl++;
                i++;
                k++;
            }
            p++;
        }
        else /* existing data */
        {            
            idx = *(source++);
            *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
            j = (mwIndex)(idx >> 32);  /* Last 32-bits */
            jcs[j+1] = ++nnz;            
            i++;
            *(sr++) = *(taildbl++); /* Copy existing data */            
        }
        k++;
    }
   
    /* copy back remaining elements of new elements list (if any) */
    while (k<i) /* source is equal to upper if we enter here */
    {     
        /* Add a new element into a matrix */
        if ((idx = list[*p]) > iprev) {
            iprev = idx;
            *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
            j = (mwIndex)(idx >> 32); /* Last 32-bits */
            jcs[j+1] = ++nnz;
            *(sr++) = headdbl[pprev=(*p)];            
        } else if (*p>pprev) /* Select the last element if draw */
            *(sr-1) = headdbl[pprev=(*p)];
        k++;
        p++;
    }
 
    /* Copy back remaining elements of existing elements (if any) */
    while (i<upper)
    {        
        idx = *(source++);
        *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
        j = (mwIndex)(idx >> 32); /* Last 32-bits */
        jcs[j+1] = ++nnz;
        *(sr++) = *(taildbl++);
        i++;       
    }
    
    return;
} /* mergedbl */

/* Mege engine for logical arrays */
void mergelog(ARRAYTYPE* list, mwIndex* pos, ARRAYTYPE* buffer, 
              mwIndex i, mwIndex upper, mwIndex *irs, mwIndex *jcs,
              mxLogical *headlog, mxLogical *taillog, mxLogical *sl)
{  
    mwIndex *p, pprev;
    ARRAYTYPE idx, iprev;
    ARRAYTYPE* source;
    mwIndex nnz;    
    mwIndex j, k; 
    
    /* p points on the pos */
    p = pos;    
    /* copy back next-greatest element at each time */
    iprev = INTMIN64; /* Smallest int64 */
    pprev = -1; /* Something smaller than 0 */
    /* k is the number of elements */
    k = 0;
    nnz = 0;
    source = buffer; 
    /* Main loop of merge */
    while (k<i && i<upper)
    {
        if (list[*p] <= (*source)) /* Copy new element in case of draw */
        {
            /* Add a new element into a matrix */
            if ((idx = list[*p]) > iprev) {            
                iprev = idx; /* Keep track old indices */
                *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
                j = (mwIndex)(idx >> 32); /* Last 32-bits */
                jcs[j+1] = ++nnz;
                *(sl++) = headlog[pprev=(*p)]; /* Copy new data */                              
            } else if ((*p)>pprev) /* Select the last element if draw */
                *(sl-1) = headlog[pprev=(*p)];              
            if (idx == (*source)) { /* Ignore data on source if conflict */
                source++;
                taillog++;
                i++;
                k++;
            }
            p++;
        }
        else /* existing data */
        {            
            idx = *(source++);
            *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
            j = (mwIndex)(idx >> 32);  /* Last 32-bits */
            jcs[j+1] = ++nnz;            
            i++;
            *(sl++) = *(taillog++); /* Copy existing data */            
        }
        k++;
    }
   
    /* copy back remaining elements of new elements list (if any) */
    while (k<i) /* source is equal to upper if we enter here */
    {     
        /* Add a new element into a matrix */
        if ((idx = list[*p]) > iprev) {
            iprev = idx;
            *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
            j = (mwIndex)(idx >> 32); /* Last 32-bits */
            jcs[j+1] = ++nnz;
            *(sl++) = headlog[pprev=(*p)];            
        } else if (*p>pprev) /* Select the last element if draw */
            *(sl-1) = headlog[pprev=(*p)];
        k++;
        p++;
    }
 
    /* Copy back remaining elements of existing elements (if any) */
    while (i<upper)
    {        
        idx = *(source++);
        *(irs++) = (mwIndex)(idx & 0xffffffff); /* First 32-bits */
        j = (mwIndex)(idx >> 32); /* Last 32-bits */
        jcs[j+1] = ++nnz;
        *(sl++) = *(taillog++);
        i++;       
    }
    
    return;
} /* mergelog */

void InitJcs(mwIndex* jcs, mwSize n)
{
    /* Fill jcs by zero */
    memset((void*)jcs, 0, sizeof(mwIndex)*(n+1));
    return;
} /* InitJcs */

void UpdateJcs(mwIndex* jcs, mwSize n)
{
    mwSize j;
    /* Fill jcs for empty columns */
    for (j=1; j<=n; j++)
        if (jcs[j]==0) jcs[j]=jcs[j-1];
    return;
} /* UpdateJcs */

/* Merge the list of existing index (tail of buffer) with a new list 
 of index. */
mxArray* MergeIdx(ARRAYTYPE* list, mwIndex* pos, mwIndex nhead,
                  ARRAYTYPE* buffer, mwIndex ntail,                  
                  mwSize m, mwSize n, mwIndex nzmax,
                  mxClassID DataClass,
                  void* Headr, void* Tailr,
                  double* Headi, double* Taili)
{
    mwIndex i, upper;
    mxArray *S;
    mwIndex *irs, *jcs;
    double *sr, *si;
    mxLogical *sl;
    mwIndex k; 
    mwSize row, maxrow, maxcol;
    int cmplx;
    double *headdbl, *taildbl;
    mxLogical *headlog, *taillog;
    
    /* Upper is to the right, i point to mid */
    upper = (i=nhead) + ntail;   
    
    /* Change here for complex matrix */
    cmplx = (Headi!=NULL || Taili!=NULL);
    
    /* Create a new sparse matrix */
    if (nzmax < upper) nzmax = upper;
    
    /* Minimum size for rows and columns of the new elements */
    if (nhead>0) {
        maxrow = 0;
        /* Find max of row */
        for (k=0; k<nhead; k++) {
            row = (mwIndex)(list[k] & 0xffffffff); /* First 32-bits */
            if (row>maxrow) maxrow=row;
        }
        /* Last element has a largest column */
        maxcol = (mwIndex)(list[pos[nhead-1]] >> 32);  /* Last 32-bits */ 
        /* One-based indexing */
        maxrow++; /* One-based indexing */
        maxcol++;            
        /* Enlarge the dimensions if needed */
        if (m<maxrow) m=maxrow;
        if (n<maxcol) n=maxcol;        
    }

    /* Create a new sparse output */
    switch (DataClass) {
        case mxLOGICAL_CLASS:
            S = mxCreateSparseLogicalMatrix(m, n, nzmax);
            if (S==NULL)
                mexErrMsgTxt("Out of memory (create new sparse matrix).");
            sl = mxGetLogicals(S);
            sr = NULL;
            si = NULL;
            
            cmplx = 0;
            
            headlog = (mxLogical*)Headr;
            taillog = (mxLogical*)Tailr;
            
            break;
        case mxDOUBLE_CLASS:
            /* Creat sparse here */
            S = mxCreateSparse(m, n, nzmax, cmplx);
            if (S==NULL)
                mexErrMsgTxt("Out of memory (create new sparse matrix).");
            sl = NULL;
            sr  = mxGetPr(S);
            si  = mxGetPi(S);
            
            headdbl = (double*)Headr;
            taildbl = (double*)Tailr;

            break;
    }
    
    irs = mxGetIr(S);
    jcs = mxGetJc(S);
    
    /* Fill jcs by zero */
    InitJcs(jcs, n);
   
    /* Create a new sparse output */
    switch (DataClass) {
        case mxLOGICAL_CLASS:
            /* merge logical */
            mergelog(list, pos, buffer, i, upper, irs, jcs, headlog, taillog, sl);
            break;
        case mxDOUBLE_CLASS:
            if (cmplx) /* merge complex */
                mergecmplx(list, pos, buffer, i, upper, irs, jcs, 
                           headdbl, taildbl, sr,
                           Headi, Taili, si);
            else /* merge double */
                mergedbl(list, pos, buffer, i, upper, irs, jcs, 
                         headdbl, taildbl, sr);
            break;
    }
    
    /* Fill jcs for empty columns */
    UpdateJcs(jcs, n);

    /* Return the result */
    return S;
    
} /* MergeIdx */


/* MainSort64, sorting interface */
void MainSort64(ARRAYTYPE *list, mwIndex l) {

    /* Work for non-empty array */
    if (l>1 && list!=NULL && pos!=NULL) {
        /* Call the recursive engine */
#if (MAINSORT==QUICKSORT)
        quicksort(0, l-1);
#elif (MAINSORT==INTROSORT)
        IntroSort(0, l-1);
#elif (MAINSORT==HEAPSORT)
        HeapSort(0, l-1);
#else
#error "Non valid MAINSORT"
#endif
    } /* if (l>0) */
   
    return;

} /* MainSort64 */

/* Free a pointer */
void* MyFree(void* Ptr)
{
    if (Ptr)
        mxFree(Ptr);
    return NULL;
} /* Free */


/* Gateway routine setspvalmex */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    mwIndex ni, nj, nv;
    mwIndex l;
    const mxArray *val, *Sin;
    mxArray *Sout;
    mxClassID DataClass, ValClass, IdxClass;
    ARRAYTYPE *buffer;
    mwIndex rows, columns, nnz, nzmax;
    
    if (nrhs!=4)
        mexErrMsgTxt("Four input arguments required.");
    
    /* Check data type of input argument  */

    Sin = prhs[0];
    if (!mxIsSparse(Sin))
        mexErrMsgTxt("First input argument must be a sparse matrix.");
    
    if (!mxIsNumeric(prhs[1]))
        mexErrMsgTxt("Second input i argument must be numeric.");
    
    if (!mxIsNumeric(prhs[2]))
        mexErrMsgTxt("Third input j argument must be numeric.");
    
    val = prhs[3];
    if (!mxIsNumeric(val) && !mxIsLogical(val))
        mexErrMsgTxt("Fourth input val argument must be numeric or logical.");    
    
    /* Get the number of elements of the list of subindexes */
    l = ni = (mwIndex)mxGetM(prhs[1])*(mwIndex)mxGetN(prhs[1]);
    nj = (mwIndex)mxGetM(prhs[2])*(mwIndex)mxGetN(prhs[2]);
    nv = (mwIndex)mxGetM(val)*(mwIndex)mxGetN(val);
    
    /* They must be equal */
    if ((ni!=l) || (ni!=l) || (nv!=l))
        mexErrMsgTxt("i and j and val must have the same number of elements.");
    
    /* They must be equal */
    if ((IdxClass = mxGetClassID(prhs[1])) != mxGetClassID(prhs[2]))
        mexErrMsgTxt("i and j must have the same class.");
     
    /* Get class, number of colmns, and number of non-zero elements */
    columns = mxGetN(Sin);
    rows =  mxGetM(Sin);
    nnz = *(mxGetJc(Sin) + columns);
    nzmax = mxGetNzmax(Sin);
    /* Pointers on data */
    DataClass =  mxGetClassID(Sin);
    /* Check with val and S have the same class */
    ValClass = mxGetClassID(val);
    /* They must be equal */
    if (ValClass != DataClass)
    {
        if ((DataClass != mxLOGICAL_CLASS) ||
            (ValClass != mxUINT8_CLASS) && (ValClass != mxINT8_CLASS))
        {
            mexErrMsgTxt("S and val must have the same class.");
        }
    }
       
    /* Vector of index */
    list = AllocatePos(l, IdxClass, mxGetData(prhs[1]), mxGetData(prhs[2]));

    /* Sort the indexes */
    MainSort64(list, l);
    
    /* Find the existing elements in Sin */
    buffer = spfind(Sin);

    /* Merge with new elements */
    Sout = MergeIdx(list, pos, l, buffer, nnz,
                    rows, columns, nzmax,
                    DataClass,
                    mxGetPr(val), mxGetPr(Sin),
                    mxGetPi(val), mxGetPi(Sin));
    
    /* Remove explicite zeros (if any) */
    SpRemoveZ(Sout);
    
    /* Assign the sparse result to first output */
    plhs[0] = Sout;

    /* Free buffer */
    buffer=MyFree(buffer);

    FreePos();
            
    return;
    
} /* setspvalmex */

