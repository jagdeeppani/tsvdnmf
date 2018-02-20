This package allows to retreive and assign values of sparse matrix.

The indexes must be provided in (row,column) form. Functions in this
package are programmed with care so no overflow of linear index is used
with large-size matrix. I try to design an engine of automatic selection
algorithms so as to provide a good runtime performance across various
cases. I'm not sure it's worth the effort, but my experience shows that it
is easy to mess with the performance when non-appropriate method is using
in sparse manipulation. This package does not warranty for best performance,
but it just has a modest pretention to prevent user to use bad algorithm.
For a better customize solution, a calibration is required on the specific
computer by runing this command: spcalib('auto')

Contains of the package:

1. getsparse, spsubsref -> Get values of sparse matrix
2. setsparse, spsubsasgn -> Set values of sparse matrix
   (optional function handle can be provided)
3. spcalib -> Calibration routine
4. defaultTref -> Function generated automatically and needed by the
                  engine to select the a good strategy for accessing
5. defaultTrefOrg.m -> factory delivery of the above file
6. testsparse -> A test script of the package
7. bugsparseindexing.m -> A script to show Matlab bug
8. maxlinind -> Compute the safe limit of linear indexing not affected
                by the bug
9. getspvalmex.c -> mex code to retrieve the sparse values from
                    sub indexes
10. setspvalmex.c -> mex code to assign the sparse values from
                    sub indexes
11. spidx_install.m Installation script
12. sparsenc.m -> Non cummulative sparse creation

Original: 31/March/2009
01/April/2009; adjust the maximum of linear index allowed according to
               bug information communicated by The Mathworks
               (applicable only for 32-bit platform, v.2009A or previous)
03/April/2009: at last, a mex implementation that overcome the Bug
	       and improve the speed by approximatively a factor of 4.
05/April/2009: a mex sparse is developped to set values of existing
               non-zero elements. It does not seem faster than adding
               two sparse matrices! So we are not using it.
               A more comprehensive calibration, however not fully used.
11/April/2009: A new mex sparse is developped from scratch with quicksort.
               Finally, we beats builtin adding routine by about 40% of speed.
	       More comprehensive calibration routine is included.
15/April/2009: Correct a bug in setspvalmex.c, which concerns 64-bit only
17/April/2009: Correct bug when calling setspvalmex.c with empty indices
16/December/2009: Workaround for BUG of LCC (related to "|" operator in 
                  expression with int64)
13/November/2010: Automatic expansion of scalar row and column

To compile mex, (mex correctly setup), type this in command window 
> build_spidxmex

Manually user can type on 32-bit plateform
> mex -v -O getspvalmex.c
> mex -v -O setspvalmex.c

On 64-bit plateform, add oprion -largeArrayDims
> mex -v -O -largeArrayDims getspvalmex.c
> mex -v -O -largeArrayDims setspvalmex.c

Author: Bruno Luong <brunoluong@yahoo.com>