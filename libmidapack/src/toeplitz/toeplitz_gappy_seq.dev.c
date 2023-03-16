

//=========================================================================
// Alternave version of the sequential routine for the gaps - in dev
// Need to change the name to gap_reduce_matrix for more explicit purpose
// fd@apc
int gstmm(double **V, int n, int m, int id0, int l, fftw_complex *T_fft,
          int lambda, fftw_complex *V_fft, double *V_rfft, fftw_plan plan_f,
          fftw_plan plan_b, int blocksize, int nfft, int *id0gap, int *lgap,
          int ngap) {

    // routine variable
    int i, j, k, p;                   // loop index
    int idf     = id0 + l - 1;
    int cfirst  = id0 / n;            // first column index
    int clast   = idf / n;            // last column index
    int clast_1 = (idf + 1) / n;
    int m_eff   = clast - cfirst + 1; // number of columns
    int rfirst  = id0 % n;
    int rlast   = idf % n;

    if (l < lambda) // test to avoid communications errors
        return print_error_message(1, __FILE__, __LINE__);

    int     nnew    = 0;
    int     lcolnew = 0;
    int     lcol;
    double *Vcol;
    int     icol, igap;
    ;
    int id0col;
    int id0out = 0;
    int lnew;

    int nfullcol;
    nfullcol = (l - (n - id0 % n) - (id0 + l) % n)
             / n; // check how many full columns input data have

    int ifirstgap, ilastgap;

    if (id0 % n != 0) {
        // first column if not full
        id0out = 0;
        id0col = id0 % n;
        Vcol   = (*V) + id0col;
        lcol   = n - id0col;

        // find the first and last gap on the column
        for (k = 0; k < ngap; k++) {
            if (lgap[k] != 0 && id0col < id0gap[k] + lgap[k]) break;
        }
        ifirstgap = k;

        for (k = ngap - 1; k >= 0; k--) {
            if (lgap[k] != 0 && id0gap[k] <= id0col) break;
        }
        ilastgap = k;


        for (i = 0; i < ngap; i++) printf("id0gap[%d]=%d\n", i, id0gap[i]);
        for (i = 0; i < ngap; i++) printf("lgap[%d]=%d\n", i, lgap[i]);

        printf("---\n");

        printf("lcol=%d\n", lcol);
        for (i = 0; i < lcol; i++) printf("Vcol[%d]=%f\n", i, Vcol[i]);

        // reduce the gap to lambda zeros on the column
        gap_reduce(&Vcol, id0col, lcol, lambda, id0gap, lgap, ngap, &lcolnew,
                   id0out);

        id0out += lcolnew;

        printf("---\n");

        for (i = 0; i < lcol; i++) printf("Vcolout[%d]=%f\n", i, Vcol[i]);
        printf("===\n");

        printf("---\n");

        for (i = 0; i < l; i++) printf("(*V)[%d]=%f\n", i, (*V)[i]);
    }

    // generic loop on the full column
    for (icol = 0; icol < nfullcol; icol++) {

        id0out = lcolnew;
        id0col = 0;
        Vcol   = (*V) + id0col;
        lcol   = n - id0col;

        gap_reduce(&Vcol, id0col, lcol, lambda, id0gap, lgap, ngap, &lcolnew,
                   id0out);

        id0out += lcolnew;

        if (icol == 0) // take the nnew directly if they is a full column
            nnew = lcolnew;
    }


    // last column if not full and more than one column
    if (idf % n != n && m > 1) {

        id0out = 0;
        id0col = 0;
        Vcol   = (*V) + id0col;
        lcol   = idf % n;

        // find the first and last gap on the column
        ifirstgap = 0; // we already know it, no need to compute

        for (k = ngap - 1; k >= 0; k--) {
            if (lgap[k] != 0 && id0gap[k] <= id0col) break;
        }
        ilastgap = k;

        // reduce the gap to lambda zeros on the column
        gap_reduce(&Vcol, id0col, lcol, lambda, id0gap, lgap, ngap, &lcolnew,
                   id0out);
    }


    // cleaning the extra terms
    if (nnew == 0) { // compute the nnew if there is no full column
        for (igap = 0; igap < ngap; igap++)
            nnew += id0gap[igap] - min(id0gap[igap], lambda);
    }

    for (i = nnew * m; i < n * m; i++) (*V)[i] = 0;

    printf("===\n");
    for (i = 0; i < l; i++) printf("(*V)[%d]=%f\n", i, (*V)[i]);


    // now do the product
    //   int id0new = ..;

    //  stmm( V, nnew, m, id0new, lnew, T_fft, lambda, V_fft, V_rfft, plan_f,
    //  plan_b, blocksize, nfft)


    // gap extend - same thing as above and reverse gap_reduce routine
    //  gap_extend()

    // maybe create gap_blockreduce to clarify this routine.

    return 0;
}

//=========================================================================

/// ...convert the data vector structure into a matrix structure optimized for
/// nfft
/** @ingroup group22
    ....Copy the data vector structure into an equivalent matrix with nfft
   column. Thus, the obtained matrix is optimize for the nfft multithreading
   algorithm use. The middle part is a direct copy of the data vector and we
   copy on the edges of each column the lambda terms needed to fullfill the
   correlation of theses data.
*/
int gap_reduce(double **V, int id0, int l, int lambda, int *id0gap, int *lgap,
               int ngap, int *newl, int id0out) {
    // routine variables
    int i, k, ii;
    int lg;
    int newid0gap_ip1;
    int newlgap_ip1;


    double *Vout;
    Vout = (*V) - (id0 - id0out);


    if (id0gap[0] > id0) {
        lg = id0gap[0] - id0;
        printf("lg=%d\n", lg);
        memmove(&(Vout)[0], &(*V)[0], lg * sizeof(double));
        // to do if the first element isn't a gap
    }

    int offset_first = id0gap[0] - id0 + min(lambda, lgap[0]);
    offset_first     = max(0, offset_first);


    //  for (k=0; k<offset_first; k++)
    //    (Vout)[k] = 0;

    printf("offset_first=%d\n", offset_first);

    for (i = 0; i < l; i++) printf("(Vout)[%d]=%f\n", i, (Vout)[i]);
    printf("---\n");

    printf("l=%d\n", l);


    for (i = 0; i < (ngap - 1); i++) {
        lg = id0gap[i + 1] - id0 - (id0gap[i] - id0 + lgap[i]);
        printf("lg=%d\n", lg);
        memmove(&(Vout)[offset_first], &(*V)[id0gap[i] - id0 + lgap[i]],
                lg * sizeof(double));

        for (ii = 0; ii < l; ii++) printf("(Vout)[%d]=%f\n", ii, (Vout)[ii]);

        newid0gap_ip1 = offset_first + lg;
        newlgap_ip1   = min(lambda, lgap[i + 1]);
        for (k = newid0gap_ip1; k < newid0gap_ip1 + newlgap_ip1; k++)
            (Vout)[k] = 0;

        printf("lg=%d ; lambda=%d ; lgap[i+1]=%d\n", lg, lambda, lgap[i + 1]);
        offset_first += lg + min(lambda, lgap[i + 1]);
    } // end gaps loop

    printf("---\n");
    printf("offset_first=%d\n", offset_first);
    for (i = 0; i < l; i++) printf("(Vout)[%d]=%f\n", i, (Vout)[i]);

    i = (ngap - 1);

    if (id0gap[i] - id0 + lgap[i] < l) {
        printf("toc\n");
        lg = l - (id0gap[i] - id0 + lgap[i]);
        memmove(&(Vout)[offset_first], &(*V)[id0gap[i] - id0 + lgap[i]],
                lg * sizeof(double));
        offset_first += lg;
        *newl = offset_first;
    } else {
        *newl = offset_first - min(lambda, lgap[ngap - 1])
              + min(lambda, (id0 + l - id0gap[ngap - 1]));
    }

    printf("---\n");
    printf("l=%d, *newl=%d, offset_first=%d\n", l, *newl, offset_first);

    for (i = offset_first; i < l; i++) (Vout)[i] = 0;

    printf("*newl=%d\n", *newl);


    return 0;
}
//=========================================================================

/// Reduce the vector and mask the defined gaps
/** @ingroup group11
   \param V input/output vector
   \param l length of the vector
   \param id0gap gap first index
   \param lgap gap length
   \param ngap number of gaps
*/
int gap_masking(double **V, int l, int *id0gap, int *lgap, int ngap) {

    int i, k;
    int lg;

    int offset_first = id0gap[0];
    for (i = 0; i < (ngap - 1); i++) {
        lg = id0gap[i + 1] - (id0gap[i] + lgap[i]);
        memmove(&(*V)[offset_first], &(*V)[id0gap[i] + lgap[i]],
                lg * sizeof(double));
        offset_first += lg;
    }

    i  = (ngap - 1);
    lg = l - (id0gap[i] + lgap[i]);

    memmove(&(*V)[offset_first], &(*V)[id0gap[i] + lgap[i]],
            lg * sizeof(double));
    offset_first += lg;


    for (i = offset_first; i < l; i++) (*V)[i] = 0;


    return 0;
}


//=========================================================================

/// Extend the vector and add zeros on the gaps locations
/** @ingroup group11
   \param V input/output vector
   \param l length of the extend vector
   \param id0gap gap first index
   \param lgap gap length
   \param ngap number of gaps
*/
int gap_filling(double **V, int l, int *id0gap, int *lgap, int ngap) {

    int i, k;
    int lg;

    int lgaptot = 0;
    for (i = 0; i < ngap; i++) lgaptot += lgap[i];

    int id0_Vmv;
    int offset_last;

    id0_Vmv     = id0gap[ngap - 1] + lgap[ngap - 1];
    offset_last = id0_Vmv - lgaptot;
    lg          = l - id0_Vmv;

    if (lg > 0)
        memmove(&(*V)[id0_Vmv], &(*V)[offset_last], lg * sizeof(double));

    for (i = ngap - 2; i >= 0; i--) {
        id0_Vmv = id0gap[i] + lgap[i];
        lg      = id0gap[i + 1] - id0_Vmv;
        offset_last -= lg;

        memmove(&(*V)[id0_Vmv], &(*V)[offset_last], lg * sizeof(double));
    }

    for (i = 0; i < ngap; i++)
        for (k = 0; k < lgap[i]; k++) (*V)[id0gap[i] + k] = 0;

    return 0;
}


// a naive version - in dev
int gap_masking_naive(double **V, int id0, int local_V_size, int m, int nrow,
                      int *id0gap, int *lgap, int ngap) {
    int i, j, k;
    int icol;

    int offsetV = id0;
    k           = 0;

    for (i = 0; i < local_V_size; i++) {
        icol = (i + id0) % nrow;

        if (icol < id0gap[k]) {
            (*V)[offsetV] = (*V)[i];
            offsetV       = offsetV + 1;
        } else if (icol >= id0gap[k] + lgap[k]) {
            k = k + 1;
            i = i - 1;
            //    (*V)[offsetV]=(*V)[i];
            //    offsetV=offsetV+1;
        } else { // do not copy, just skip
        }


        return 0;
    }
