#ifdef HAVE_SSE2
inline
void fast_prod24 (boost::numeric::ublas::vector <double> *y,
                  const boost::numeric::ublas::matrix <double> &M,
                  const boost::numeric::ublas::vector <double> &x) {

    const double *Matptr = &M(0,0);
    const double *vecptr = &x(0);;
    double *outptr = const_cast <double *> (&(*y)[0]);

    asm (
         "movq %%rax, %%rcx;"
         "addq $4608, %%rcx;"
         "1:"
         // packed moves, loading the matrix row
         // packed multiplies, elementwise product with x
         // and start summing up all of it.
         "movapd    (%%rax), %%xmm0;"
         "mulpd     (%%rsi), %%xmm0;"
         "movapd  16(%%rax), %%xmm1;"
         "mulpd   16(%%rsi), %%xmm1;"
         "addpd      %%xmm0, %%xmm1;"
         "movapd  32(%%rax), %%xmm2;"
         "mulpd   32(%%rsi), %%xmm2;"
         "movapd  48(%%rax), %%xmm3;"
         "mulpd   48(%%rsi), %%xmm3;"
         "addpd      %%xmm2, %%xmm3;"
         "addpd      %%xmm1, %%xmm3;"
         "movapd  64(%%rax), %%xmm4;"
         "mulpd   64(%%rsi), %%xmm4;"
         "movapd  80(%%rax), %%xmm5;"
         "mulpd   80(%%rsi), %%xmm5;"
         "addpd      %%xmm4, %%xmm5;"
         "movapd  96(%%rax), %%xmm6;"
         "mulpd   96(%%rsi), %%xmm6;"
         "movapd 112(%%rax), %%xmm7;"
         "mulpd  112(%%rsi), %%xmm7;"
         "addpd      %%xmm6, %%xmm7;"
         "addpd      %%xmm5, %%xmm7;"
         "movapd 128(%%rax), %%xmm8;"
         "mulpd  128(%%rsi), %%xmm8;"
         "movapd 144(%%rax), %%xmm9;"
         "mulpd  144(%%rsi), %%xmm9;"
         "addpd      %%xmm8, %%xmm9;"
         "movapd 160(%%rax), %%xmm10;"
         "mulpd  160(%%rsi), %%xmm10;"
         "movapd 176(%%rax), %%xmm11;"
         "mulpd  176(%%rsi), %%xmm11;"
                                             // done in parallel
                                             "addq $192, %%rax;"
         "addpd     %%xmm10, %%xmm11;"
         // collect everything
         "addpd      %%xmm9, %%xmm11;"
         "addpd      %%xmm3, %%xmm7;"
         "addpd      %%xmm7, %%xmm11;"
         // unpack data
         "movq      %%xmm11, %%xmm12;"
         "unpckhpd  %%xmm11, %%xmm11;"
         "addsd     %%xmm11, %%xmm12;"
         // store away
         "movq      %%xmm12, (%%rdi);"
         "addq $8, %%rdi;"
         // loop test
         "cmpq        %%rax, %%rcx;"
         "jne 1b;"
         // restore pointers
         "addq $-192, %%rdi;"
         "addq $-4608, %%rax;"
          : // no output
          : // input
            "a"(Matptr), "S"(vecptr), "D"(outptr)
          : // clobbered registers
            "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
            "%xmm5", "%xmm6", "%xmm7", "%xmm8", "%xmm9",
            "%xmm10", "%xmm11", "%xmm12", "%rcx"
    );
}

#else // HAVE_SSE2

inline
void fast_prod24 (boost::numeric::ublas::vector <double> *y,
                  const boost::numeric::ublas::matrix <double> &M,
                  const boost::numeric::ublas::vector <double> &x) {
    noalias (*y) = prod (M, x);
}

#endif // HAVE_SSE2

inline
void   dgemv_24  (boost::numeric::ublas::vector <double> *y,
                  const boost::numeric::ublas::matrix <double> &M,
                  const boost::numeric::ublas::vector <double> &x) {
    fast_prod24 (y, M, x);
}
