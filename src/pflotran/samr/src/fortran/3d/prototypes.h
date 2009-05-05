extern "C"{

   int samrccellmatmult3d_(
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&,
      const int&,
      const double*,
      const int&,
      const double*,
      const int&,
      const double*);

   int samrccellmatdiagscale3d_(
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&,
      const int&,
      const double*,
      const int&,
      const double*,
      const int&,
      const double*);


   int samrccellmatdiagscalelocal3d_(
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&,
      const int*,
      const int&,
      const double*,
      const int&,
      const double*);

   int samrapply7ptstencil3d_(
      const int&, const int&, const int&, const int&, const int&, const int&,
      const double*,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const double*,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const double*,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const double*);

   int pflotranpcflux3d_(
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&,
      const int&,
      const double*,
      const double*,
      const double*,
      const double*,
      const double*);


   int pflotranpcapply3d_(
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&,
      const int&,
      const double*,
      const double*,
      const double*, const double*, const double*,
      const double*,
      const double*,
      const double*);

   int samrsetjacobiansrccoeffs3d_(      
      const int&, const int&, const int&, const int&, const int&, const int&,
      const int&,
      const double*,
      const int&, const int&, const int&, const int&, const int&, const int&,
      const double*);
}
