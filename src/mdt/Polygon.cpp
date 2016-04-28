#include "Polygon.h"

Polygon::Polygon() {
  num_points = 0;
  x = NULL;
  y = NULL;
}

void Polygon::createRiverEdgePolygon() {

  num_points = 8;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 594811.; y[n++] = 114864.;
  x[n] = 594636.; y[n++] = 115636.;
  x[n] = 594578.; y[n++] = 115817.;
  x[n] = 594413.; y[n++] = 116531.;
  x[n] = 594351.; y[n++] = 116922.;
  x[n] = 594317.; y[n++] = 117316.;
  x[n] = 592637.; y[n++] = 117114.;
  x[n] = 593288.; y[n++] = 114309.;

}

void Polygon::createNorthPondWestTrPolygon() {

  num_points = 4;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 594149.71; y[n++] = 116442.83;
  x[n] = 594169.08; y[n++] = 116447.74;
  x[n] = 594164.31; y[n++] = 116467.14;
  x[n] = 594144.91; y[n++] = 116462.3;

}

void Polygon::createNorthPondEastTrPolygon() {

  num_points = 4;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 594178.86; y[n++] = 116450.12;
  x[n] = 594198.18; y[n++] = 116455.04;
  x[n] = 594193.44; y[n++] = 116474.39;
  x[n] = 594173.93; y[n++] = 116469.47;

}

void Polygon::createIFCPolygon() {

  num_points = 3;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 594234.56; y[n++] = 116094.49;
  x[n] = 594239.81; y[n++] = 116034.77;
  x[n] = 594287.74; y[n++] = 116068.8;

}

void Polygon::createSPPPolygon() {

  num_points = 6;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 594270.; y[n++] = 116219.;
  x[n] = 594192.; y[n++] = 116173.;
  x[n] = 594229.; y[n++] = 116059.;
  x[n] = 594327.; y[n++] = 116057.;
  x[n] = 594344.; y[n++] = 116084.;
  x[n] = 594338.; y[n++] = 116198.;

}


void Polygon::createSPPPolygonCorrect() {

  num_points = 6;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 594298.; y[n++] = 116214.;
  x[n] = 594192.; y[n++] = 116155.;
  x[n] = 594242.; y[n++] = 116022.;
  x[n] = 594356.; y[n++] = 116022.;
  x[n] = 594379.; y[n++] = 116053.;
  x[n] = 594376.; y[n++] = 116184.;

}

void Polygon::createPlumePolygon() {

  num_points = 4;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
#if 0
  // old 50x50m plume
  x[n] = 594260.9; y[n++] = 115996.65;
  x[n] = 594309.51; y[n++] = 116008.59;
  x[n] = 594297.43; y[n++] = 116057.06;
  x[n] = 594248.82; y[n++] = 116045.26;
#else
  // new 60x60m plume at center ifc
  x[n] = 594241.65; y[n++] = 116032.95;
  x[n] = 594299.80; y[n++] = 116047.45;
  x[n] = 594285.40; y[n++] = 116105.75;
  x[n] = 594227.15; y[n++] = 116091.20;
#endif


}

void Polygon::createMidIFCPolygon() {

  num_points = 6;
  x = new PetscReal[num_points];
  y = new PetscReal[num_points];

  PetscInt n = 0;
  x[n] = 1.; y[n++] = 1.;

}

PetscInt Polygon::pointInPolygon(PetscReal x_, PetscReal y_) {

  PetscInt inside = 0;
  PetscInt j = 0;
  for (PetscInt i=0; i<num_points; i++) {
    if (++j == num_points) j = 0;
    if ((y[i] < y_ && y[j] >= y_) || (y[j] < y_ && y[i] >= y_)) {
      if (x[i] + (y_-y[i])/(y[j]-y[i])*(x[j]-x[i]) < x_)
        inside = !inside;
    }
  }
  return inside;
}

Polygon::~Polygon() {
  if (x) delete [] x;
  if (y) delete [] y;
}
