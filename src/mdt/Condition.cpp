#include "Condition.h"

Condition::Condition() {

  nullify();
  addToList();

}

Condition::Condition(char *filename) {

  nullify();
  FileIO *file = new FileIO(filename);
  file->getLine();
  file->readInt(&max_time_index);
  file->getLine();
  for (PetscInt i=0; i<3; i++)
    file->readDouble(&datum[i]);
  file->getLine();
  for (PetscInt i=0; i<3; i++)
    file->readDouble(&gradient[i]);
  times = new PetscReal[max_time_index];
  scalars = new PetscReal[max_time_index];
  for (PetscInt i=0; i<max_time_index; i++) {
    file->getLine();
    file->readDouble(&(times[i]));
    file->readDouble(&(scalars[i]));
  }
  delete file;
  PetscInt len = strcspn(filename,".");
  strncpy(name,filename,len);
  name[len] = '\0';
  addToList();

}

Condition::Condition(Condition *old_condition) {
  nullify();
  itype = old_condition->itype;
  ctype = new char[MAXWORDLENGTH];
  strcpy(old_condition->ctype,ctype);
  cur_value = old_condition->cur_value;
  for (PetscInt i=0; i<3; i++) {
    datum[i] = old_condition->datum[i];
    gradient[i] = old_condition->gradient[i];
  }
  id = old_condition->id;
  max_time_index = old_condition->max_time_index;
  times = new PetscReal[max_time_index];
  scalars = new PetscReal[max_time_index];
  for (PetscInt i=0; i<max_time_index; i++) {
    times[i] = old_condition->times[i];
    scalars[i] = old_condition->scalars[i];
  }
}

void Condition::nullify() {

  itype = -1;
  ctype = new char[MAXWORDLENGTH];
  cur_time_index = 0;
  cur_value = -999.;
  for (PetscInt i=0; i<3; i++) {
    datum[i] = 0.;
    gradient[i] = 0.;
  }
  max_time_index = 0;
  times = NULL;
  scalars = NULL;
  next = NULL;
  name[0] = '\0';

}

void Condition::addToList() {
  id = num_conditions;
  if (!list) list = this;
  if (end_of_list) end_of_list->next = this;
  end_of_list = this;
  num_conditions++;
}

void Condition::setId(PetscInt i) { id = i; }
void Condition::setType(PetscInt i) { itype = i; }
void Condition::setType(char *str) { 
  if (!ctype) ctype = new char[32]; 
  strcpy(ctype,str);
}
void Condition::setTimeUnit(char *str) { 
  strcpy(time_unit,str);
}
void Condition::setDatum(PetscReal *d) { 
  datum[0] = d[0];
  datum[1] = d[1];
  datum[2] = d[2];
}
void Condition::setGradient(PetscReal *d) { 
  gradient[0] = d[0];
  gradient[1] = d[1];
  gradient[2] = d[2];
}

PetscInt Condition::getId() { return id; }
PetscInt Condition::getType() { return itype; }
char *Condition::getTypePtr() { return &ctype[0]; }
char *Condition::getTimeUnitPtr() { return &time_unit[0]; }
PetscReal *Condition::getDatumPtr() { return &datum[0]; }
PetscReal *Condition::getGradientPtr() { return &gradient[0]; }
char *Condition::getName() { return &name[0]; }
Condition *Condition::getNext() { return next; };
 
void Condition::convertListToArray() {

  PetscInt count = 0;
  _array = new Condition *[num_conditions];

  Condition *cur_condition = list;
  while (cur_condition) {
    _array[count++] = cur_condition;
    cur_condition = cur_condition->next;
  }
}

void Condition::updateConditions(PetscReal time) {

  for (PetscInt icond=0; icond<num_conditions; icond++) {
    
    PetscInt cur_time_index = _array[icond]->cur_time_index;
    PetscInt next_time_index = min(cur_time_index+1,
                              _array[icond]->max_time_index);

    // ensure that condition has started
    if (time >= _array[icond]->times[cur_time_index] ||
        PetscAbsReal(time-_array[icond]->times[cur_time_index]) < 1.e-40) {
      // find appropriate time interval
      while (time >= _array[icond]->times[next_time_index] &&
             cur_time_index != next_time_index) {
        cur_time_index = next_time_index;
        // ensure that time index doesnot go beyond end of array
        if (next_time_index < _array[icond]->max_time_index)
          next_time_index++;
      }

      // interpolate value based on time
      if (cur_time_index < _array[icond]->max_time_index) {
        PetscReal time_fraction = 
          (time-_array[icond]->times[cur_time_index])/
          (_array[icond]->times[next_time_index]-
           _array[icond]->times[cur_time_index]);
        _array[icond]->cur_value = 
          (1.-time_fraction)*_array[icond]->scalars[cur_time_index] +
          time_fraction*_array[icond]->scalars[next_time_index];
      }
      else {
        _array[icond]->cur_value = 
          _array[icond]->scalars[_array[icond]->max_time_index];
      }
    }
  }
}

void Condition::initializeConditions() {

  for (PetscInt icond=0; icond<num_conditions; icond++) {

    if (!strcmp("dirichlet",_array[icond]->ctype))
      _array[icond]->itype = 1;
    else if (!strcmp("neumann",_array[icond]->ctype))
      _array[icond]->itype = 2;
    else if (!strcmp("hydraulic gradient",_array[icond]->ctype))
      _array[icond]->itype = 3;
    else if (!strcmp("seepage face",_array[icond]->ctype))
      _array[icond]->itype = 4;
    else {
      if (!myrank) cout << "Condition type: " << _array[icond]->ctype
        << " not supported." << endl;
      exit(0);
    }
  }
}

PetscReal Condition::computeHydrostaticPressure(PetscReal *coord) {
  PetscReal dx = coord[0]-datum[0];
  PetscReal dy = coord[1]-datum[1];
  PetscReal dz = coord[2]-datum[2];
  return cur_value + dx*gradient[0] + dy*gradient[1] + dz*gradient[2];
}

void Condition::printInfo() {
/*
    PetscInt idlocal;
  PetscReal area;
  PetscReal center[3],dist[3];
  PetscReal normal[3];
  char *type;
  PetscReal scalar;
  Condition *next;
*/  
  printf("\n  id: %5d\n",id);
  printf("type: %s scalar: %9.3e\n",ctype,scalars[0]);
  printf("  datumx: %9.3f datumy: %9.3f datumz: %9.3f\n",datum[0],datum[1],datum[2]);
  printf("   gradx: %9.3f  grady: %9.3f  gradz: %9.3f\n",gradient[0],gradient[1],gradient[2]);
}

void Condition::printConditions() {
  PetscErrorCode ierr;
 
  ierr = PetscSequentialPhaseBegin(PETSC_COMM_WORLD,1);
#include "Globals.h" 
  PetscPrintf(PETSC_COMM_WORLD,"\nBoundary Connections\n");
  printf("Processor[%d]: %d connections\n",myrank,num_conditions);
  Condition *cur_cond = list;
  while (cur_cond) {
    cur_cond->printInfo();
    cur_cond = cur_cond->next;
  }

  ierr = PetscSequentialPhaseEnd(PETSC_COMM_WORLD,1);
}

Condition::~Condition() {
  if (ctype) delete [] ctype;
  ctype = NULL;
  if (times) delete [] times;
  times = NULL;
  if (scalars) delete [] scalars;
  scalars = NULL;
}
