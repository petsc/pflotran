#ifndef FILEIO_H_
#define FILEIO_H_

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#include "petscsys.h"

#define MAXCARDLENGTH 5 // add 1 to account of end of line \0
#define MAXWORDLENGTH 33
#define MAXSTRINGLENGTH 1025

class FileIO {
public:

  FileIO(char *filename);
  virtual ~FileIO();


  PetscInt getLine();
  PetscInt getInputLine();
  PetscInt readDouble(PetscReal *d);
  PetscReal readDoubleFast();
  PetscInt readInt(PetscInt *i);
  PetscInt readWord(char *word);
  PetscInt readQuotedWords(char *words);
  static PetscInt removeQuotes(char *str);
  PetscInt comparesTo(char *str);
  PetscInt startsWith(char *str);
  PetscInt findStringInFile(char *card);
  static void checkDefaultMessage(char *word, PetscErrorCode *ierr);
  static void checkErrorMessage(char *word1, char *word2, PetscErrorCode ierr);
  static void checkLineErrorMessage(char *word, PetscErrorCode ierr);
  static void toLower(char *word);
  static void toUpper(char *word);

  fstream file;
  stringstream *buffer;

};

#endif /*FILEIO_H_*/
