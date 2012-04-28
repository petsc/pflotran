#include "FileIO.h"

FileIO::FileIO(char *filename) {

  PetscPrintf(PETSC_COMM_WORLD,"%s\n",filename);
  file.open(filename,fstream::in);
  if (!file.is_open()) cout << "ERROR opening file " << filename << ".\n";
  buffer = NULL;

}

PetscInt FileIO::getLine() {

  PetscInt ierr = 0;
  delete buffer;

  while(1) {
    string s;
    getline(file,s,'\n');
    if (s.compare(0,1,":") && s.compare(0,1,"!")) {
      buffer = new stringstream(s.c_str());
      break;
    }
  }
  return file.eof() ? 0 : 1;

}

PetscInt FileIO::getInputLine() {

  return getLine();

}

PetscInt FileIO::readDouble(PetscReal *d) {

  *buffer >> *d;
  return buffer->fail() ? 0 : 1;
  
}

PetscInt FileIO::readInt(PetscInt *i) {

  *buffer >> *i;

  return buffer->fail() ? 0 : 1;
  
}

PetscInt FileIO::readWord(char *word) {

  /* Remove any preceding spaces(32), tabs(9), or commas(44) etc */
  char c;
  *buffer >> noskipws >> c;
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;

  /* Copy next group of chars to string */
  string str;
  while (c != 32 && c != 44 && c != 9 && c != '\0') {
    str.append(1,c);
    *buffer >> noskipws >> c;
  }
  strcpy(word,str.c_str());

  /* Remove any trailing spaces or commas etc */
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;
  buffer->unget();

  if (strlen(word) == 0) return 1;
  else return 0;

}

PetscInt FileIO::readQuotedWords(char *words) {

  /* Remove any preceding spaces(32), tabs(9), or commas(44) etc */
  char c;
  *buffer >> noskipws >> c;
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;

  string str;
  if (c == 34) { // quote found
    while (c != 34) {
      str.append(1,c);
      *buffer >> noskipws >> c;
    }
  }
  else {
    while (c != 32 && c != 44 && c != 9 && c != '\0') {
      str.append(1,c);
      *buffer >> noskipws >> c;
    }
  }
  strcpy(words,str.c_str());

  /* Remove any trailing spaces or commas etc */
  while ((c == 32 || c == 44 || c == 9) && c != '\0')
    *buffer >> noskipws >> c;
  buffer->unget();

  if (strlen(words) == 0) return 1;
  else return 0;

}

PetscInt FileIO::removeQuotes(char *str) {

  /* Remove all quotes */
  string str2;
  while (1) {
    size_t found = str2.find("\"");
    if (found == string::npos) break;
    else str2.erase(found);
  }
  strcpy(str,str2.c_str());

  return 0;

}

PetscInt FileIO::findStringInFile(char *card) {

  PetscInt ierr = 0;

  file.seekg(0,ios::beg);
  size_t len = strlen(card);
  size_t found = 0;
  while((ierr = getLine()) != 1) {
    string str = buffer->str();
    found = str.find(card);
    if (found!=string::npos){
      break;
    }
  }
  //return found ? 0 : 1;
  return found==string::npos ? 0 : 1;
}

PetscInt FileIO::comparesTo(char *str) {
  string str2 = buffer->str();
  return str2.compare(str);
}

PetscInt FileIO::startsWith(char *str) {
  return comparesTo(str);
}

void FileIO::checkDefaultMessage(char *word, PetscErrorCode *ierr) {

  if (ierr) cout << "\"" << word << "\" set to default value" << endl;
  *ierr = 0;

}

void FileIO::checkErrorMessage(char *word1, char *word2, PetscErrorCode ierr) {

  if (ierr) {
    cout << "Error reading \"" << word1 << "\" under keyword \"" << word2 << 
            "\"." << endl;
    exit(1);
  }

}

void FileIO::checkLineErrorMessage(char *word, PetscErrorCode ierr) {

  if (ierr) {
    cout << "Error reading in string in \"" << word << "\"." << endl;
    exit(1);
  }

}

void FileIO::toUpper(char *str) {

  PetscInt len = (PetscInt)strlen(str);
  for (PetscInt i=0; i < len; i++) {
    str[i] = toupper(str[i]);
  }

}

void FileIO::toLower(char *str) {

  PetscInt len = (PetscInt)strlen(str);
  for (PetscInt i=0; i < len; i++) {
    str[i] = tolower(str[i]);
  }

}

FileIO::~FileIO() {
  file.close();
}
