#include <R.h>
#include "methas.h"

void fexitc(const char *msg);

extern Cifns AreaIntCifns, BadGeyCifns, DgsCifns, DiggraCifns, 
             GeyerCifns, LookupCifns, SoftcoreCifns, 
             StraussCifns, StraussHardCifns, 
             MultiStraussCifns, MultiStraussHardCifns;

Cifns NullCifns = NULL_CIFNS;

typedef struct CifPair {
  char *name;
  Cifns *p;
} CifPair;

CifPair CifTable[] = { 
  {"areaint",  &AreaIntCifns},
  {"badgey",   &BadGeyCifns},
  {"dgs",      &DgsCifns},
  {"diggra",   &DiggraCifns},
  {"geyer",    &GeyerCifns},
  {"lookup",   &LookupCifns},
  {"sftcr",    &SoftcoreCifns},
  {"strauss",  &StraussCifns},
  {"straush",  &StraussHardCifns},
  {"straussm", &MultiStraussCifns},
  {"straushm", &MultiStraussHardCifns},
  {(char *) NULL, (Cifns *) NULL}
};

Cifns getcif(cifname) 
     char *cifname;
{
  int i;
  CifPair cp;
  for(i = 0; CifTable[i].name; i++) {
    cp = CifTable[i];
    if(strcmp(cifname, cp.name) == 0)
      return(*(cp.p));
  }
  fexitc("Unrecognised cif name; bailing out.\n");
  /* control never passes to here, but compilers don't know that */
  return(NullCifns);
}

/* R interface function, to check directly whether cif is recognised */

void knownCif(cifname, answer) 
     char** cifname;
     int* answer;
{
  int i;
  CifPair cp;
  for(i = 0; CifTable[i].name; i++) {
    cp = CifTable[i];
    if(strcmp(*cifname, cp.name) == 0) {
      *answer = 1;
      return;
    }
  }
  *answer = 0;
  return;
}
