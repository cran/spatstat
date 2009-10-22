#include <R.h>
#include "methas.h"

void fexitc(const char *msg);

extern Cifns AreaIntCifns, BadGeyCifns, DgsCifns, DiggraCifns, 
             GeyerCifns, LookupCifns, SoftcoreCifns, 
             StraussCifns, StraussHardCifns, 
             MultiStraussCifns, MultiStraussHardCifns;

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

Cifns getcifold(cifname) 
     char *cifname;
{
  if(strcmp(cifname, "strauss") == 0) 
    return(StraussCifns);
  if(strcmp(cifname, "straush") == 0) 
    return(StraussHardCifns);
  if(strcmp(cifname, "straussm") == 0) 
    return(MultiStraussCifns);
  if(strcmp(cifname, "straushm") == 0) 
    return(MultiStraussHardCifns);
  if(strcmp(cifname, "areaint") == 0)
    return(AreaIntCifns);
  if(strcmp(cifname, "geyer") == 0)
    return(GeyerCifns);
  if(strcmp(cifname, "dgs") == 0)
    return(DgsCifns);
  if(strcmp(cifname, "diggra") == 0)
    return(DiggraCifns);
  if(strcmp(cifname, "sftcr") == 0)
    return(SoftcoreCifns);
  if(strcmp(cifname, "lookup") == 0)
    return(LookupCifns);
  if(strcmp(cifname, "badgey") == 0)
    return(BadGeyCifns);
  fexitc("Unrecognised cif name; bailing out.\n");
}
