#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

FILE *out;

void R_CheckUserInterrupt(void);

struct Point{ long int No; float X; float Y; float R; struct Point *next; }; 

struct Point2{ long int No; float X; float Y; 
  char InLower[2]; 
  double Beta; double TempBeta; struct Point2 *next; }; 

struct Point3{ char Case; char XCell; char YCell; struct Point3 *next; }; 

const float Pi=3.141593;

double slumptal(void){
  return(runif(0,1));
}

void normal2(double *x, double *y, double sigma){

  double v1,v2,rsq,fac;

  do{
    v1=2*slumptal()-1.0;
    v2=2*slumptal()-1.0;
    rsq=v1*v1+v2*v2;
  }while (rsq >= 1.0 || rsq == 0.0);
  fac=sqrt(-2.0*log(rsq)/rsq);

  *x = v2*fac*sigma;
  *y = v1*fac*sigma;
}


long int poisson(double lambda){
  return((long int)rpois(lambda));
}

class Point2Pattern {
public:
  long int UpperLiving[2];
  long int MaxXCell, MaxYCell, NoP;
  double XCellDim, YCellDim, Xmin, Xmax, Ymin, Ymax;
  struct Point2 *headCell[10][10],*dummyCell;
  char DirX[10], DirY[10];
 
  Point2Pattern(double xmin, double xmax,
		double ymin, double ymax, 
		long int mxc, long int myc){
    long int i,j;
    UpperLiving[0] = 0;
    UpperLiving[1] = 0;
    Xmin = xmin; Xmax = xmax;
    Ymin = ymin; Ymax = ymax;
    DirX[1] = 1; DirY[1] = 0;
    DirX[2] = 1; DirY[2] = -1;
    DirX[3] = 0; DirY[3] = -1;
    DirX[4] = -1; DirY[4] = -1;
    DirX[5] = -1; DirY[5] = 0;
    DirX[6] = -1; DirY[6] = 1;
    DirX[7] = 0; DirY[7] = 1;
    DirX[8] = 1; DirY[8] = 1;    
    NoP = 0;
    dummyCell = (struct Point2 *) malloc(sizeof *dummyCell);
    if(dummyCell == NULL) printf("error: malloc returned NULL\n");
    dummyCell->next = dummyCell;
    dummyCell->No = 0;
    MaxXCell = mxc; MaxYCell = myc;
    if(MaxXCell>9) MaxXCell = 9;
    if(MaxYCell>9) MaxYCell = 9;
    for(i=0;i<=MaxXCell;i++){
      for(j=0;j<=MaxYCell;j++){
	headCell[i][j] = 
	  (struct Point2 *) malloc(sizeof *headCell[i][j]);
	if(headCell[i][j] == NULL) printf("error: malloc returned NULL\n");
	headCell[i][j]->next=dummyCell;
      }
    }
    XCellDim = (Xmax-Xmin)/((double)(MaxXCell+1));
    YCellDim = (Ymax-Ymin)/((double)(MaxYCell+1));
  };
  ~Point2Pattern(){}
  void Print();
  void Return(double *X, double *Y, int *num, int maxnum);
  long int Count();
  void Empty();
  void Clean();
  void DumpToFile(char FileName[100]);
  void ReadFromFile(char FileName[100]);
};

void Point2Pattern::Print(){
  long int i,j,k;
  k =0;
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	k++;
	printf("%f %f %d %d %d=%d %d=%d UL0 %d UL1 %d %f\n",
	       TempCell->X,TempCell->Y,k,
	       TempCell->No,
	       i,int(TempCell->X/XCellDim),
	       j,int(TempCell->Y/YCellDim),
	       TempCell->InLower[0],TempCell->InLower[1],
	       TempCell->Beta);
	TempCell = TempCell->next;
      }
    }
  }
  printf("Printed %d points.\n",k);
}

void Point2Pattern::Return(double *X, double *Y, int *num, int maxnum){
  long int i,j,k;
  k =0; *num = 0;
  if(UpperLiving[0]<=maxnum){
    struct Point2 *TempCell;
    for(i=0;i<=MaxXCell;i++){
      for(j=0;j<=MaxYCell;j++){
	//printf("%d %d:\n",i,j);
	TempCell = headCell[i][j]->next;
	while(TempCell->next != TempCell){
	  X[k] = TempCell->X;
	  Y[k] = TempCell->Y;	
	  k++;
	  TempCell = TempCell->next;
	}
      }
    }    
    *num = k;
  }else{
    *num = -1;
  }
}



long int Point2Pattern::Count(){
  long int i,j,k;
  k =0;
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	k++;
	TempCell = TempCell->next;
      }
    }
  }
  //printf("Printed %d points.\n",k);
  return(k);
}


void Point2Pattern::Empty(){
  struct Point2 *TempCell, *TempCell2;
  long int i,j;
  
  //k=0;

  for(i=0; i<=this->MaxXCell; i++){
    for(j=0; j<=this->MaxYCell; j++){
      TempCell = headCell[i][j]->next;
      while(TempCell!=TempCell->next){	
	//k++; printf("%d %d %d\n",i,j,k);
	TempCell2 = TempCell->next;
	free(TempCell);
	TempCell = TempCell2;
      }
      headCell[i][j]->next = dummyCell;
    }
  }
}

void Point2Pattern::Clean(){
  struct Point2 *TempCell, *TempCell2;
  long int i,j;
  
  for(i=0; i<=MaxXCell; i++){
    for(j=0; j<=MaxYCell; j++){
      TempCell = headCell[i][j];
      TempCell2 = headCell[i][j]->next;
      while(TempCell2!=TempCell2->next){
	TempCell2->No = 0;
	if(TempCell2->InLower[0]==0){
	  TempCell->next = TempCell2->next;
	  free(TempCell2);
	  TempCell2 = TempCell->next;
	}
	else{
	  TempCell2 = TempCell2->next;
	  TempCell = TempCell->next;
	}
      }
    }
  }
}

void Point2Pattern::DumpToFile(char FileName[100]){
  FILE *out;
  long int i,j;
  out = fopen(FileName,"w");
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	fprintf(out,"%f\t%f\t%d\n",
	       TempCell->X,TempCell->Y,TempCell->No);
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}


void Point2Pattern::ReadFromFile(char FileName[100]){
  FILE *out;
  long int i,j,k,XCell,YCell;
  float f1,xs,ys;
  out = fopen(FileName,"r");
  struct Point2 *TempCell;
  k=0;
  while(feof(out)==0){
    k++;
    fscanf(out,"%f%f\n",&xs,&ys);
    //printf("%f %f\n",xs,ys);
    TempCell = (struct Point2 *) malloc(sizeof *TempCell);
    if(TempCell==NULL) printf("error: malloc returned NULL\n");
    TempCell->No = k;
    TempCell->X = xs;
    TempCell->Y = ys;
    TempCell->InLower[0] = 1;
    TempCell->InLower[1] = 1;

    f1 = (xs-Xmin)/XCellDim;  XCell = int(f1);
    if(XCell>MaxXCell) XCell = MaxXCell;
    f1 = (ys-Ymin)/YCellDim;  YCell = int(f1);
    if(YCell>MaxYCell) YCell = MaxYCell;

    TempCell->next = headCell[XCell][YCell]->next;
    headCell[XCell][YCell]->next = TempCell;

  }
  fclose(out);
  printf("%d points loaded.\n",k);

}





class PointProcess {
 public:
  double Xmin, Xmax, Ymin, Ymax, TotalBirthRate, InteractionRange;
  PointProcess(double xmin, double xmax, double ymin, double ymax){
    Xmin = xmin; Xmax = xmax;
    Ymin = ymin; Ymax = ymax;
  }
  ~PointProcess(){}
  virtual void NewEvent(double *x, double *y, char *InWindow)=0;
  virtual void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP)=0;
  virtual double Interaction(double r)=0;
  virtual void CalcBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
    printf("Define CalcBeta...\n");
  }
  virtual void CheckBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
    printf("Define CheckBeta...\n");
  }
  virtual double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p)
  { return(0.0);};
  virtual double lnDens(Point2Pattern *p2p);
  virtual void Beta(struct Point2 *TempCell){
    TempCell->Beta = 0;
    printf("Define Beta...\n");};
};

double PointProcess::lnDens(Point2Pattern *p2p){  
  double f1;
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx;
  double dy,dx, lnDens,dst;
  struct Point2 *TempCell, *TempCell2;

  dx = (Xmax-Xmin)/(double(p2p->MaxXCell+1));
  dy = (Ymax-Ymin)/(double(p2p->MaxYCell+1));
  rx = int(InteractionRange/dx+1.0);
  ry = int(InteractionRange/dy+1.0);
  
  //printf("1:%f 2:%f 3:%d 4:%d 5:%f 6:%f\n",dx,dy,rx,ry,
  // this->InteractionRange,InteractionRange);
  //printf("mx:%d my:%d\n",p2p->MaxXCell,p2p->MaxYCell);

  lnDens = 0;

  //printf("lnDens: %f (0)\n",lnDens);
  
  for(xc = 0; xc <= p2p->MaxXCell; xc++){
    for(yc = 0; yc <= p2p->MaxYCell; yc++){
      //if(xc==1) printf("%d %d\n",xc,yc);
      TempCell = p2p->headCell[xc][yc]->next;
      while(TempCell != TempCell->next){
	lnDens += log(TempCell->Beta);
	//printf("lnDens: %f (1) %d %d %d %d Beta %f\n",lnDens,xc,yc,
	//       p2p->MaxXCell,p2p->MaxYCell,TempCell->Beta);
	//if(lnDens<(-100000)){printf("%f",lnDens); scanf("%f",&f1);}
	if(InteractionRange>0){
	  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
	  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
	  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
	  if((yc-ry)>=0) fy=yc-ry; else fy = 0;
	  for(xco = fx; xco <= tx; xco++){
	    for(yco = fy; yco <= ty; yco++){
	      //if(xc==1) printf("%d %d %d %d %d %d\n",xco,yco,fx,tx,fy,ty);
	      TempCell2 = p2p->headCell[xco][yco]->next;
	      while(TempCell2!=TempCell2->next){
		if(TempCell2 != TempCell){
		  dst = sqrt(pow(TempCell->X-TempCell2->X,2)+
			     pow(TempCell->Y-TempCell2->Y,2));
		  lnDens += log(Interaction(dst));
		}
		TempCell2 = TempCell2->next; 
	      }
	    }
	  }
	  //printf("lnDens: %f\n",lnDens);
	}
	TempCell = TempCell->next;
      }
    }
  }
  return(lnDens);

}




class StraussProcess : public PointProcess {
 public:
  double beta, gamma, R;
  StraussProcess(double xmin, double xmax, double ymin, double ymax, 
		double b, double g, double Ri);
  ~StraussProcess(){}
  void NewEvent(double *x, double *y, char *InWindow);
  void GeneratePoisson(Point *headPoint, 
			       long int *GeneratedPoints,
			       long int *LivingPoints,
			       long int *NoP);
  double Interaction(double r);
  void CalcBeta(long int xsidepomm, long int ysidepomm, 
	   double *betapomm);
  void CheckBeta(long int xsidepomm, long int ysidepomm, 
		 double *betapomm);
  double lnCondInt(struct Point2 *TempCell, Point2Pattern *p2p);
  void Beta(struct Point2 *TempCell);
  void CalcBeta(Point2Pattern *p2p);
};

StraussProcess::StraussProcess(double xmin, double xmax, 
			      double ymin, double ymax, 
			      double b, double g, double Ri) :
  PointProcess(xmin, xmax, ymin, ymax){
  beta = b; gamma = g; R = Ri; 
  InteractionRange = R;
  TotalBirthRate = beta*(xmax-xmin)*(ymax-ymin);
  }  

double StraussProcess::Interaction(double r)
{
  double rtn;
  rtn = 1;
  if(r<R) rtn = gamma;
  return(rtn);
}

void StraussProcess::NewEvent(double *x, double *y, char *InWindow)
{
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = slumptal()*Xdim+Xmin;
  *y = slumptal()*Ydim+Ymin;
  *InWindow = 1;
}

void StraussProcess::GeneratePoisson(Point *headPoint, 
			      long int *GeneratedPoints,
			      long int *LivingPoints,
			      long int *NoP)
{
  int i;
  double xtemp, ytemp, L, Xdim, Ydim;
  struct Point *TempPoint;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  L = beta*Xdim*Ydim;
  *GeneratedPoints = poisson(L);
  *LivingPoints = *GeneratedPoints;
  for (i=1; i<=*GeneratedPoints ; i++){
    //printf("Generating StraussProcess Poisson 3\n");
    //scanf("%f",&f1);
    xtemp = slumptal()*Xdim+Xmin;
    ytemp = slumptal()*Ydim+Ymin;
    //printf("Generating StraussProcess Poisson 3.2\n");
    TempPoint = (struct Point *) malloc(sizeof *TempPoint);
    if(TempPoint == NULL ) printf("error: malloc returned NULL\n");
    TempPoint->X = xtemp;
    TempPoint->Y = ytemp;
    TempPoint->No = i;
    TempPoint->R = slumptal();
    //printf("Generating StraussProcess Poisson 3.6\n");
    TempPoint->next = headPoint->next;
    headPoint->next = TempPoint;
    *NoP = *NoP + 1;
  }
};

void StraussProcess::CalcBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
  long int i,j,k;
  k=0;
  printf("\ndiagnostic message: Strauss CalcBeta... %d %d\n",xsidepomm,ysidepomm);
  for(i=0; i<xsidepomm; i++){
    for(j=0; j<ysidepomm; j++){
      *(betapomm + i*ysidepomm + j) = this->beta;
      k++;
    }
  } 
}

void StraussProcess::CheckBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
  long int i,j,k;
  double d1;
  k=0;
  printf("\ndiagnostic message: Strauss CalcBeta... %d %d\n",xsidepomm,ysidepomm);
  for(i=0; i<xsidepomm; i++){
    for(j=0; j<ysidepomm; j++){
      if((fabs(*(betapomm + i*ysidepomm + j)- beta)>0.001) && (k==0)){
	printf("%f %f %f %d %d\n",fabs(*(betapomm + i*ysidepomm + j)- beta),
	       *(betapomm + i*ysidepomm + j),beta,i,j);
	k++;
	scanf("%f",&d1);
      }
    }
  } 
}

double StraussProcess::lnCondInt(struct Point2 *TempCell, 
				 Point2Pattern *p2p){
  double f1;
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx,k;
  double dy,dx, lnCI,dst;
  struct Point2 *TempCell2;
  
  f1 = (TempCell->X-p2p->Xmin)/p2p->XCellDim;  xc = int(f1);
  if(xc>p2p->MaxXCell) xc = p2p->MaxXCell;
  f1 = (TempCell->Y-p2p->Ymin)/p2p->YCellDim;  yc = int(f1);
  if(yc>p2p->MaxYCell) yc = p2p->MaxYCell;
  
  dx = (Xmax-Xmin)/(double(p2p->MaxXCell+1));
  dy = (Ymax-Ymin)/(double(p2p->MaxYCell+1));
  rx = int(this->InteractionRange/dx+1.0);
  ry = int(this->InteractionRange/dy+1.0);
  
  lnCI = log(TempCell->Beta);

  k = 0;
  
  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
  if((yc-ry)>=0) fy=yc-ry; else fy = 0;

  //printf("MCI! %d %d %d %d\n",fx,tx,fy,ty);

  for(xco = fx; xco <= tx; xco++){
    for(yco = fy; yco <= ty; yco++){
      TempCell2 = p2p->headCell[xco][yco]->next;
      while(TempCell2!=TempCell2->next){
	if(TempCell2 != TempCell){
	  k++;
	  dst = sqrt(pow(TempCell->X-TempCell2->X,2)+
		     pow(TempCell->Y-TempCell2->Y,2));
	  lnCI += log(Interaction(dst));
	}
	TempCell2 = TempCell2->next; 
	
      }
    }
  }
  return(lnCI);
}

void StraussProcess::Beta(struct Point2 *TempCell){
  TempCell->Beta = beta;
}

void StraussProcess::CalcBeta(Point2Pattern *p2p){
  long int xco,yco,xcm,ycm,fx,tx,fy,ty,ry,rx;
  double dy,dx;
  struct Point2 *TempMother;

  for(xco = 0; xco <= p2p->MaxXCell; xco++){
    for(yco = 0; yco <= p2p->MaxYCell; yco++){
      TempMother = p2p->headCell[xco][yco]->next;
      while(TempMother!=TempMother->next){
	TempMother->Beta = this->beta;
	TempMother = TempMother->next;
      }
    }
  }
}


class Sampler{
 public:
  PointProcess *PP;
  Point2Pattern *P2P;
  long int GeneratedPoints, LivingPoints, NoP;
  //long int UpperLiving[2];
  Sampler(PointProcess *p){ PP = p;}
  ~Sampler(){}
  void Sim(Point2Pattern *p2p);
  long int BirthDeath(long int TimeStep,
		      struct Point *headLiving,
		      struct Point *headDeleted,
		      struct Point3 *headTransition);
  void Sampler::Forward(long int TS, long int TT, char TX, char TY,
		      struct Point *Proposal, long int *DDD);
};


void Sampler::Forward(long int TS, long int TT, char TX, char TY,
		      struct Point *Proposal, long int *DDD){

  long int i,j,k, NoN[2], XCell, YCell, XC, YC, DirectionN;
  double dtmp,dtmpx,dtmpy, tmpR, TempGamma[2], TempI;
  struct Point2 *TempCell, *TempCell2;
  float f1;

  /* Birth */
  if(TT==1){
    f1 = (Proposal->X-P2P->Xmin)/P2P->XCellDim;  XCell = int(f1);
    if(XCell>P2P->MaxXCell) XCell = P2P->MaxXCell;
    f1 = (Proposal->Y-P2P->Ymin)/P2P->YCellDim;  YCell = int(f1);
    if(YCell>P2P->MaxYCell) YCell = P2P->MaxYCell;

    TempCell = (struct Point2 *) malloc(sizeof *TempCell);
    if(TempCell==NULL) printf("error: malloc returned NULL\n");
    TempCell->No = Proposal->No;
    TempCell->X = Proposal->X;
    TempCell->Y = Proposal->Y;

    tmpR = Proposal->R;
    TempCell->next = P2P->headCell[XCell][YCell]->next;
    P2P->headCell[XCell][YCell]->next = TempCell;
    TempCell->InLower[0]=0;
    TempCell->InLower[1]=0;

    TempGamma[0] = 1.0; TempGamma[1] = 1.0;    

    /*same cell*/
    TempCell2 = TempCell->next; 
    while(TempCell2 != TempCell2->next){
      dtmpx = TempCell->X - TempCell2->X;
      dtmpy = TempCell->Y - TempCell2->Y;
      dtmp  = sqrt(dtmpx*dtmpx+dtmpy*dtmpy);      
      TempI = PP->Interaction(dtmp);
      if(TempCell2->InLower[0]==1) TempGamma[0] = TempGamma[0]*TempI;
      if(TempCell2->InLower[1]==1) TempGamma[1] = TempGamma[1]*TempI;
      TempCell2=TempCell2->next;
    }
    /*eight other cells*/
    for(DirectionN=1;DirectionN<=8;DirectionN++){
      if(((XCell+P2P->DirX[DirectionN])>=0) &&
	 ((XCell+P2P->DirX[DirectionN])<=P2P->MaxXCell) &&
	 ((YCell+P2P->DirY[DirectionN])>=0) &&
	 ((YCell+P2P->DirY[DirectionN])<=P2P->MaxYCell)){
	TempCell2 = 
	  P2P->headCell[XCell+P2P->DirX[DirectionN]]
	  [YCell+P2P->DirY[DirectionN]]->next;
	while(TempCell2!=TempCell2->next){
	  dtmpx = TempCell->X - TempCell2->X;
	  dtmpy = TempCell->Y - TempCell2->Y;
	  dtmp  = sqrt(dtmpx*dtmpx+dtmpy*dtmpy);      
	  TempI = PP->Interaction(dtmp);
	  if(TempCell2->InLower[0]==1) 
	    TempGamma[0] = TempGamma[0]*TempI;
	  if(TempCell2->InLower[1]==1) 
	    TempGamma[1] = TempGamma[1]*TempI;
	  TempCell2=TempCell2->next;
	}
      }
    }

    if(tmpR <= TempGamma[1] ){ 
      TempCell->InLower[0]=1;
      P2P->UpperLiving[0] = P2P->UpperLiving[0] +1;
    }
    if(tmpR <= TempGamma[0] ){ 
      TempCell->InLower[1]=1;
      P2P->UpperLiving[1] = P2P->UpperLiving[1] +1;
    }
  }
  /* Death */
  if(TT==0){
    TempCell=P2P->headCell[TX][TY];
    while(TempCell->next->No != *DDD){
      TempCell = TempCell->next;
      if(TempCell->next == TempCell) {
	printf("error: unexpected self-reference. Dumping...\n"); 
	P2P->Print(); break;
      }
    };
    if(*DDD!=TempCell->next->No) 
      printf("diagnostic message: multi cell:  !!DDD:%d TempUpper->No:%d ",
	     *DDD,TempCell->No);
    if(TempCell->next->InLower[0]==1)
      P2P->UpperLiving[0] = P2P->UpperLiving[0] -1;
    if(TempCell->next->InLower[1]==1) 
      P2P->UpperLiving[1] = P2P->UpperLiving[1] -1;
    TempCell2 = TempCell->next;
    TempCell->next = TempCell2->next;
    free(TempCell2);
    /* Common stuff */
    //KillCounter ++;
    *DDD = *DDD - 1;
  }
}


long int Sampler::BirthDeath(long int TimeStep,
		      struct Point *headLiving,
		      struct Point *headDeleted,
		      struct Point3 *headTransition){
  long int i,k,j,n;
  float f1,f2,f3,f4;
  double min,tmp,xtemp,ytemp;
  char InWindow, Success;
  struct Point *TempPoint, *TempPoint2;
  struct Point3 *TempTransition;

  R_CheckUserInterrupt();

  f1 = LivingPoints; f2 = PP->TotalBirthRate; f3 = f2/(f1+f2);
  f4 = slumptal();
  n = 0;
  Success = 0;

  //printf("LivingPoints: %d TotalBirthRate %f GeneratedPoints %d\n",
  // LivingPoints,PP->TotalBirthRate,GeneratedPoints);
  
  /* Birth */
  while(Success==0){
  if(f4<f3){
    //printf("Ping 1 (BD)\n");
    PP->NewEvent(&xtemp, &ytemp, &InWindow);
    //printf("Ping 2 (BD)\n");
    if(InWindow==1){
      Success = 1;
      TempTransition = (struct Point3 *) malloc(sizeof *TempTransition);
      if(TempTransition == NULL) printf("error: malloc returned NULL\n");
      //printf("Ping 3 (BD)\n");
      TempTransition->Case = 0;
      LivingPoints ++;
      GeneratedPoints ++;
      TempPoint = (struct Point *) malloc(sizeof *TempPoint);
      if(TempPoint == NULL) printf("error: malloc returned NULL\n");
      TempPoint->X = xtemp;
      TempPoint->Y = ytemp;
      TempPoint->No = GeneratedPoints;
      TempPoint->R = slumptal();
      TempPoint->next = headLiving->next;
      headLiving->next = TempPoint;
      NoP ++;
      f1 = (TempPoint->X-P2P->Xmin)/P2P->XCellDim;  
      TempTransition->XCell = int(f1);
      f1 = (TempPoint->Y-P2P->Ymin)/P2P->YCellDim;  
      TempTransition->YCell = int(f1);
      //printf("X %f XCell %d\n",TempPoint->X,TempTransition->XCell);
      if(TempTransition->XCell>P2P->MaxXCell){
	TempTransition->XCell=P2P->MaxXCell;
	printf("diagnostic message: random X value exceeded limit\n");
      }
      if(TempTransition->YCell>P2P->MaxYCell){ 
	TempTransition->YCell=P2P->MaxYCell;
	printf("diagnostic message: random Y value exceeded limit\n");
      }
      TempTransition->next = headTransition->next;
      headTransition->next = TempTransition;
    }
  }
  /* Death */
  else{
    Success = 1;
    TempTransition = (struct Point3 *) malloc(sizeof *TempTransition);
    if(TempTransition == NULL) printf("error: malloc returned NULL\n");
    TempTransition->Case = 1;
    f1 = LivingPoints; f2 = f1*slumptal()+1.0;
    n = int(f2);
    if(n>LivingPoints){
      printf("diagnostic message: random integer n=%d > %d = number of living points\n", n,LivingPoints);
      n=LivingPoints;
    }
    TempPoint = headLiving;
    for(i=1; i<=n; i++){ 
      TempPoint2 = TempPoint;
      TempPoint = TempPoint->next;
      }
    TempPoint2->next = TempPoint->next;
    
    TempPoint->next = headDeleted->next;  
    headDeleted->next = TempPoint;

    LivingPoints --;
    NoP --;
    TempTransition->next = headTransition->next;
    headTransition->next = TempTransition;
  }
  }
  return(n);
}

void Sampler::Sim(Point2Pattern *p2p) {

  P2P = p2p;
  long int StartTime, EndTime, TimeStep, D0Time, D0Living;
  long int XCell, YCell, DDD, i;
  float f1;
  
  /* Initialising linked listed for backward simulation */
  struct Point *headDeleted, *headLiving, *dummyDeleted, *dummyLiving;
  struct Point *TempPoint;
  headLiving = (struct Point *) malloc(sizeof *headLiving);
  if(headLiving == NULL ) printf("error: malloc returned NULL!!\n");
  dummyLiving = (struct Point *) malloc(sizeof *dummyLiving);  
  if(dummyLiving == NULL ) printf("error: malloc returned NULL!!\n");
  headLiving->next = dummyLiving; dummyLiving->next = dummyLiving;

  headDeleted = (struct Point *) malloc(sizeof *headDeleted);
  if(headDeleted == NULL ) printf("error: malloc returned NULL!!\n");
  dummyDeleted = (struct Point *) malloc(sizeof *dummyDeleted);  
  if(dummyDeleted == NULL ) printf("error: malloc returned NULL!!\n");
  headDeleted->next = dummyDeleted; dummyDeleted->next = dummyDeleted;

  struct Point2 *TempCell, *TempCell2;

  struct Point3 *headTransition, *dummyTransition;
  headTransition = (struct Point3 *) malloc(sizeof *headTransition);
  if(headTransition == NULL) printf("error: malloc returned NULL\n");
  dummyTransition = (struct Point3 *) malloc(sizeof *dummyTransition);
  if(dummyTransition == NULL) printf("error: malloc returned NULL\n");
  headTransition->next = dummyTransition; 
  dummyTransition->next = dummyTransition;
  
  PP->GeneratePoisson(headLiving, &GeneratedPoints,
			      &LivingPoints,
			      &NoP);  
    
  StartTime=1;
  EndTime=1;

  TimeStep = 0; D0Time = 0;
  D0Living = GeneratedPoints;

  long int tmp, D0;
  
  do{
    tmp=BirthDeath(TimeStep, headLiving,
		      headDeleted,
		      headTransition);
    if(tmp>0){ 
      if(tmp>(LivingPoints+1-D0Living)){
	D0Living --;
      }
    }
    D0Time++;
  }while(D0Living>0);
  tmp=BirthDeath(TimeStep, headLiving,
		      headDeleted,
		      headTransition); 
  StartTime=1; EndTime=D0Time+1; D0 = 0;

  do{	 
    if(D0==1){
      for(TimeStep=StartTime;TimeStep<=EndTime;TimeStep ++){
	tmp=BirthDeath(TimeStep, headLiving,
		       headDeleted,
		       headTransition);      
      }
    }
    D0 = 1;
    P2P->Empty();
    
    /*
    headUpper->next = dummyUpper; dummyUpper->next = dummyUpper;
    for(XCell=0;XCell<=P2P->MaxXCell;XCell++){
      for(YCell=0;YCell<=P2P->MaxYCell;YCell++){
	headUpperCell[XCell][YCell]->next=dummyUpper;
      }
    }
    */
    
    P2P->UpperLiving[0] = LivingPoints;
    P2P->UpperLiving[1] = 0;
    
    P2P->NoP = 0;
    i=0;
    TempPoint = headLiving->next;
    while(TempPoint!=TempPoint->next){
      i++;
      TempCell2 = (struct Point2 *) malloc(sizeof *TempCell2);
      if(TempCell2 == NULL) printf("error: malloc returned NULL\n");
      TempCell2->No = TempPoint->No;
      TempCell2->X = TempPoint->X;
      TempCell2->Y = TempPoint->Y;
      TempCell2->InLower[0] = 1;
      TempCell2->InLower[1] = 0;
      f1 = (TempPoint->X-P2P->Xmin)/P2P->XCellDim;  XCell = int(floor(f1));
      if(XCell>P2P->MaxXCell) XCell = P2P->MaxXCell;
      f1 = (TempPoint->Y-P2P->Ymin)/P2P->YCellDim;  YCell = int(floor(f1));
      if(YCell>P2P->MaxYCell) YCell = P2P->MaxYCell;
      TempCell2->next = P2P->headCell[XCell][YCell]->next;
      P2P->headCell[XCell][YCell]->next = TempCell2;
      
      TempPoint = TempPoint->next;
    }
    
    //P2P->DumpToFile("temp0.dat");
    
    struct Point3 *TempTransition;
    struct Point *Proposal;
    
    TempTransition = headTransition->next;
    Proposal = headDeleted->next;
    DDD = GeneratedPoints;
    
    for(TimeStep=EndTime;TimeStep>=1;TimeStep--){
      R_CheckUserInterrupt();
      Forward(TimeStep,TempTransition->Case,
	      TempTransition->XCell,TempTransition->YCell,
	      Proposal,&DDD);
      if(TempTransition->Case == 1) Proposal = Proposal ->next;
      TempTransition = TempTransition->next;
    }
    
    /* Doubling strategy used!*/
    StartTime = EndTime+1;
    EndTime=EndTime*2;
    
    //P2P->DumpToFile("temp.dat");
    
  }while(P2P->UpperLiving[0]!=P2P->UpperLiving[1]);
  P2P->Clean();
  i=0;
  struct Point *TempPoint2;
  TempPoint = headLiving;
  TempPoint2 = headLiving->next;
  while(TempPoint!=TempPoint->next){
    i++;
    free(TempPoint);
    TempPoint = TempPoint2;
    TempPoint2 = TempPoint2->next;
  }
  free(TempPoint);
  
  i = 0;
  TempPoint = headDeleted;
  TempPoint2 = headDeleted->next;
  while(TempPoint!=TempPoint->next){
    i++;
    free(TempPoint);
    TempPoint = TempPoint2;
    TempPoint2 = TempPoint2->next;
  }
  free(TempPoint);
  //printf("%d ",i);

  struct Point3 *TempTransition,*TempTransition2;

  i = 0;
  TempTransition = headTransition;
  TempTransition2 = headTransition->next;
  while(TempTransition!=TempTransition->next){
    i++;
    free(TempTransition);
    TempTransition = TempTransition2;
    TempTransition2 = TempTransition2->next;
  }
  free(TempTransition);
  //printf("%d ST: %d ET: %d\n",i,StartTime,EndTime);
  //scanf("%f",&f1);

};

extern "C" {
  void PerfectStrauss(double *beta,
		      double *gamma,
		      double *r,
		      double *xmin,
		      double *xmax,
		      double *ymin,
		      double *ymax,
		      int *noutmax,
		      double *xout,
		      double *yout,
		      int *nout){
    int xcells, ycells;
    xcells = (int) fmin(floor((*xmax-*xmin)/ *r),9.0);
    ycells = (int) fmin(floor((*ymax-*ymin)/ *r),9.0);
    //printf("xcells %d   ycells %d\n",xcells,ycells);
    // Initalise point Strauss point process
    StraussProcess ExampleProcess(*xmin,*xmax,*ymin,*ymax,*beta,*gamma,*r);  
    // Initalise point pattern
    Point2Pattern ExamplePattern(*xmin,*xmax,*ymin,*ymax, xcells, ycells);
    // parameters: min x, max x, min y, max y, "cells" in x and y direction
    // used for speeding up neighbour counting, 9 is max here
    
    // Initalise perfect sampler
    Sampler PerfectSampler(&ExampleProcess);
    
    // Perform perfect samling
    PerfectSampler.Sim(&ExamplePattern);
    
    ExamplePattern.Return(xout, yout, nout, *noutmax);
    
  }
}
