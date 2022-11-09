/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2022 Victor Matveev
 *                    LBM/NIDDK/NIH and DMS/NJIT
 *
 *              Calcium Calculator.cpp / calc.cpp
 *
 *   The main() routine and the highest-order ADI engine procedures;
 *   tortuosity function interpreter/initializer
 *
 ************************************************************************
 
    This file is part of Calcium Calculator (CalC).

    CalC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CalC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CalC.  If not, see <https://www.gnu.org/licenses/>

 ************************************************************************/

#define _CRT_SECURE_NO_WARNINGS
#define GL_SILENCE_DEPRECATION

#include "PlatformSpecific.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "syntax.h"
#include "box.h"
#include "grid.h"
#include "vector.h"
#include "field.h"
#include "table.h"
#include "peak.h"
#include "markov.h"
#include "interpol.h"
#include "simulation.h"
#include "fplot.h"
#include "gate.h"
#include "loop.h"
#include "time.h"

double  CHARGE_LOSS;

int     DIMENSIONALITY = 3;
int     VERBOSE        = 3;
int     GEOMETRY       = CARTESIAN3D;
int     GEOMETRY1      = 0;
int     GEOMETRY2      = 0;
int     GEOMETRY3      = 0;
char    LABEL_DIM1[2];
char    LABEL_DIM2[6];
char    LABEL_DIM3[4];

//********************************************************************************************

PlotArray* GluPlotArray = NULL;
char*      globalLabelX = NULL;
char*      versionStr = StrCpy("6.10.1");
char       scriptFileName[1024];

#define EXTRA_PARAM_STRING  "; pA=5.182134 ; pi=4 atan(1) ; "

#define RETURN_NORMAL 0     // return types used by ADI step methods
#define RETURN_END    1

//********************************************************************************************
const char *getMethod(CaMethod &, BufMethod &);
double get_sim_time(TokenString &params);

void setFieldByFunction(TokenString &TS, long p, double *Diff, const char *errStr);
void setFieldByFunction(TokenString &TS, long p, bool   *Diff, const char *errStr);

int Buf1Dstep(BufferObj &Buf, VectorObj &Buf_next, FieldObj &Ca, VectorObj &, double dt, double);
int Ca1Dstep(FieldObj &Ca, VectorObj &Canext, BufferArray &, VectorObj &, double dt, double);

int Buf2Dstep(BufferObj &Buf, VectorObj &Buf_next, FieldObj &Ca, VectorObj &, double dt, double);
int Ca2Dstep(FieldObj &Ca, VectorObj &Canext, BufferArray &, VectorObj &, double dt, double);

int Ca3Dstep(FieldObj &Ca, VectorObj &Canext, BufferArray &, VectorObj &, double dt, double);
int Buf3Dstep(BufferObj &Buf, VectorObj &Buf_next, FieldObj &Ca, VectorObj &, double dt, double);

void getRun( TokenString &TS, int i, bool *adaptive, double *time=0, 
             double *accuracy=0, double *dtMax=0, double *dt=0, double *dtStretch=0, double *ODEaccuracy=0);


//**************************************************************************************************

void header() {
     fprintf(stderr,"\n******************************************************************");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*  Calcium Calculator (CalC)  *  version 6.10.1  *  Oct 10, 2022 *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*                Victor Matveev (C) 2001-2022                    *");
	 fprintf(stderr,"\n*   CalC is distributed under GPLv3: see attached license file   *");
     fprintf(stderr,"\n*                                                                *");
	 fprintf(stderr,"\n*  Dept of Math Sciences, New Jersey Institute of Technology     *");
     fprintf(stderr,"\n*                     and LBM, NIDDK, NIH                        *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*  Supported in part by NSF DMS0417416, DMS0817703, DMS517085    *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n****************************************************************\n");
}

 //**************************************************************************************************
 //                                           M A I N
 //**************************************************************************************************


 int main(int argc, char **argv) {


 char fname[1024];

 try {

   if (argc < 2) {
	 header();
     fprintf(stderr, "\n\n Enter the CalC script file name: ");
	 fflush(stderr);
     scanf("%s", fname);
	 unsigned long i = strlen(argv[0]);
	 while ( argv[0][i] != '/' && argv[0][i] != '\\' && i > 0 ) i--;
     strncpy(scriptFileName, argv[0], i+1);
	 scriptFileName[i+1] = 0;
	 fprintf(stderr, "\n Full path = %s \n", scriptFileName);
	 strcat (scriptFileName, fname);
     fprintf(stderr, " File name = %s \n\n", scriptFileName );
   }
   else strcpy(scriptFileName, argv[1]);

   TokenString *TS;
   TS = new TokenString(scriptFileName, EXTRA_PARAM_STRING, "", argc, argv);
   TS->get_int_param("verbose", &VERBOSE);
#ifndef _NO_GLUT_
   if ( TS->Assert("plot.method", "gl") )  glutInit(&argc, argv);
#endif
   if (VERBOSE) header();

   VectorObj result( getTrackVarNum(*TS) );
   LoopObj   vary( *TS, &result );
   long seed = long(time(NULL));
   seed =  TS->get_long_param("seed");
   srand(seed);
   delete TS;

   // ++++++++++++++++++++++++++++ MAIN LOOP ++++++++++++++++++++++++++++

   for (int step = 0; step < vary.steps; step++) {
     vary.step();
     TS = new TokenString(scriptFileName, EXTRA_PARAM_STRING, vary.loopVarString, argc, argv); 
     if (TS->token_count("exit"))     { if (VERBOSE) fprintf(stderr, "\n > Exit command: breaking execution <\n ");     delete TS; break;    }
     if (TS->token_count("continue")) { if (VERBOSE) fprintf(stderr, "\n > Continue command: breaking execution <\n "); delete TS; continue; }
     TS->get_int_param("verbose", &VERBOSE);
     SimulationObj *Simulation;
     double *trackPtr;

     if ( !TS->token_count("Run") ) {                  //######## CALCULATOR MODE #########
       if (VERBOSE) fprintf(stderr, "***** No simulation runs specified; running in calculator mode *****\n");
       if ( vary() ) 
		   for (int i = 0; i < result.size; i++) result[i] = TS->get_double( getTrackVar(*TS, i) );
       TS->printResults(0);
     }        
     else {
       if ( !(TS->token_count("Ca.D") || TS->token_count("grid")) || equal(TS->get_string("mode"), "ODE") ) {
            if (VERBOSE) fprintf(stderr, "***** Running in ODE-only mode *****\n");
	        Simulation = new ODESimulationObj(*TS);    //########  ODE solver mode: ########
       } 
       else Simulation = new SimulationObj(*TS);       //#########  FULL PDE MODE  #########
	   Simulation->Run();                              //#### RUN THE DIFFERENCE SCHEME ####
       if ( vary() ) {
		   for (int i = 0; i < result.size; i++) 
				if ( (trackPtr = Simulation->ResolveID( vary.trackIDs[i] )) ) result[i] = *trackPtr;
				else  TS->errorMessage( TS->token_index(vary.trackIDs[i]), 0, "Cannot track an undefined variable" );
	   }
       TS->printResults(Simulation);
#ifndef _NO_GLUT_
       if (Simulation->Plots->gl_on && Simulation->Plots->plot_num)  glutMainLoop();
#endif
       delete Simulation;
	 }

     vary.draw(); 
#ifndef _NO_GLUT_
     if (vary.plots->gl_on && step == vary.steps - 1) {
         GluPlotArray = vary.plots;
         glutMainLoop();
     }
#endif
     delete TS;
   } // ****************************** Loop over steps

 }
 catch (char *str) {  if ( !equal(str,"") ) {
                        fprintf(stderr, "\n\n*** Error: %s\n", str);
						fflush(stderr); 
                        delete [] str;  }
                   }
 catch (int i) {		fprintf(stderr, "\n\n*** CalC: breaking execution on exception #%d ***\n", i); 
                        perror(" System error (if any): ");
						fflush(stderr); }

 if (VERBOSE > 3)  {    fprintf(stderr, "\n\n*** Enter any string to exit CalC ***\n"); 
                        fflush(stderr);
		      			char s[100] = "0";
                        scanf("%s", s); 
                        }
 return 0;
 }


//**********************************************************************************************
 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 //**********************************************************************************************

//**************************************************************************
//                C A L C I U M   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) Ca*  = (1 + Ax/2 + Ay + Az) Ca + dt H(Ca,B)
//  (1 - Ay/2) Ca** = -Ay/2 Ca + Ca*
//  (1 - Az/2) Ca^  = -Az/2 Ca + Ca** + dt / 2 * { H(Ca^,B) - H(Ca,B) }
//
//  H(Ca,B) = Sum_B{ - kplus * Ca * B + kminus * (B_total - B) }
//
//  Ca = Ca(t_n), Ca^ = Ca(t_{n+1}), B = B(t_{n+1/2})
//**************************************************************************

int Ca3Dstep(FieldObj &Ca, VectorObj &Next, BufferArray &Buf, VectorObj &Lin, double dt, double T)
{
  double    t0 = Ca.Time;
  int       rvalue = RETURN_NORMAL;
 
  if ( (T > 0.0) && (t0 + dt >= T))  {  
	  rvalue = RETURN_END; dt = T - t0; 
  }

  double nu = Ca.getD() * dt;

  Next = Ca + dt * Ca.bgr * *Ca.kuptake;  

  Lin  = - dt *  *Ca.kuptake;
  for (int b = 0; b < Buf.buf_num; b++)
    {
    Lin  -= (dt * Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];
    Next += ( dt * Buf.array[b]->kminus->Evaluate()) * (*(Buf.array[b]->total) - *Buf.array[b]);
    }

  Next += Lin % Ca;       
  Ca.Time += 0.5 * dt;
  Ca.add_sources(Next, dt);  
 
	  Ca.Run3Dz(-0.5 * nu,        nu,   nu,  0.5*nu,  0,  Next.elem);  
	  Ca.Run3Dx(-0.5 * nu, -0.5 * nu,   0.,      0.,  0,  Next.elem);  
	  Lin *= 0.5;
	  Next -= (Lin % Ca);
	  Ca.Run3Dy(-0.5 * nu,   0., -0.5 * nu, 0., Lin.elem, Next.elem);  


  if (T > 0.0) { Ca = Next; Ca.Time += 0.5 * dt; } /*Ca.cleanBoundaries();*/ 
  else Ca.Time = t0;

  return rvalue;
}

//**************************************************************************
//                  B U F F E R   D - G   S T E P
//**************************************************************************

int Buf3Dstep(BufferObj &Buf, VectorObj &Next, FieldObj &Ca, VectorObj &Lin, double dt, double T)
{ 
int    rvalue = RETURN_NORMAL;

  if ( (T > 0.0) && (Buf.Time + dt >= T)) 
    { dt = T - Buf.Time; /*printf("Buf: dt reduced to %g\n", dt);*/ rvalue = RETURN_END; }
  double nu = Buf.getD() * dt;
  
  double tstore = Buf.Time;
  Buf.Time += dt*0.5;
  double kplus  = Buf.kplus->Evaluate();
  double kminus = Buf.kminus->Evaluate();
  Buf.Time = tstore;

  Lin = -dt * (kplus * Ca + kminus);
  Next = Buf % ( 1.0 + Lin ) + dt * kminus * *(Buf.total);

	  Buf.Run3Dz(-0.5 * nu, nu, nu, 0.5 * nu, 0, Next.elem);
	  Buf.Run3Dx(-0.5 * nu, -0.5 * nu, 0., 0., 0, Next.elem);
	  Lin *= 0.5;
	  Next -= (Lin % Buf);
	  Buf.Run3Dy(-0.5 * nu, 0., -0.5 * nu, 0., Lin.elem, Next.elem);

 if (T > 0.0) {  Buf.Time += dt; Buf = Next; } /*Buf.cleanBoundaries();*/ 
  
  return rvalue;
}

//**************************************************************************
//                C A L C I U M   2 D   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) Ca*  = (1 + Ay/2) Ca  + dt/2 H(Ca,B)
//  (1 - Ay/2) Ca^  = (1 + Ax/2) Ca* + dt/2 H(Ca^,B)
//
//  H(Ca,B) = Sum_B{ - kplus * Ca * B + kminus * (B_total - B) }
//
//  Ca = Ca(t_n), Ca^ = Ca(t_{n+1}), B = B(t_{n+1/2})
//**************************************************************************

int Ca2Dstep(FieldObj &Ca, VectorObj &Next, BufferArray &Buf, VectorObj &Lin, double dt, double T)
{
  double    t0 = Ca.Time;
  int       b, rvalue = RETURN_NORMAL;
  FieldObj  Temp(Ca);

  if ( (T > 0.0) && (t0 + dt >= T))  { rvalue = RETURN_END; dt = T - t0; }

  double nu = Ca.getD() * dt;

  Next = Ca + 0.5 * dt * Ca.bgr * *Ca.kuptake;
  Lin = - 0.5 * dt *  *Ca.kuptake; 

  for (b = 0; b < Buf.buf_num; b++)
    {
    Lin  -= (0.5 * dt * Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];
    Next += (*(Buf.array[b]->total) - *Buf.array[b]) * (0.5 * dt * Buf.array[b]->kminus->Evaluate());
    }

  Next += Lin % Ca;
  Ca.Time += 0.5 * dt;
  Ca.add_sources(Next, 0.5 * dt);

  Ca.Run2Dx(    -0.5 * nu,     0.0,    0.5 * nu,  0,  Next.elem);

  Temp = Next;  // IMPORTANT! cf. 3D methods: Ca* in second split operator, not Ca
  Next += 0.5 * dt * Ca.bgr * *Ca.kuptake; 
  for (b = 0; b < Buf.buf_num; b++)    
    Next += (*(Buf.array[b]->total) - *Buf.array[b]) * (0.5 * dt * Buf.array[b]->kminus->Evaluate());
  Ca.add_sources(Next, 0.5 * dt);

  Temp.Run2Dy(-0.5 * nu, 0.5 * nu,     0.0,     Lin.elem,  Next.elem);

  if (T > 0.0) { Ca = Next; Ca.Time += 0.5 * dt; } /*Ca.cleanBoundaries();*/ 
  else Ca.Time = t0;
  return rvalue;
}

//**************************************************************************
//               B U F F E R   2 D   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) B*  = (1 + Ay/2) B  + dt/2 H(Ca,B)
//  (1 - Ay/2) B^  = (1 + Ax/2) B* + dt/2 H(Ca,B^)
//
//  H(Ca,B) = - kplus * Ca * B + kminus * (B_total - B)
//
//  B = B(t_n), B^ = B(t_{n+1}), Ca = Ca(t_{n+1/2})
//**************************************************************************

int Buf2Dstep(BufferObj &Buf, VectorObj &Next, FieldObj &Ca, VectorObj &Lin, double dt, double T)
{ 
int       rvalue = RETURN_NORMAL;
FieldObj  Temp(Buf);

  if ( (T > 0.0) && (Buf.Time + dt >= T)) 
    { dt = T - Buf.Time; /*printf("Buf: dt reduced to %g\n", dt);*/ rvalue = RETURN_END; }
  double nu = Buf.getD() * dt;
  
  double tstore = Buf.Time;
  Buf.Time += dt*0.5;
  double kplus  = Buf.kplus->Evaluate();
  double kminus = Buf.kminus->Evaluate();
  Buf.Time = tstore;

  Lin = -0.5 * dt * (kplus * Ca + kminus);

  Next = Buf + 0.5 * dt * kminus * *(Buf.total) + Lin % Buf;

  Buf.Run2Dx(    -0.5 * nu,     0.0,    0.5 * nu,  0,  Next.elem);

  Temp = Next; // IMPORTANT! cf. 3D methods.
  Next += 0.5 * dt * kminus * *Buf.total;
  
  Temp.Run2Dy(-0.5 * nu, 0.5 * nu,     0.0,     Lin.elem,  Next.elem);

  if (T > 0.0) { Buf.Time += dt; Buf = Next; } /*Buf.cleanBoundaries();*/ 
  return rvalue;
}


/**************************************************************************
     C A L C I U M   C R A N K - N I C H O L S O N   1 D   S T E P
***************************************************************************

  (1 - Ax/2) Ca^ = (1 + Ax/2) Ca + dt/2 ( H(Ca, B) + H(Ca^,B) ) + dt sources(t_{n+1/2})
                 = (1 + Ax/2) Ca - dt/2 Sum_B{ kplus * B } ( Ca + Ca^ )
                 + dt { Sum_B{ kminus * (B_total - B) + sources(t_{n+1/2}) }

  ( 1 - Ax/2 - Lin ) Ca^ = Ax/2 Ca + (1 + Lin ) Ca
                         + dt { Sum_B{ kminus * (B_total - B) + sources(t_{n+1/2}) }

  Lin = - dt/2 Sum_B{ kplus * B }


  H(Ca,B) = Sum_B{ - kplus * Ca * B + kminus * (B_total - B) }

  Ca = Ca(t_n), Ca^ = Ca(t_{n+1}), B = B(t_{n+1/2})
****************************************************************************/

int Ca1Dstep(FieldObj &Ca, VectorObj &Next, BufferArray &Buf, VectorObj &Lin, double dt, double T)
{
  int       rvalue = RETURN_NORMAL;
  double    t0 = Ca.Time;
 
  if ( (T > 0.0) && (t0 + dt >= T))  { rvalue = RETURN_END; dt = T - t0; }

  double nu = Ca.getD() * dt;

  Next = Ca  + dt * Ca.bgr * *Ca.kuptake;
  Lin = -0.5 * dt * *Ca.kuptake; 

  for (int b = 0; b < Buf.buf_num; b++) {
    Lin  -= (0.5 * dt * Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];
    Next += ( dt * Buf.array[b]->kminus->Evaluate() ) * (*(Buf.array[b]->total) - *Buf.array[b]);
    }

  Next += Ca % Lin;
  Ca.Time += 0.5 * dt;
  Ca.add_sources(Next, dt);
  Ca.Run1D( -0.5 * nu, 0.5 * nu, Lin.elem, Next.elem);

  if (T > 0.0) { Ca = Next; Ca.Time += 0.5 * dt; /*Ca.cleanBoundaries();*/ }
    else Ca.Time = t0;

  return rvalue;
}

/**************************************************************************
       B U F F E R    C R A N K - N I C H O L S O N   1 D   S T E P
**************************************************************************

  (1 - Ax/2) B^ = (1 + Ax/2) B + dt/2 ( H(Ca, B) + H(Ca,B^) )
                = (1 + Ax/2) B - dt/2 (kplus Ca + kminus) ( B + B^ ) + dt kminus B_total 

  ( 1 - Ax/2 - Lin ) B^ = Ax/2 Ca + (1 + Lin ) B + dt kminus B_total

  Lin = - dt/2 (kplus Ca + kminus)

  H(Ca,B) = - kplus * Ca * B + kminus * (B_total - B)

  B = B(t_n), B^ = B(t_{n+1}), Ca = Ca(t_{n+1/2})
****************************************************************************/

int Buf1Dstep(BufferObj &Buf, VectorObj &Next, FieldObj &Ca, VectorObj &Lin, double dt, double T)
{
  int    rvalue = RETURN_NORMAL;

  if ( (T > 0.0) && (Buf.Time + dt >= T)) { dt = T - Buf.Time; rvalue = RETURN_END; }
  double nu = Buf.getD() * dt;

  double tstore = Buf.Time;
  Buf.Time += dt*0.5;
  double kplus  = Buf.kplus->Evaluate();
  double kminus = Buf.kminus->Evaluate();
  Buf.Time = tstore;

  Lin  = -0.5 * dt * (kplus * Ca + kminus );
  Next = Buf % ( 1.0 + Lin ) + (dt * kminus) * *Buf.total;

  Buf.Run1D( -0.5 * nu, 0.5 * nu, Lin.elem, Next.elem);

  if (T > 0.0) { Buf = Next; Buf.Time += dt; /*Buf.cleanBoundaries();*/ }

  return rvalue;
}


//**************************************************************************

const char *getMethod( CaMethod &CaStep, BufMethod &BufStep) {

 switch (DIMENSIONALITY)  { 

 case 1: CaStep = &Ca1Dstep; BufStep = &Buf1Dstep; 
	     return makeMessage("1D (%s) Crank-Nicholson implicit", LABEL_DIM1);

 case 2: CaStep = &Ca2Dstep; BufStep = &Buf2Dstep; 
	     return  makeMessage("2D (%s,%s) Crank-Nicholson ADI", LABEL_DIM1, LABEL_DIM2);

 case 3: CaStep = &Ca3Dstep; BufStep = &Buf3Dstep;       
	     return  makeMessage("3D (%s,%s,%s) Douglas-Gunn ADI", LABEL_DIM1, LABEL_DIM2, LABEL_DIM3);
		 
 default: throw  makeMessage("Can't interpret DIMENSIONALITY=%d", DIMENSIONALITY);
		  return NULL;
 }
 }

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


void getRun( TokenString &TS, int i, bool *adaptive, double *time, double *accuracy,
			 double *dtMax, double *dt, double *dtStretch, double *ODEaccuracy)
{
  char firstArg[MAX_TOKEN_LENGTH];

  if (adaptive) *adaptive = false;             // non-adaptive is the default
  long pos = TS.token_index("Run", i);

  TS.trail_pars(pos, 's', firstArg, 'E');
  pos ++;
  int pmax = TS.tokens_to_eol(pos);

  if ( pmax > 0 && !TS.isConst(firstArg) ) { 
    pos++; 
    if (adaptive) {
		if ( equal(firstArg,"adaptive") ) *adaptive = true;
	    else if ( !equal(firstArg, "nonadaptive") ) TS.errorMessage(pos, 0, "Undefined run method");
	}
  }

  if (time) {
    *time = TS.get_double(pos);
    if (*time <= 0) TS.errorMessage(pos, makeMessage("Simulation time has to be positive (T = %lf)", *time) );
  } 

  if ( ! isLineEnd(TS[pos + 1]) ) {
	if (*adaptive) 
		TS.trail_pars(pos, 'd', accuracy, 'd', dtMax, 'd', dt, 'd', dtStretch, 'd', ODEaccuracy, 'E');
	else 
		TS.trail_pars(pos, 'd', dtMax, 'd', accuracy, 'E');
  }
  if (dt) if (*dt <= 0) 
    TS.errorMessage(pos+1, makeMessage("Time step must be positive (dt = %lf)", *dt) );
  if (accuracy) if (*accuracy <= 0 || *accuracy > 1.0) // this argument is actually ODEaccuracy for non-adaptive run
    TS.errorMessage(pos+2, makeMessage("Accuracy must be between 0 and 1 (accuracy = %lf)", *accuracy) );
  if (dtMax) if (*dtMax <= 0 ) 
    TS.errorMessage(pos+1, makeMessage("Max time step must be positive (dtMax = %lf)", *dtMax) );
  if (dtStretch) if (*dtStretch <= 1 || *dtStretch >= 2) 
    TS.errorMessage(pos+4, makeMessage("Time step stretch must be between 1 and 2 (dtStretch = %lf)", *dtStretch) );
  if (ODEaccuracy) if (*ODEaccuracy <= 0 || *ODEaccuracy > 1.0) 
    TS.errorMessage(pos+5, makeMessage("ODE accuracy must be between 0 and 1 (ODEaccuracy = %lf)", *ODEaccuracy) );
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//**************************************************************************

double get_sim_time(TokenString &params)
  {
  double tt, t = 0.0;
  bool adaptive;

  for (int i = 1; i<= params.token_count("Run"); i++) {
    getRun(params, i, &adaptive, &tt); 
    t += tt;
  }
  return t;
  }

//**************************************************************************


void setFieldByFunction(TokenString &TS, long p, double *Field, const char *errStr) {

  int xsize = FieldObj::Grid->xsize;  double *xcoord = FieldObj::Grid->xcoord;
  int ysize = FieldObj::Grid->ysize;  double *ycoord = FieldObj::Grid->ycoord;
  int zsize = FieldObj::Grid->zsize;  double *zcoord = FieldObj::Grid->zcoord;

  int  ix, iy, iz;
  long l = 0;
  double x, y, z;
  ExpressionObj *T = 0;

  switch (DIMENSIONALITY) {

  case 1:
    T = new ExpressionObj(TS, p, makeMessage("%s (function of %s)", errStr, LABEL_DIM1), 0, 0, &x, LABEL_DIM1); 
    for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate();
    }
    break;  

  case 2:
    T = new ExpressionObj(TS, p, 
            makeMessage("%s (function of %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2), 0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2); 
    for (iy = 0; iy < ysize; iy++) {
        y = ycoord[iy];
        for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate();
	}
    }
    break;  
 
  case 3:
    T = new ExpressionObj(TS, p, 
         makeMessage("%s (function of %s, %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2, LABEL_DIM3), 
         0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2, &z, LABEL_DIM3); 
    for (iz = 0; iz < zsize; iz++) {
      z = zcoord[iz];
      for (iy = 0; iy < ysize; iy++) {
        y = ycoord[iy];
        for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate();
	}
      }
    }
    break;  
 }

 if (T && VERBOSE) T->print(stderr); 
 if (T) delete T;
 return;
}

//**************************************************************************


void setFieldByFunction(TokenString &TS, long p, bool *Field, const char *errStr) {

  int xsize = FieldObj::Grid->xsize;  double *xcoord = FieldObj::Grid->xcoord;
  int ysize = FieldObj::Grid->ysize;  double *ycoord = FieldObj::Grid->ycoord;
  int zsize = FieldObj::Grid->zsize;  double *zcoord = FieldObj::Grid->zcoord;

  int  ix, iy, iz;
  long l = 0;
  double x, y, z;
  ExpressionObj *T = 0;

  switch (DIMENSIONALITY) {

  case 1:
    T = new ExpressionObj(TS, p, makeMessage("%s (function of %s)", errStr, LABEL_DIM1), 0, 0, &x, LABEL_DIM1); 
    for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate() > 0;
    }
    break;  

  case 2:
    T = new ExpressionObj(TS, p, 
            makeMessage("%s (function of %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2), 0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2); 
    for (iy = 0; iy < ysize; iy++) {
        y = ycoord[iy];
        for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate() > 0;
	}
    }
    break;  
 
  case 3:
    T = new ExpressionObj(TS, p, 
         makeMessage("%s (function of %s, %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2, LABEL_DIM3), 
         0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2, &z, LABEL_DIM3); 

	for (iz = 0; iz < zsize; iz++) 
	{
		z = zcoord[iz]; 
		for (iy = 0; iy < ysize; iy++) 
		{
			y = ycoord[iy]; 
			for (ix = 0; ix < xsize; ix++)  
			{
				x = xcoord[ix];
				double val = T->Evaluate();
				Field[l++] = (val > 0);
			}
		}
	}
	break;  
 }

  if (T && VERBOSE) { T->print(stderr); fprintf(stderr, "\n"); }
 if (T) delete T;
 return;
}

//**************************************************************************
//**************************************************************************
