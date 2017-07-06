#pragma once

#include "lib.h"
#include <windows.h>

class CTime
{
public:
	CTime(){
    setup();
	}
void start(){
     QueryPerformanceCounter(&nBeginTime); 

}	

void end(){
     QueryPerformanceCounter(&nEndTime);
}

void setup(){
 QueryPerformanceFrequency(&nFreq);
}

void print(char *token){
	if(debug==0) return;
     time=(double)(nEndTime.QuadPart-nBeginTime.QuadPart)/(double)nFreq.QuadPart;
     printf("%s %f\n",token,time);
}

     LARGE_INTEGER nFreq;

     LARGE_INTEGER nBeginTime;

     LARGE_INTEGER nEndTime;

     double time;
	 int debug;

};

