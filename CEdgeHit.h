#pragma once

#include "lib.h"


class CEdgeHit
{
public:
	CEdgeHit(){
	prev=NULL;
	next=NULL;
	used=0;
	}
int x;
int y;
int used;
CEdgeHit *prev;
CEdgeHit *next;

void print(){
cout<<x<<" "<<y<<" "<<used<<endl;
}

};

