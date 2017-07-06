#pragma once

#include "lib.h"


class CTagMatch
{
private:	
	
public:
	CTagMatch(){

	}
	~CTagMatch(){

	}	

int match(CPolyTag &pt1,CPolyTag &pt2,double error=0.2){
	m1.clear();
	m2.clear();
	double t;
	for(int i=0;i<pt1.ang.size();i++){
		for(int j=0;j<pt2.ang.size();j++){
			t=fabs(pt1.ang[i]-pt2.ang[j]);
			if(t>error&&t<=PI) continue;
			if(t<2*PI-error&&t>PI) continue;
			if(fabs(pt1.range[i]-pt2.range[j])>error) continue;
			m1.push_back(i);
			m2.push_back(j);
			break;
		}
	}
	
	return m1.size();
}	

int ref(CPolyTag &pt,int x,double error=PI/32){
	double m=-0.75*PI+x*PI/2;
		double t;
	for(int i=0;i<m1.size();i++){
			t=fabs(m-pt.ang[i]);
			if(t>error&&t<=PI) continue;
			if(t<2*PI-error&&t>PI) continue;
			cout<<pt.ang[i]<<endl;
			for(int j=0;j<m1.size();j++){
			cout<<pt.ang[j]<<" |||   ";
			}
			cout<<endl;

			return i;
	}
	return -1;
}
	
private:

public:
vector <int> m1;
vector <int> m2;


};


