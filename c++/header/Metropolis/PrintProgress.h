#pragma once

#include <iostream>

	static void print_progress (double percentage) ; 
	void print_progress (double percentage) {
	
		#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
		#define PBWIDTH 60

    	int val = (int) (percentage * 100);
    	int lpad = (int) (percentage * PBWIDTH);
    	int rpad = PBWIDTH - lpad;
    	printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    	fflush (stdout);
	
		#undef PBSTR
		#undef PBWIDTH
	}

