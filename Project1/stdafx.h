// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <conio.h>
#include <tchar.h>


// TODO: reference additional headers your program requires here

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include <time.h>
#include <omp.h>
//#include <math.h>
#include <cmath>

#include "mkl.h"

#include <locale.h>
#include <Windows.h>
#include <algorithm>
#include <iomanip>


int user_pardiso(const vector<int>&,    const vector<int>& _jg,
                 const vector<double>&, const vector<double>& _al,
                 const vector<double>&, const vector<double>& _b,
                 vector<double>& _solution);