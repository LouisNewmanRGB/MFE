%module SimulationCpp

%{
#define SWIG_FILE_WITH_INIT
#include "Simulation.h"
%}

%include "windows.i"
%include "numpy.i"

%init %{
import_array();
%}

%include "Export.h"
//%apply (double* IN_ARRAY1, int DIM1) {(double* vector1, int length1), (double* input_vector1, int input_length1), (double* input_vector2, int input_length2)};
%apply (double* IN_ARRAY1, int DIM1) {(double a[3], int len_a), (double b[3], int len_b), (double c[3], int len_c), (double* sequenceTimes, int len_sequenceTimes)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double *sequenceVectors, int rows_sequenceVectors, int cols_sequenceVectors)};
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** output_vector1, int* output_length1)};
%include "Simulation.h"