/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */

#ifndef SRC_DEBUG_H_
#define SRC_DEBUG_H_

#include <stdarg.h>
#include <iostream>
#include "./basic.h"

class Debug {
 public:
  template<typename TYPE>
static void PrintArray(TYPE *A, size_t _length) {
  if (_length == 0) {
    std::cout << "[]";
    return;
  }
  std::cout << "[";
  for (size_t i = 0; i < _length-1; i++) {
    // printf("%lu|", A[i]);
    std::cout << A[i] << "|";
  }
  std::cout << A[_length-1] << "]" << std::endl;
}

static void PrintLine(size_t _length) {
  std::cout << "-";
  for (size_t i = 0; i < _length; i++) {
    std::cout << "--";
  }
  std::cout << "-" << std::endl;
}

template<typename TYPE>
static bool equalArrays(TYPE *A, TYPE  * B, size_t _length) {
  for (size_t i = 0; i < _length; i++) {
    if (A[i] != B[i]) {
      PrintArray(A, _length);
      printf("\n");
      PrintArray(B, _length);
      printf("Diff in pos %u \n", (uint)i);
      return false;
    }
  }
  return true;
}

static void AbortPrint(const char* format, ...) {
  va_list args;
  char * buffer = new char[1024];
  va_start(args, format);
  vsprintf(buffer, format, args);
  va_end(args);

  fprintf(stderr, "\n***********\nABORTING::  %s ::\n***********\n", buffer);
  fflush(stderr);
  delete[]buffer;
  exit(-1);
}
};
#endif  // SRC_DEBUG_H_
