/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_
#include "./basic.h"
#include "./debug.h"

class Utils {
 public:
  static char RandomChar();
  // Genetic:
  static void CreateTrio(char * sequence,
                         size_t seq_len,
                         double error_ratio,
                         double short_indels,
                         double long_indels,
                         double repeat_ratio,
                         double str_ratio,
                         double parents_pm_ratio,
                         double child_pm_ratio,
                         char ** _M1, char ** _M2, size_t * _M_len,
                         char ** _F1, char ** _F2, size_t * _F_len,
                         char ** _C1, char ** _C2, size_t * _C_len);
  static void CreateGuy(char * seq,
                        size_t seq_len,
                        double pm_ratio,
                        char **g1,
                        char **g2,
                        int * n_dels);
  static void CreateOffspring(char * M1, char * M2, size_t M_len,
                              char * F1, char * F2, size_t F_len,
                              double pm_ratio,
                              char **C1, char **C2, size_t *C_len,
                              int  * non_inherited_dels, int * non_inherited_snps);
  static int MutateSTR(char** A ,
                       char** B,
                       char** C,
                       char** D,
                       size_t *len,
                       double mutation_ratio);
  static int MutateREPEATS(char** A ,
                           char** B,
                           char** C,
                           char** D,
                           size_t *len,
                           double mutation_ratio);
  static int MutateLONGINS(char** A ,
                           char** B,
                           char** C,
                           char** D,
                           size_t *len,
                           double mutation_ratio);
  static int MutateSHORTINS(char** A ,
                            char** B,
                            char** C,
                            char** D,
                            size_t *len,
                            double mutation_ratio);
  static int MutateGENERIC_INSERT(char **A ,
                                  char **B,
                                  char **C,
                                  char **D,
                                  size_t * len,
                                  double mutation_ratio,
                                  size_t min_size,
                                  size_t max_size,
                                  bool do_repeat,
                                  bool do_tandem);
  static int MutateLONGDEL(char *A , size_t len, double mutation_ratio);
  static int MutateSHORTDEL(char *A , size_t len, double mutation_ratio);
  static int MutateGENERIC_DEL(char *A ,
                               size_t len,
                               double mutation_ratio,
                               size_t min_size,
                               size_t max_size);
  static int MutateSNP(char *A , size_t len, double mutation_ratio);
  static int MutatePointDEL(char *A , size_t len, double mutation_ratio);
  static int MutateNoise(char *A , size_t len, double ratio);
  static char MutateChar(char a);
  /////
  static int PhaseErrors(char * C1, char* C2, char* phase1, char* phase2, size_t len);
  static int CountNonGaps(char* seq, size_t len);
  static char * CopySeq(char * seq, size_t len);
  static void Truncate(char ** seq, size_t new_len);
  static void Compact(char ** A,
                      char ** B,
                      size_t *len);
  static void Scramble(char *A , char *B, size_t len, size_t chunk_len);
  static char * MeioticShuffle(char *A , char *B, size_t len, size_t n_cuts);
  static size_t * GenerateCuts(size_t len, size_t n_cuts);
  // General:
  static void StartClock();
  static double StopClock();
  static void ReadFastaFile(char * file_name, char **ans, size_t * len);
  static void SaveFastaFile(char * out_file_name, char *seq, size_t  len);
  static void SaveChar(char * seq, size_t length, char * file_name);
};
#endif  // SRC_UTILS_H_
