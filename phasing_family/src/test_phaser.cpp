/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */


#include <sys/times.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cassert>
#include "./phaser.h"
#include "./utils.h"
#include "./debug.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

score_t SCORE_GAP = -6;
score_t SCORE_MISMATCH = -10;
score_t SCORE_MATCH = 100;

/*
score_t SCORE_GAP = -2;
score_t SCORE_MISMATCH = -4;
score_t SCORE_MATCH = 5;
*/

bool global_fail = false;
int global_success_count = 0;
void Fail();
void Success();
void Summary();

void TestPhaserShortK_A();
void TestPhaserShortK_B();
void TestPhaserShortK_C();
void TestPhaserShortK_D();

void TestPhaserSimple();

void TestPhaserPartialDeleteChild_A();
void TestPhaserPartialDeleteChild_B();
void TestPhaserPartialDeleteChild_C();

void TestPhaserIdentical_A();
void TestPhaserIdentical_B();
void TestPhaserIdentical_C();
void TestPhaserIdentical_D();
void TestPhaserIdentical_E();
void TestPhaserIdentical_F();
void TestPhaserIdentical_G();
void TestPhaserIdentical_H();

void TestPhaserOnlyMotherRecomb_A();
void TestPhaserOnlyMotherRecomb_B();
void TestPhaserOnlyMotherRecomb_C();
void TestPhaserOnlyMotherRecomb_D();
void TestPhaserOnlyMotherRecomb_E();
void TestPhaserOnlyMotherRecomb_F();
void TestPhaserOnlyMotherRecomb_G();
void TestPhaserOnlyMotherRecomb_H();

void TestPhaserFatherMotherRecomb_A();
void TestPhaserFatherMotherRecomb_B();
void TestPhaserFatherMotherRecomb_C();
void TestPhaserFatherMotherRecomb_D();
void TestPhaserFatherMotherRecomb_E();
void TestPhaserFatherMotherRecomb_F();
void TestPhaserFatherMotherRecomb_G();
void TestPhaserFatherMotherRecomb_H();

void TestPhaserAllRecomb_A();
void TestPhaserAllRecomb_B();
void TestPhaserAllRecomb_C();
void TestPhaserAllRecomb_D();
void TestPhaserAllRecomb_E();
void TestPhaserAllRecomb_F();
void TestPhaserAllRecomb_G();
void TestPhaserAllRecomb_H();

void TestPhaserMismatchIndel();
void TestPhaserParentsGapped_A();
void TestPhaserParentsGapped_B();
void TestPhaserParentsGapped_C();
void TestPhaserParentsGapped_D();

void TestPhaserChildGapped_A();
void TestPhaserChildGapped_B();
void TestPhaserChildGapped_C();

void TestPhaserExhaustiveVarLength();
void TestPhaserExhaustiveSameLength();

void TestFasta();


void Fail() {
  printf(ANSI_COLOR_RED "\tTest Failed\n\n" ANSI_COLOR_RESET);
  global_fail = true;
}
void Success() {
  global_success_count++;
  printf(ANSI_COLOR_BLUE "\t\tSuccess\n" ANSI_COLOR_RESET);
}
bool equalPhases(char * p1, char * p2, size_t len);
bool equalPhases(char * p1, char * p2, size_t len) {
  bool strict = true;
  for (size_t i = 0; i < len; i++) {
    if (p1[i] != p2[i]) {
      if (strict || (p1[i] != '?' && p2[i] != '?')) {
        Debug::PrintArray(p1, len);
        Debug::PrintArray(p2, len);
        return false;
      }
    }
  }
  return true;
}

void Summary() {
  if (!global_fail) {
    printf(ANSI_COLOR_GREEN "\n\t%i TEST RUN\n" ANSI_COLOR_RESET, global_success_count);
    printf(ANSI_COLOR_GREEN "\tALL UNIT TEST PASSED\n" ANSI_COLOR_RESET);
  } else {
    printf(ANSI_COLOR_RED "\n\tSOME UNIT TESTS FAILED\n" ANSI_COLOR_RESET);
  }
}


void TestFasta() {
  printf("Running TestFasta:\n");
  size_t len = 1000;
  char * seq = new char[len];
  for (size_t i = 0; i < len; i++)
    seq[i] = 'A';

  Utils::SaveFastaFile((char*)"tmp_file.fasta", seq, len);
  size_t len2;
  char * seq2;
  Utils::ReadFastaFile((char*)"tmp_file.fasta", &seq2, &len2);
  if (len != len2) {
    Fail();
    return;
  }
  for (size_t i = 0; i < len; i++) {
    if (seq[i] != seq2[i]) {
      Fail();
      return;
    }
  }
  Success();
}

////// thoe were for dev / test.
void TestPhaserShortK_A() {
  printf("Running TestPhaserShortK_A:\n");
  char M1[2] = {'1', '1'};
  char M2[2] = {'2', '2'};
  size_t M_len = 2;

  char F1[2] = {'8', '8'};
  char F2[2] = {'9', '9'};
  size_t F_len = 2;

  char C1[2] =        {'1', '1' };
  char C2[2] =        {'8', '8'};
  char phase_real[2] = {'0', '0'};
  size_t C_len = 2;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (4*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserShortK_B() {
  printf("Running TestPhaserShortK_B:\n");
  char M1[1] = {'1' };
  char M2[1] = {'2' };
  size_t M_len = 1;

  char F1[1] = {'8'};
  char F2[1] = {'9'};
  size_t F_len = 1;

  char C1[1] =        {'1'};
  char C2[1] =        {'8'};
  char phase_real[1] = {'0' };
  size_t C_len = 1;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserShortK_C() {
  printf("Running TestPhaserShortK_C:\n");
  char M1[1] = {'1'};
  char M2[1] = {'2'};
  size_t M_len = 1;

  char F1[2] = {'8', '8'};
  char F2[2] = {'9', '7'};
  size_t F_len = 2;

  char C1[1] =        {'7'};
  char C2[1] =        {'1'};
  char phase_real[1] = {'1' };
  size_t C_len = 1;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*SCORE_MATCH + 1* SCORE_GAP)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}


void TestPhaserShortK_D() {
  printf("Running TestPhaserShortK_D:\n");
  char M1[2] = {'1', '3' };
  char M2[2] = {'2', '4'};
  size_t M_len = 2;

  char F1[2] = {'6', '8'};
  char F2[2] = {'7', '9'};
  size_t F_len = 2;

  char C1[1] =        {'2'};
  char C2[1] =        {'9'};
  char phase_real[1] = {'0' };
  size_t C_len = 1;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*SCORE_MATCH + 2*SCORE_GAP)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}


void TestPhaserSimple() {
  printf("Running TestPhaserSimple:\n");
  char M1[10] = {'1', '1', '1', '1', '1', '1', '1', '1', '1', '1'};
  char M2[10] = {'2', '2', '2', '2', '2', '2', '2', '2', '2', '2'};
  size_t M_len = 10;

  char F1[10] = {'8', '8', '8', '8', '8', '8', '8', '8', '8', '8'};
  char F2[10] = {'9', '9', '9', '9', '9', '9', '9', '9', '9', '9'};
  size_t F_len = 10;

  char C1[10] =        {'5', '1', '1', '1', '1', '1', '1', '1', '1', '1'};
  char C2[10] =        {'8', '8', '8', '8', '8', '8', '8', '8', '8', '8'};
  char phase_real[10] = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};
  size_t C_len = 10;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (19*SCORE_MATCH + 1* SCORE_MISMATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_A() {
  printf("Running TestPhaserIdentical_A:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 5, 7};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;


  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_B() {
  printf("Running TestPhaserIdentical_B:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {2, 2, 2, 2};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 0
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_C() {
  printf("Running TestPhaserIdentical_C:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 5, 7};
  char C2[4] = {9, 9, 9, 9};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 0

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_D() {
  printf("Running TestPhaserIdentical_D:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {2, 2, 2, 2};
  char C2[4] = {9, 9, 9, 9};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_E() {
  printf("Running TestPhaserIdentical_E:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {8, 8, 8, 8};
  char C2[4] = {1, 3, 5, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_F() {
  printf("Running TestPhaserIdentical_F:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {8, 8, 8, 8};
  char C2[4] = {2, 2, 2, 2};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_G() {
  printf("Running TestPhaserIdentical_G:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 9, 9, 9};
  char C2[4] = {1, 3, 5, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 1

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserIdentical_H() {
  printf("Running TestPhaserIdentical_H:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 9, 9, 9};
  char C2[4] = {2, 2, 2, 2};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}
//////
void TestPhaserOnlyMotherRecomb_A() {
  printf("Running TestPhaserOnlyMotherRecomb_A:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 2, 2};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_B() {
  printf("Running TestPhaserOnlyMotherRecomb_B:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 2, 5, 2};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 0
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_C() {
  printf("Running TestPhaserOnlyMotherRecomb_C:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 2, 2};
  char C2[4] = {9, 9, 9, 9};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 0

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_D() {
  printf("Running TestPhaserOnlyMotherRecomb_D:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 2, 5, 2};
  char C2[4] = {9, 9, 9, 9};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_E() {
  printf("Running TestPhaserOnlyMotherRecomb_E:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {8, 8, 8, 8};
  char C2[4] = {1, 2, 2, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_F() {
  printf("Running TestPhaserOnlyMotherRecomb_F:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {8, 8, 8, 8};
  char C2[4] = {1, 2, 2, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_G() {
  printf("Running TestPhaserOnlyMotherRecomb_G:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 9, 9, 9};
  char C2[4] = {1, 2, 5, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 1

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserOnlyMotherRecomb_H() {
  printf("Running TestPhaserOnlyMotherRecomb_H:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 9, 9, 9};
  char C2[4] = {2, 2, 2, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

//////
void TestPhaserFatherMotherRecomb_A() {
  printf("Running TestPhaserFatherMotherRecomb_A:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 2, 2};
  char C2[4] = {8, 9, 9, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_B() {
  printf("Running TestPhaserFatherMotherRecomb_B:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 2, 5, 2};
  char C2[4] = {9, 8, 9, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 0
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_C() {
  printf("Running TestPhaserFatherMotherRecomb_C:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 2, 2};
  char C2[4] = {9, 9, 8, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 0

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_D() {
  printf("Running TestPhaserFatherMotherRecomb_D:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 2, 5, 2};
  char C2[4] = {8, 8, 8, 9};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_E() {
  printf("Running TestPhaserFatherMotherRecomb_E:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 8, 9, 8};
  char C2[4] = {1, 2, 2, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_F() {
  printf("Running TestPhaserFatherMotherRecomb_F:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {8, 9, 8, 9};
  char C2[4] = {1, 2, 2, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_G() {
  printf("Running TestPhaserFatherMotherRecomb_G:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 8, 9, 9};
  char C2[4] = {1, 2, 5, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 1

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserFatherMotherRecomb_H() {
  printf("Running TestPhaserFatherMotherRecomb_H:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 9, 8, 9};
  char C2[4] = {2, 2, 2, 7};
  char phase_real[4] = {'1', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_A() {
  printf("Running TestPhaserAllRecomb_A:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 3, 9, 8};
  char C2[4] = {8, 9, 2, 2};
  char phase_real[4] = {'0', '0', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_B() {
  printf("Running TestPhaserAllRecomb_B:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 2, 9, 8};
  char C2[4] = {9, 8, 5, 2};
  char phase_real[4] = {'0', '0', '1', '1'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 0
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_C() {
  printf("Running TestPhaserAllRecomb_C:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 9, 8, 2};
  char C2[4] = {9, 3, 2, 8};
  char phase_real[4] = {'0', '1', '1', '0'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 0

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_D() {
  printf("Running TestPhaserAllRecomb_D:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 8, 5, 9};
  char C2[4] = {8, 2, 8, 2};
  char phase_real[4] = {'0', '1', '0', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_E() {
  printf("Running TestPhaserAllRecomb_E:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {1, 8, 9, 8};
  char C2[4] = {9, 2, 2, 7};
  char phase_real[4] = {'0', '1', '1', '1'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_F() {
  printf("Running TestPhaserAllRecomb_F:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {8, 9, 2, 9};
  char C2[4] = {1, 2, 8, 7};
  char phase_real[4] = {'1', '1', '0', '1'};
  size_t C_len = 4;
  // mf = 1
  // ff = 0
  // cc = 1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_G() {
  printf("Running TestPhaserAllRecomb_G:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {9, 8, 5, 7};
  char C2[4] = {1, 2, 9, 9};
  char phase_real[4] = {'1', '1', '0', '0'};
  size_t C_len = 4;
  // mf = 0
  // ff = 1
  // cc = 1

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserAllRecomb_H() {
  printf("Running TestPhaserAllRecomb_H:\n");
  char M1[4] = {1, 3, 5, 7};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {2, 9, 8, 7};
  char C2[4] = {9, 2, 2, 9};
  char phase_real[4] = {'0', '1', '1', '0'};
  size_t C_len = 4;

  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  if (score != (2*((int)C_len)*SCORE_MATCH)) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserMismatchIndel() {
  printf("Running TestPhaserMismatchIndel:\n");
  char M1[4] = {1, 1, 1, 1};
  char M2[4] = {2, 2, 2, 2};
  size_t M_len = 4;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[4] = {5, '-', 1, 1};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0', '0'};
  size_t C_len = 4;
  // mf = 0
  // ff = 0
  // cc = 0
  // 1 mismatch + 1 deletion
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 2)*SCORE_MATCH + 1*SCORE_MISMATCH + 1*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserParentsGapped_A() {
  printf("Running TestPhaserParentsGapped_A:\n");

  char M1[7] = {1,  1 , 1, 1, 1, '-', 1};
  char M2[7] = {2,  2 , 2, 2, 3,  5 , 2};
  size_t M_len = 7;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[6] = {1, 1, 1, 1, 1, 1};
  char C2[6] = {8, 8, 8, 8, 8, 8};
  char phase_real[6] = {'0', '0', '0', '0' , '0', '0'};
  size_t C_len = 6;

  // only 2 deletions to  align C2 to F
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 2)*SCORE_MATCH + 2*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserParentsGapped_B() {
  printf("Running TestPhaserParentsGapped_B:\n");
  char M1[7] = {1, '-', 1, 1, 3, 3  , 1};
  char M2[7] = {2,   2, 2, 2, 3, '-', 2};
  size_t M_len = 7;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[5] = {1, 1, 1, 3, 1};
  char C2[5] = {8, 8, 8, 8, 8};
  char phase_real[5] = {'0', '0', '0' , '0', '0'};
  size_t C_len = 5;
  // only 1 deletions to  align C2 to F
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 1)*SCORE_MATCH + 1*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserParentsGapped_C() {
  printf("Running TestPhaserParentsGapped_C:\n");
  char M1[2] = {1, '-'};
  char M2[2] = {2,   2};
  size_t M_len = 2;

  char F1[2] = {8, 8};
  char F2[2] = {9, 9};
  size_t F_len = 2;

  char C1[2] = {1, 1};
  char C2[2] = {8, 8};
  size_t C_len = 2;
  char phase_real[2] = {'0', '0'};
  // only 1 deletions to  align C2 to F
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 1)*SCORE_MATCH + 1*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserParentsGapped_D() {
  printf("Running TestPhaserParentsGapped_D:\n");
  char M1[7] = {1, '-', 1, 1, 3, 3  , 1};
  char M2[7] = {2,   2, 2, 2, 3, '-', 2};
  size_t M_len = 7;

  char F1[4] = {8, 8, 8, 8};
  char F2[4] = {9, 9, 9, 9};
  size_t F_len = 4;

  char C1[6] = {1, 2, 8, 1, 8, 8};
  char C2[6] = {8, 8, 1, 8, 3, 1};
  char phase_real[6] = {'0', '0', '1', '0' , '1', '1'};
  size_t C_len = 6;


  // only 2 deletions to  align C2 to F
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 2)*SCORE_MATCH + 2*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserChildGapped_A() {
  printf("Running TestPhaserChildGapped_A:\n");
  char M1[2] = {1, 1 };
  char M2[2] = {2, 2 };
  size_t M_len = 2;

  char F1[3] = {8, 8, 8};
  char F2[3] = {9, 9, 9};
  size_t F_len = 3;

  char C1[3] = {1, '-', 1};
  char C2[3] = {8,  8,  8};
  char phase_real[3] = {'0', '0', '0'};
  size_t C_len = 3;
  // 2 deletions from C2.
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 1)*SCORE_MATCH;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserChildGapped_B() {
  printf("Running TestPhaserChildGapped_B:\n");
  char M1[2] = {1, 1 };
  char M2[2] = {2, 2 };
  size_t M_len = 2;

  char F1[2] = {8, 8};
  char F2[2] = {9, 9};
  size_t F_len = 2;

  char C1[3] = {1, '-', 1};
  char C2[3] = {8,  8,  8};
  char phase_real[3] = {'0', '0', '0'};
  size_t C_len = 3;
  // 2 deletions from C2.
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 2)*SCORE_MATCH + 1*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}


void TestPhaserChildGapped_C() {
  printf("Running TestPhaserChildGapped_C:\n");
  const size_t M_len = 4;
  char M1[M_len] = {1, 1, 1, 1 };
  char M2[M_len] = {2, 2, 2, 2 };

  const size_t F_len = 4;
  char F1[F_len] = {8, 8, 8, 8};
  char F2[F_len] = {9, 9, 9, 9};

  const size_t C_len = 6;
  char C1[C_len] = {1, 1, '-', '-', 1, 1};
  char C2[C_len] = {8, 8,  8 ,  8 , 8, 8};
  char phase_real[C_len] = {'0', '0', '0' , '0', '0', '0'};
  // 2 deletions from C2.
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected =
      (2*(int)C_len - 4)*SCORE_MATCH + 2*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}



void TestPhaserPartialDeleteChild_A() {
  printf("Running TestPhaserPartialDeleteChild_A:\n");
  char M1[2] = {2, 1};
  char M2[2] = {4, 4};
  size_t M_len = 2;

  char F1[6] = {8, 8, 7, 7, 8, 8 };
  char F2[6] = {9, 9, 6, 6, 9, 9 };
  size_t F_len = 6;

  char C1[4] = {2, 3, 3, 1};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0' , '0'};
  size_t C_len = 4;
  // no recombination.
  // two deletions  in C1,
  // two insertions in F1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected = 6*SCORE_MATCH + 4*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserPartialDeleteChild_B() {
  printf("Running TestPhaserPartialDeleteChild_B:\n");
  char M1[2] = {2, 1};
  char M2[2] = {4, 4};
  size_t M_len = 2;

  char F1[6] = {8, 8, 8, 8, 7, 7 };
  char F2[6] = {9, 9, 6, 6, 9, 9 };
  size_t F_len = 6;

  char C1[4] = {2, 3, 3, 1};
  char C2[4] = {8, 8, 8, 8};
  char phase_real[4] = {'0', '0', '0' , '0'};
  size_t C_len = 4;
  // no recombination.
  // two deletions  in C1,
  // two insertions in F1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected = 6*SCORE_MATCH + 4*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}

void TestPhaserPartialDeleteChild_C() {
  printf("Running TestPhaserPartialDeleteChild_C:\n");
  char M1[1] = {1};
  char M2[1] = {4};
  size_t M_len = 1;

  char F1[2] = {8, 8};
  char F2[2] = {9, 9};
  size_t F_len = 2;

  char C1[2] = {3, 1};
  char C2[2] = {8, 8};
  char phase_real[4] = {'0', '0'};
  size_t C_len = 2;
  // no recombination.
  // two deletions  in C1,
  // two insertions in F1
  Phaser * tmp =  new Phaser(M1, M2, M_len,
                             F1, F2, F_len,
                             C1, C2, C_len);
  tmp->SetScoreGap(SCORE_GAP);
  tmp->SetScoreMismatch(SCORE_MISMATCH);
  tmp->SetScoreMatch(SCORE_MATCH);
  score_t score = tmp->similarity_and_phase();
  char * phase_algor = tmp->GetPhaseString();
  score_t expected = 3*SCORE_MATCH + 1*SCORE_GAP;
  if (score != expected) {
    Fail();
    return;
  }
  if (!equalPhases(phase_real, phase_algor, C_len)) {
    Fail();
    return;
  }
  delete(tmp);
  Success();
}


void TestPhaserExhaustiveSameLength() {
  printf("Running TestPhaserExhaustiveSameLength:\n");
  size_t M_len = 10;
  size_t F_len = 10;
  size_t C_len = 10;
  size_t n_repeats = 5;

  char * M1 = new char[M_len];
  char * M2 = new char[M_len];
  for (size_t i = 0; i < M_len; i++) {
    M1[i] = '1';
    M2[i] = '2';
  }
  char * F1 = new char[F_len];
  char * F2 = new char[F_len];
  for (size_t i = 0; i < F_len; i++) {
    F1[i] = '8';
    F2[i] = '9';
  }

  char * C1 = new char[C_len];
  char * C2 = new char[C_len];
  char * phase_real = new char[C_len];

  for (size_t dels = 0; dels <= C_len/3; dels++) {
    for (size_t mismatches = 0; mismatches <= C_len/2; mismatches++) {
      for (size_t counter = 0; counter < n_repeats; counter++) {
        for (size_t i = 0; i < C_len; i++) {
          bool coin_1 = rand()%2;
          bool coin_2 = rand()%2;
          bool coin_3 = rand()%2;
          char a = coin_1 ? '1' : '2';
          char b = coin_2 ? '8' : '9';

          C1[i] = coin_3 ? a : b;
          C2[i] = coin_3 ? b : a;
          phase_real[i] = coin_3 ? '0' : '1';
        }

        for (size_t i = 0; i < mismatches; i++) {
          bool coin_1 = rand()%2;
          if (coin_1)
            C1[2*i] = '5';
          else
            C2[2*i] = '5';
        }
        for (size_t i = 0; i < dels; i++) {
          bool coin_1 = rand()%2;
          if (coin_1)
            C1[2*i + 1] =  '-';
          else
            C2[2*i + 1] =  '-';
        }

        Phaser * tmp =  new Phaser(M1, M2, M_len,
                                   F1, F2, F_len,
                                   C1, C2, C_len);
        tmp->SetScoreGap(SCORE_GAP);
        tmp->SetScoreMismatch(SCORE_MISMATCH);
        tmp->SetScoreMatch(SCORE_MATCH);
        score_t score = tmp->similarity_and_phase();
        char * phase_algor = tmp->GetPhaseString();
        score_t expected = (2*(int)C_len - (int)mismatches - (int)dels)*SCORE_MATCH +
            (int)mismatches* SCORE_MISMATCH + (int)dels * SCORE_GAP;
        if (score != expected) {
          printf("Wrong score\n");
          Debug::PrintArray(M1, M_len);
          Debug::PrintArray(M2, M_len);
          Debug::PrintLine(C_len);
          Debug::PrintArray(F1, F_len);
          Debug::PrintArray(F2, F_len);
          Debug::PrintLine(C_len);
          Debug::PrintArray(C1, C_len);
          Debug::PrintArray(C2, C_len);
          Fail();
          return;
        }
        if (!equalPhases(phase_real, phase_algor, C_len)) {
          printf("Wrong phases\n");
          Debug::PrintArray(M1, M_len);
          Debug::PrintArray(M2, M_len);
          Debug::PrintLine(C_len);
          Debug::PrintArray(F1, F_len);
          Debug::PrintArray(F2, F_len);
          Debug::PrintLine(C_len);
          Debug::PrintArray(C1, C_len);
          Debug::PrintArray(C2, C_len);
          Fail();
          return;
        }
        delete(tmp);
      }
    }
  }
  delete[] M1;
  delete[] M2;
  delete[] F1;
  delete[] F2;
  delete[] C1;
  delete[] C2;
  delete[] phase_real;
  Success();
}

void TestPhaserExhaustiveVarLength() {
  printf("Running TestPhaserExhaustiveVarLength:\n");
  size_t M_len = 40;
  size_t F_len = 20;
  bool M_is_larger = M_len > F_len ? 1 : 0;
  size_t C_len = 30;
  size_t n_repeats = 5;

  char * M1 = new char[M_len];
  char * M2 = new char[M_len];
  for (size_t i = 0; i < M_len; i++) {
    M1[i] = '1';
    M2[i] = '2';
  }
  char * F1 = new char[F_len];
  char * F2 = new char[F_len];
  for (size_t i = 0; i < F_len; i++) {
    F1[i] = '8';
    F2[i] = '9';
  }

  char * C1 = new char[C_len];
  char * C2 = new char[C_len];
  char * phase_real = new char[C_len];

  for (size_t mismatches = 0; mismatches <= C_len/4; mismatches++) {
    for (size_t counter = 0; counter < n_repeats; counter++) {
      for (size_t i = 0; i < C_len; i++) {
        bool coin_1 = rand()%2;
        bool coin_2 = rand()%2;
        bool coin_3 = rand()%2;
        char a = coin_1 ? '1' : '2';
        char b = coin_2 ? '8' : '9';

        C1[i] = coin_3 ? a : b;
        C2[i] = coin_3 ? b : a;
        phase_real[i] = coin_3 ? '0' : '1';
      }
      size_t real_mismatches = 0;
      for (size_t i = 0; i < mismatches; i++) {
        bool coin_1 = rand()%2;
        char original;
        if (coin_1) {
          original = C1[2*i];
          C1[2*i] =  '5';
        } else {
          original = C2[2*i];
          C2[2*i] =  '5';
        }
        if (M_is_larger) {
          if (original == '1' || original == '2')
            real_mismatches++;
        } else {
          if (original == '8' || original == '9')
            real_mismatches++;
        }
      }
      Phaser * tmp =  new Phaser(M1, M2, M_len,
                                 F1, F2, F_len,
                                 C1, C2, C_len);
      tmp->SetScoreGap(SCORE_GAP);
      tmp->SetScoreMismatch(SCORE_MISMATCH);
      tmp->SetScoreMatch(SCORE_MATCH);

      size_t m_match = std::min(C_len, M_len);
      size_t f_match = std::min(C_len, F_len);
      size_t t_match = m_match + f_match - real_mismatches;
      size_t gaps = std::max(C_len, M_len) + std::max(C_len, F_len) - m_match - f_match;
      score_t expected = ((int)t_match)*SCORE_MATCH +
                         (int)gaps*SCORE_GAP + (int)real_mismatches*SCORE_MISMATCH;
      score_t score = tmp->similarity_and_phase();
      char * phase_algor = tmp->GetPhaseString();

      if (score != expected) {
        printf("Wrong score\n");
        Debug::PrintArray(M1, M_len);
        Debug::PrintArray(M2, M_len);
        Debug::PrintLine(C_len);
        Debug::PrintArray(F1, F_len);
        Debug::PrintArray(F2, F_len);
        Debug::PrintLine(C_len);
        Debug::PrintArray(C1, C_len);
        Debug::PrintArray(C2, C_len);
        Fail();
        return;
      }
      if (!equalPhases(phase_real, phase_algor, C_len)) {
        printf("Wrong phases\n");
        Debug::PrintArray(M1, M_len);
        Debug::PrintArray(M2, M_len);
        Debug::PrintLine(C_len);
        Debug::PrintArray(F1, F_len);
        Debug::PrintArray(F2, F_len);
        Debug::PrintLine(C_len);
        Debug::PrintArray(C1, C_len);
        Debug::PrintArray(C2, C_len);
        Fail();
        return;
      }
      delete(tmp);
    }
  }
  delete[] M1;
  delete[] M2;
  delete[] F1;
  delete[] F2;
  delete[] C1;
  delete[] C2;
  delete[] phase_real;
  Success();
}


int main() {
  // The following asseertions are not necessary in general,
  // but they are the sensible option, and we use them to
  // calculate expeted answers in the tests:
  assert(SCORE_MISMATCH <= SCORE_GAP);  // Prefer gap over mismatch
  assert(2*SCORE_GAP <= SCORE_MISMATCH);  // Prefer mismatch over two gaps

  bool special = 0;
  bool exhaustive = 1;
  if (special) {
  } else {
    TestFasta();
    
    TestPhaserShortK_A();
    TestPhaserShortK_B();
    TestPhaserShortK_C();
    TestPhaserShortK_D();
    TestPhaserChildGapped_A();
    TestPhaserChildGapped_B();
    TestPhaserChildGapped_C();

    TestPhaserSimple();

    TestPhaserIdentical_A();
    TestPhaserIdentical_B();
    TestPhaserIdentical_C();
    TestPhaserIdentical_D();
    TestPhaserIdentical_E();
    TestPhaserIdentical_F();
    TestPhaserIdentical_G();
    TestPhaserIdentical_H();

    TestPhaserOnlyMotherRecomb_A();
    TestPhaserOnlyMotherRecomb_B();
    TestPhaserOnlyMotherRecomb_C();
    TestPhaserOnlyMotherRecomb_D();
    TestPhaserOnlyMotherRecomb_E();
    TestPhaserOnlyMotherRecomb_F();
    TestPhaserOnlyMotherRecomb_G();
    TestPhaserOnlyMotherRecomb_H();

    TestPhaserFatherMotherRecomb_A();
    TestPhaserFatherMotherRecomb_B();
    TestPhaserFatherMotherRecomb_C();
    TestPhaserFatherMotherRecomb_D();
    TestPhaserFatherMotherRecomb_E();
    TestPhaserFatherMotherRecomb_F();
    TestPhaserFatherMotherRecomb_G();
    TestPhaserFatherMotherRecomb_H();

    TestPhaserAllRecomb_A();
    TestPhaserAllRecomb_B();
    TestPhaserAllRecomb_C();
    TestPhaserAllRecomb_D();
    TestPhaserAllRecomb_E();
    TestPhaserAllRecomb_F();
    TestPhaserAllRecomb_G();
    TestPhaserAllRecomb_H();

    TestPhaserMismatchIndel();
    TestPhaserPartialDeleteChild_A();
    TestPhaserPartialDeleteChild_B();
    TestPhaserPartialDeleteChild_C();

    TestPhaserParentsGapped_A();
    TestPhaserParentsGapped_B();
    TestPhaserParentsGapped_C();
    TestPhaserParentsGapped_D();

    if (exhaustive) {
      TestPhaserExhaustiveSameLength();
      TestPhaserExhaustiveVarLength();
    }
  }
  Summary();
}
