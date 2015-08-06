/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */

#include "./phaser.h"
#include <stdio.h>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include "./basic.h"
#include "./debug.h"

Phaser::Phaser(char * _M1,
               char * _M2,
               size_t  _M_len,
               char * _F1,
               char * _F2,
               size_t _F_len,
               char * _C1,
               char * _C2,
               size_t _C_len) {
  M1 = _M1;
  M2 = _M2;
  M_len = _M_len;
  F1 = _F1;
  F2 = _F2;
  F_len = _F_len;
  C1 = _C1;
  C2 = _C2;
  C_len = _C_len;

  assert(C_len > 0);
  assert(M_len > 0);
  assert(F_len > 0);

  I_len = M_len;
  J_len = F_len;
  phase_string = new char[C_len];
  for (size_t i = 0; i < C_len; i++) {
    phase_string[i] = '?';
  }
  verbose = false;
}


// Basic DPA algorithm;
// O(n^3) time
// O(n^2) space
score_t Phaser::similarity() {
  size_t dumb_i;
  size_t dumb_j;
  size_t dumb_k;
  score_t answer = partial_aligner(0, 0, 0,
                                   M_len-1, F_len-1, C_len-1,
                                   &dumb_i, &dumb_j, &dumb_k);
  return answer;
}

score_t Phaser::similarity_and_phase() {
  score_t ans = aligner(0, 0, 0, M_len-1, F_len-1, C_len-1);
  PrintPhaseString();
  return ans;
}

score_t Phaser::aligner(size_t i_ini,
                        size_t j_ini,
                        size_t k_ini,
                        size_t i_end,
                        size_t j_end,
                        size_t k_end) {
  size_t i_med, j_med, k_med;
  score_t ans = partial_aligner(i_ini,
                                j_ini,
                                k_ini,
                                i_end,
                                j_end,
                                k_end,
                                &i_med,
                                &j_med,
                                &k_med);
  if (verbose) {
    printf("%lu, %lu, %lu \n", k_ini, k_med, k_end);
  }

  bool p1 = false;
  bool p2 = false;
  score_t ans_1, ans_2;
  if (k_ini <= k_med && k_med != k_end &&
      k_ini < k_med+1) {
    ans_1= aligner(i_ini,
                   j_ini,
                   k_ini,
                   i_med,
                   j_med,
                   k_med);
    p1 = true;
  }

  if (k_med+1 <= k_end &&
      k_ini < k_med+1) {
    ans_2= aligner(i_med + 1,
                   j_med + 1,
                   k_med + 1,
                   i_end ,
                   j_end,
                   k_end);
    p2 = true;
  }
  if (p1 && p2) {
    if (ans != ans_1 + ans_2) {
      fprintf(stderr, "This should never happen.\n");
      fprintf(stderr, "Inconsistency in recursive call, Phaser::aligner.\n");
      fprintf(stderr, "Please send us a report.\n");
      exit(-1);
    }
  } else {
    if (p1 || p2) {
      printf("im curious:\n");
      printf("%lu, %lu, %lu \n", k_ini, k_med, k_end);
    }
  }
  return ans;
}

// Checkpoint algorithm.
// O(n^3) time
// O(n^2) space
score_t Phaser::partial_aligner(size_t i_ini,
                                size_t j_ini,
                                size_t k_ini,
                                size_t i_end,
                                size_t j_end,
                                size_t k_end,
                                size_t *i_med,
                                size_t *j_med,
                                size_t *k_med) {
  assert(j_end >= j_ini || j_end + 1 == j_ini);
  assert(i_end >= i_ini || i_end + 1 == i_ini);
  assert(k_end >= k_ini);
  score_t * prev_face[8];
  score_t * curr_face[8];
  my_pair * prev_check[8];
  my_pair * curr_check[8];

  // Obs: we will use local i (resp. j, k) from 0 to I_len (J_len, K_len).
  // characters are stracted from i_ini+i (j, k resp.).
  // checkpoint answer is shifted back at the end.

  I_len = i_end - i_ini + 1;
  J_len = j_end - j_ini + 1;
  size_t K_len = k_end - k_ini + 1;
  size_t mid_k = K_len/2;
  for (size_t m = 0; m < 8; m++) {
    prev_face[m] = new score_t[(I_len+1) * (J_len+1)];
    prev_check[m] = new my_pair[(I_len+1) * (J_len+1)];
  }

  // 8 points:
  for (bool cf : {false, true}) {
    for (bool ff : {false, true}) {
      for (bool mf : {false, true}) {
        size_t m = m_index(mf, ff, cf);
        prev_face[m][IJ(0, 0)] = 0;
        prev_check[m][IJ(0, 0)] = my_pair(0, 0);
      }
    }
  }

  // 8 lines (j=0):
  for (size_t i = 1; i <= I_len; i++) {
    for (bool cf : {false, true}) {
      for (bool ff : {false, true}) {
        for (bool mf : {false, true}) {
          size_t m = m_index(mf, ff, cf);
          char m_char = mf ? M2[i_ini + i-1] : M1[i_ini + i-1];
          prev_face[m][IJ(i, 0)] = std::max(prev_face[m_index(0, ff, cf)][IJ(i-1, 0)] + score(m_char, '-'),  //  NOLINT
                                            prev_face[m_index(1, ff, cf)][IJ(i-1, 0)] + score(m_char, '-'));  //  NOLINT
          prev_check[m][IJ(i, 0)] = my_pair(i, 0);
        }
      }
    }
  }

  // 8 faces:
  for (size_t j = 1; j <= J_len; j++) {
    for (size_t i = 0; i <= I_len; i++) {
      for (bool cf : {false, true}) {
        for (bool ff : {false, true}) {
          for (bool mf : {false, true}) {
            size_t m = m_index(mf, ff, cf);
            char f_char = ff ? F2[j_ini + j-1] : F1[j_ini + j-1];
            prev_face[m][IJ(i, j)] = std::max(prev_face[m_index(mf, 0, cf)][IJ(i, j-1)] + score(f_char, '-'),  //  NOLINT
                                              prev_face[m_index(mf, 1, cf)][IJ(i, j-1)] + score(f_char, '-'));  //  NOLINT
            prev_check[m][IJ(i, j)] = my_pair(i, j);
          }
        }
      }
    }
  }
  PrintFace(prev_face);

  bool malloc_opt = true;
  if (malloc_opt) {
    for (size_t m = 0; m < 8; m++) {
      curr_face[m] = new score_t[(I_len+1) * (J_len+1)];
      curr_check[m] = new my_pair[(I_len+1) * (J_len+1)];
    }
  }

  // the rest of the faces:
  for (size_t k = 1; k <= K_len; k++) {
    if (!malloc_opt) {
      for (size_t m = 0; m < 8; m++) {
        curr_face[m] = new score_t[(I_len+1) * (J_len+1)];
        curr_check[m] = new my_pair[(I_len+1) * (J_len+1)];
      }
    }
    for (size_t j = 0; j <= J_len; j++) {
      for (size_t i = 0; i <= I_len; i++) {
        for (bool cf : {false, true}) {
          for (bool ff : {false, true}) {
            for (bool mf : {false, true}) {
              // m_char anf f_char should not be used, at least i > 0  (j > 0).
              // If that is not the case, we initilize with a non-accepted character that
              // will trig an error if used.
              char m_char, f_char;
              if (i > 0) {
                m_char = mf ? M2[i_ini + i-1] : M1[i_ini + i-1];
              } else {
                m_char = 'J';  // Not in gen alphabet, will trigger an error if used.
              }
              if (j > 0) {
                f_char = ff ? F2[j_ini + j-1] : F1[j_ini + j-1];
              } else {
                f_char = 'J';  // Not in gen alphabet, will trigger an error if used.
              }

              char c_1    = cf ? C2[k_ini + k-1] : C1[k_ini + k-1];
              char c_2    = cf ? C1[k_ini + k-1] : C2[k_ini + k-1];

              UpdateGeneral(curr_face,
                            prev_face,
                            curr_check,
                            prev_check,
                            i,
                            j,
                            k,
                            mid_k,
                            mf,
                            ff,
                            cf,
                            m_char,
                            f_char,
                            c_1,
                            c_2);
            }
          }
        }
      }
    }

    for (size_t m = 0; m < 8; m++) {
      if (malloc_opt) {
        score_t * tmp_face = prev_face[m];
        my_pair * tmp_check = prev_check[m];
        prev_face[m] = curr_face[m];
        prev_check[m] = curr_check[m];
        curr_face[m] = tmp_face;
        curr_check[m] = tmp_check;
      } else {        delete[] (prev_face[m]);
        delete[] (prev_check[m]);
        prev_face[m] = curr_face[m];
        prev_check[m] = curr_check[m];
      }
    }
    if (k >= mid_k) {
      // verbose = true;
    }
    PrintFace(prev_face);
    PrintCheck(prev_check);
  }

  // we use char_i = M[i-1]
  for (int i = 0; i < 8; i++) {
    assert(prev_check[i][IJ(I_len, J_len)].first <= I_len);
    assert(prev_check[i][IJ(I_len, J_len)].second <= J_len);

    prev_check[i][IJ(I_len, J_len)].first--;
    prev_check[i][IJ(I_len, J_len)].second--;
  }
  mid_k--;


  // extract max:
  score_t ans;
  bool flip_ans;
  ExtractMax(prev_face, prev_check, &ans, i_med, j_med, &flip_ans);
  *k_med = k_ini + mid_k;
  *i_med = i_ini + (*i_med);
  *j_med = j_ini + (*j_med);

  assert(CorrectIniMedEnd(i_ini, *i_med, i_end));
  assert(CorrectIniMedEnd(j_ini, *j_med, j_end));
  assert(CorrectIniMedEnd(k_ini, *k_med, k_end));


  char phase_char = flip_ans ? '1' : '0';
  if (phase_string[k_end] != '?') {
    assert(phase_string[k_end] == phase_char);
  }
  if (phase_string[k_end] == '?') {
    phase_string[k_end] = phase_char;
  }

  for (size_t m = 0; m < 8; m++) {
    delete[] (prev_face[m]);
    delete[] (prev_check[m]);
    if (malloc_opt) {
      delete[] (curr_face[m]);
      delete[] (curr_check[m]);
    }
  }
  return ans;
}

// TODO(possible optimization):
// To have two versions, one that keep track,
// and one who does not, to avoid some work for k < mid_k ?

// With the current scheme we store the larger i, j
// that can be aligned to mid_k in the optimal alignment.
void Phaser::UpdateGeneral(score_t ** curr_face,
                           score_t ** prev_face,
                           my_pair ** curr_check,
                           my_pair ** prev_check,
                           size_t i,
                           size_t j,
                           size_t k,
                           size_t mid_k,
                           bool mf,
                           bool ff,
                           bool cf,
                           char m_char,
                           char f_char,
                           char c_1,
                           char c_2) {
  size_t m = m_index(mf, ff, cf);
  score_t max_score;
  my_pair max_check;

  // only k decreases. Two deletions from C.
  score_t c_ins_1 =  prev_face[m_index(mf, ff, 0)][IJ(i, j)] + score(c_1, '-') + score(c_2, '-');
  score_t c_ins_2 =  prev_face[m_index(mf, ff, 1)][IJ(i, j)] + score(c_1, '-') + score(c_2, '-');
  my_pair check_1 = prev_check[m_index(mf, ff, 0)][IJ(i, j)];
  my_pair check_2 = prev_check[m_index(mf, ff, 1)][IJ(i, j)];

  if (c_ins_1 >= c_ins_2) {
    max_score = c_ins_1;
    max_check = check_1;
  } else {
    max_score = c_ins_2;
    max_check = check_2;
  }

  if (i > 0) {
    if (m_char == '-') {
      score_t p1 = curr_face[m_index(0, ff, cf)][IJ(i-1, j)];
      score_t p2 = curr_face[m_index(1, ff, cf)][IJ(i-1, j)];
      curr_face[m][IJ(i, j)] = std::max(p1, p2);
      if (k == mid_k) {
        curr_check[m][IJ(i, j)] = my_pair(i, j);
      } else if (k > mid_k) {
        curr_check[m][IJ(i, j)] = (p1 > p2) ?
            curr_check[m_index(0, ff, cf)][IJ(i-1, j)] :
            curr_check[m_index(1, ff, cf)][IJ(i-1, j)];
      }
      return;
    }
    // only i decreases. Single deletion from M.
    // TODO(Readability): This might be a one-level for over pre_mf.
    score_t ins_val;
    my_pair ins_check;
    ins_val   =  curr_face[m_index(0, ff, cf)][IJ(i-1, j)] + score(m_char, '-');
    ins_check = curr_check[m_index(0, ff, cf)][IJ(i-1, j)];
    UpdateVals(ins_val, ins_check, &max_score, &max_check);

    ins_val  =   curr_face[m_index(1, ff, cf)][IJ(i-1, j)] + score(m_char, '-');
    ins_check = curr_check[m_index(1, ff, cf)][IJ(i-1, j)];
    UpdateVals(ins_val, ins_check, &max_score, &max_check);

    // k and i decreases: single deletions from C, M aligns.
    for (bool pre_cf : {false, true}) {
      for (bool pre_mf : {false, true}) {
        score_t del_val  =   prev_face[m_index(pre_mf, ff, pre_cf)][IJ(i-1, j)] + score(c_1, m_char) + score(c_2, '-');  //  NOLINT
        my_pair del_check = prev_check[m_index(pre_mf, ff, pre_cf)][IJ(i-1, j)];
        UpdateVals(del_val, del_check, &max_score, &max_check);
      }
    }
  }

  if (j > 0) {
    if (f_char == '-') {
      score_t p1 = curr_face[m_index(mf, 0, cf)][IJ(i, j-1)];
      score_t p2 = curr_face[m_index(mf, 1, cf)][IJ(i, j-1)];
      curr_face[m][IJ(i, j)] = std::max(p1, p2);
      if (k == mid_k) {
        curr_check[m][IJ(i, j)] = my_pair(i, j);
      } else if (k > mid_k) {
        curr_check[m][IJ(i, j)] = (p1 > p2) ?
            curr_check[m_index(mf, 0, cf)][IJ(i, j-1)] :
            curr_check[m_index(mf, 1, cf)][IJ(i, j-1)];
      }
      return;
    }
    score_t ins_val;
    my_pair ins_check;
    // only j decreases. Single deletion from F.
    ins_val = curr_face[m_index(mf, 0, cf)][IJ(i, j-1)] + score(f_char, '-');
    ins_check = curr_check[m_index(mf, 0, cf)][IJ(i, j-1)];
    UpdateVals(ins_val, ins_check, &max_score, &max_check);

    ins_val = curr_face[m_index(mf, 1, cf)][IJ(i, j-1)] + score(f_char, '-');
    ins_check = curr_check[m_index(mf, 1, cf)][IJ(i, j-1)];
    UpdateVals(ins_val, ins_check, &max_score, &max_check);


    // k and j decreases: single deletions from C, F aligns.
    for (bool pre_cf : {false, true}) {
      for (bool pre_ff : {false, true}) {
        score_t del_val =   prev_face[m_index(mf, pre_ff, pre_cf)][IJ(i, j-1)] + score(c_1, '-') + score(c_2, f_char);  //  NOLINT
        my_pair del_check = prev_check[m_index(mf, pre_ff, pre_cf)][IJ(i, j-1)];
        UpdateVals(del_val, del_check, &max_score, &max_check);
      }
    }
  }

  if (i > 0 && j > 0) {
    for (bool pre_cf : {false, true}) {
      for (bool pre_ff : {false, true}) {
        for (bool pre_mf : {false, true}) {
          score_t aln_val = prev_face[m_index(pre_mf, pre_ff, pre_cf)][IJ(i-1, j-1)] + score(c_1, m_char) + score(c_2, f_char);  //  NOLINT
          my_pair aln_check = prev_check[m_index(pre_mf, pre_ff, pre_cf)][IJ(i-1, j-1)];
          UpdateVals(aln_val, aln_check, &max_score, &max_check);
        }
      }
    }  }

  curr_face[m][IJ(i, j)] = max_score;

  if (k == mid_k) {
    curr_check[m][IJ(i, j)] = my_pair(i, j);
  } else if (k > mid_k) {
    curr_check[m][IJ(i, j)] = max_check;
  }
}

Phaser::~Phaser() {
  delete[] phase_string;
}


///// Some debug functions:

void Phaser::PrintPhaseString() {
  if (verbose) {
    printf("\n");
    Debug::PrintArray(phase_string, C_len);
    printf("\n");
  }
}

void Phaser::PrintSequences() {
  printf("Phase:\n");
  Debug::PrintArray(phase_string, C_len);
  printf("\n");
  Debug::PrintLine(C_len);
  printf("M:\n");
  Debug::PrintArray(M1, M_len);
  Debug::PrintArray(M2, M_len);
  Debug::PrintLine(C_len);
  printf("F:\n");
  Debug::PrintArray(F1, F_len);
  Debug::PrintArray(F2, F_len);
  Debug::PrintLine(C_len);
  printf("C:\n");
  Debug::PrintArray(C1, C_len);
  Debug::PrintArray(C2, C_len);
  Debug::PrintLine(C_len);
  Debug::PrintLine(C_len);
  return;
}

void Phaser::PrintFace(score_t ** face) {
  const char separator    = ' ';
  const int width   = 5;
  if (verbose) {
    Debug::PrintLine(2*F_len/3);
    for (bool cf : {false, true}) {
      for (bool ff : {false, true}) {
        for (bool mf : {false, true}) {
          if (!cf && !ff && !mf) {
            Debug::PrintLine(2*F_len);
            printf("%i %i %i\n", mf, ff, cf);
            Debug::PrintLine(2*F_len);
            size_t m = m_index(mf, ff, cf);
            for (size_t i = 0; i <= I_len; i++) {
              for (size_t j = 0; j <= J_len; j++) {
                std::cout << std::left << std::setw(width) << std::setfill(separator) <<
                    face[m][IJ(i, j)];
              }
              std::cout << std::endl;
            }
          }
        }
      }
    }
    printf("\n");
  }
}

void Phaser::PrintCheck(my_pair ** check) {
  const char separator    = ' ';
  const int width   = 5;
  if (0) {
    Debug::PrintLine(2*F_len/3);
    for (bool cf : {false, true}) {
      for (bool ff : {false, true}) {
        for (bool mf : {false, true}) {
          if (!cf && !ff && !mf) {
            Debug::PrintLine(2*F_len);
            printf("%i %i %i\n", mf, ff, cf);
            Debug::PrintLine(2*F_len);
            size_t m = m_index(mf, ff, cf);
            for (size_t i = 0; i <= I_len; i++) {
              for (size_t j = 0; j <= J_len; j++) {
                std::cout << std::left << std::setw(width) << std::setfill(separator) <<
                    "(" << check[m][IJ(i, j)].first << "," << check[m][IJ(i, j)].second<< ")";
              }
              std::cout << std::endl;
            }
          }
        }
      }
    }
    printf("\n");
  }
}
