/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi

Input: Unphased father, unphased mother, unphased child.

Output: Phased child. Alignment among:
    Child's mother chromosome and the infered recombintation from the mather.
    Child's father chromosome and the infered recombintation from the father.

    If there were no mutations, but only recombination in the parents,
    the distance should be zero and the method is exact.
    If there were mutations our method looks for the feasible alignment
    that minimize the distance. (It might not be the actual one, if the actual
    process had more mutations that the minimum we obtained):
Inputs:

    M1  = AAAAAAA
    M2  = BBBBBBB

    F1  = CCCCCCC
    F2  = DDDDDDD

    C1  = AyAACCDD
    C2  = CCDDBBxB

    We infere: (2 mutations)
    CM  = AyAABBBX
    CF  = CCDDCCDD

    But the real might perfectly be (4 mutations):
    CM* = CyAABBBX
    CF* = AXDDCCDD

 */

#ifndef SRC_PHASER_H_
#define SRC_PHASER_H_

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cassert>
#include <utility>
#include "./basic.h"

typedef std::pair<size_t, size_t> my_pair;

class Phaser {
 protected:
  char * M1;
  char * M2;
  size_t M_len;
  char * F1;
  char * F2;
  size_t F_len;
  char * C1;
  char * C2;
  size_t C_len;

  char * phase_string;

  // local variables
  size_t I_len;
  size_t J_len;

  score_t SCORE_GAP;
  score_t SCORE_MISMATCH;
  score_t SCORE_MATCH;
  bool verbose;

 public:
  // constructor receive the input data.
  Phaser(char * _M1,
         char * _M2,
         size_t  _M_len,
         char * _F1,
         char * _F2,
         size_t _F_len,
         char * _C1,
         char * _C2,
         size_t _C_len);

  // Only computes similarity distance.
  score_t similarity();

  // Computes similarity distance, and
  // build the phase_string using the checkpoint method.
  score_t similarity_and_phase();

  // compute through partial_aligner and calls itself recursively.
  score_t aligner(size_t i_ini,
                  size_t j_ini,
                  size_t k_ini,
                  size_t i_end,
                  size_t j_end,
                  size_t k_end);

  // computes partial alignment and keep track of checkpoint.
  score_t partial_aligner(size_t i_ini,
                          size_t j_ini,
                          size_t k_ini,
                          size_t i_end,
                          size_t j_end,
                          size_t k_end,
                          size_t *i_med,
                          size_t *j_med,
                          size_t *k_med);
  // Auxiliar functions:

  void UpdateGeneral(score_t ** curr_face,
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
                     char c_2);

  // Inline methods :
  inline void ExtractMax(score_t ** prev_face,
                         my_pair ** prev_check,
                         score_t * ans,
                         size_t* i_med,
                         size_t* j_med,
                         bool * flip) {
    *ans = prev_face[0][IJ(I_len, J_len)];
    *i_med = ((prev_check[0][IJ(I_len, J_len)]).first);
    *j_med = ((prev_check[0][IJ(I_len, J_len)]).second);
    *flip = 0;
    for (bool cf : {false, true}) {
      for (bool ff : {false, true}) {
        for (bool mf : {false, true}) {
          size_t m = m_index(mf, ff, cf);
          if (prev_face[m][IJ(I_len, J_len)] > (*ans)) {
            *ans = prev_face[m][IJ(I_len, J_len)];
            *i_med = ((prev_check[m][IJ(I_len, J_len)]).first);
            *j_med = ((prev_check[m][IJ(I_len, J_len)]).second);
            *flip = cf;
          }
        }
      }
    }
  }

  inline void UpdateVals(score_t candidate,
                        my_pair check,
                        score_t * max,
                        my_pair * real_check) {
    if (candidate > *max) {
      *max = candidate;
      *real_check = check;
    }
  }

  inline char * GetPhaseString() {
    return phase_string;
  }

  inline bool CorrectIniMedEnd(size_t ini, size_t med, size_t end) {
    if (!(med <= end || med+1 <= end))
        return false;
    if (!(med+1 >= ini))
        return false;
    return true;
  }

  inline size_t m_index(bool m, bool f, bool c) {
    return (size_t)m + ((size_t)f << 1) + ((size_t)c << 2);
    /*
    size_t ans = 0;
    if (m) ans+=1;
    if (f) ans+=2;
    if (c) ans+=4;
    return ans;
    */
  }


  inline score_t score(char a, char b) {
    if ((a == '-') && (b == '-'))
      return 0;
    if (a == b)
      return SCORE_MATCH;
    if ((a == '-') || (b == '-'))
      return SCORE_GAP;
    else
      return SCORE_MISMATCH;
  }

  inline size_t IJ(size_t x, size_t y) {
    assert(x <= I_len);
    assert(y <= J_len);
    return (y * (I_len+1)) + x;
  }


  // Accesors and mutators:
  inline void SetScoreGap(score_t val) {
    assert(val < 0);
    SCORE_GAP = val;
  }
  inline void SetScoreMismatch(score_t val) {
    assert(val < 0);
    SCORE_MISMATCH = val;
  }
  inline void SetScoreMatch(score_t val) {
    assert(val > 0);
    SCORE_MATCH = val;
  }

  // Debug:
  void PrintSequences();
  void PrintFace(score_t ** face);
  void PrintCheck(my_pair ** check);
  void PrintPhaseString();


  ~Phaser();
};

#endif  // SRC_PHASER_H_
