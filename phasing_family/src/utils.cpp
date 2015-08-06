/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */

#include <sys/times.h>
#include <stdarg.h>
#include <unistd.h>
#include <inttypes.h>
#include <ctime>
#include <climits>
#include <cstdio>
#include <cassert>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include "./utils.h"
#include "./fasta.h"
#include "./basic.h"

using std::max;
using std::min;
using std::vector;
using std::pair;
using std::cout;
using std::left;
using std::endl;
using std::sort;

char gen_alph[4] = {'A', 'C', 'G', 'T'};
/* Time meassuring */
double ticks = (double)sysconf(_SC_CLK_TCK);
struct tms t1, t2;
void Utils::StartClock() {
  times(&t1);
}

double Utils::StopClock() {
  times(&t2);
  return (t2.tms_utime-t1.tms_utime)/ticks;
}

void Utils::CreateTrio(char * sequence,
                       size_t seq_len,
                       double error_ratio,
                       double long_indels,
                       double short_indels,
                       double repeat_ratio,
                       double str_ratio,
                       double parents_pm_ratio,
                       double child_pm_ratio,
                       char ** _M1, char ** _M2, size_t * _M_len,
                       char ** _F1, char ** _F2, size_t * _F_len,
                       char ** _C1, char ** _C2, size_t * _C_len) {
  // size_t n_cuts = 5;
  char* F1, *F2, *M1, *M2;
  size_t F_len, M_len;
  int m_dels,  f_dels;
  CreateGuy(sequence,
            seq_len,
            parents_pm_ratio,
            &M1,
            &M2,
            &m_dels);
  CreateGuy(sequence,
            seq_len,
            parents_pm_ratio,
            &F1,
            &F2,
            &f_dels);
  // Long INS need to be coordinated among M & F
  // in orfer to make the offsprings aligned by construction.
  // For order, we will do long DELS here:

  // double delete_fraction = 0.4 +  (0.2*rand() / (double)RAND_MAX);
  // assert(delete_fraction <= 0.61);
  double long_dels = long_indels/2;
  double long_ins  = long_indels/2;
  m_dels += Utils::MutateLONGDEL(M1, seq_len, long_dels);
  m_dels += Utils::MutateLONGDEL(M2, seq_len, long_dels);
  f_dels += Utils::MutateLONGDEL(F1, seq_len, long_dels);
  f_dels += Utils::MutateLONGDEL(F2, seq_len, long_dels);
  Utils::MutateLONGINS(&M1, &M2, &F1, &F2, &seq_len, long_ins);

  double short_dels = short_indels/2;
  double short_ins  = short_indels/2;
  m_dels += Utils::MutateSHORTDEL(M1, seq_len, short_dels);
  m_dels += Utils::MutateSHORTDEL(M2, seq_len, short_dels);
  f_dels += Utils::MutateSHORTDEL(F1, seq_len, short_dels);
  f_dels += Utils::MutateSHORTDEL(F2, seq_len, short_dels);
  Utils::MutateSHORTINS(&M1, &M2, &F1, &F2, &seq_len, short_ins);

  Utils::MutateSTR(&M1, &M2, &F1, &F2, &seq_len, str_ratio);
  Utils::MutateREPEATS(&M1, &M2, &F1, &F2, &seq_len, repeat_ratio);

  F_len = seq_len;
  M_len = seq_len;
  char *C1, *C2;
  size_t C_len;
  int non_inherited_dels, non_inherited_snps;
  CreateOffspring(M1, M2, M_len,
                  F1, F2, F_len,
                  child_pm_ratio,
                  &C1, &C2, &C_len,
                  &non_inherited_dels,
                  &non_inherited_snps);

  Utils::Compact(&M1, &M2, &M_len);
  Utils::Compact(&F1, &F2, &F_len);
  Utils::Compact(&C1, &C2, &C_len);

  MutateNoise(M1, M_len, error_ratio);
  MutateNoise(M2, M_len, error_ratio);
  MutateNoise(F1, F_len, error_ratio);
  MutateNoise(F2, F_len, error_ratio);
  MutateNoise(C1, C_len, error_ratio);
  MutateNoise(C2, C_len, error_ratio);


  *_M1 = M1;
  *_M2 = M2;
  *_M_len = M_len;
  *_F1 = F1;
  *_F2 = F2;
  *_F_len = F_len;
  *_C1 = C1;
  *_C2 = C2;
  *_C_len = C_len;
}

void Utils::CreateGuy(char * seq,
                      size_t seq_len,
                      double pm_ratio,
                      char ** g1,
                      char **g2,
                      int *n_dels) {
  int n_deletions = 0;
  // homozygous:
  (*g1) = Utils::CopySeq(seq, seq_len);
  n_deletions += Utils::MutatePointDEL(*g1, seq_len, pm_ratio/2);
  Utils::MutateSNP(*g1, seq_len, pm_ratio/2);
  (*g2) = Utils::CopySeq(*g1, seq_len);

  // heterozygous
  n_deletions += Utils::MutatePointDEL(*g1, seq_len, pm_ratio/4);
  Utils::MutateSNP(*g1, seq_len, pm_ratio/4);

  n_deletions += Utils::MutatePointDEL(*g2, seq_len, pm_ratio/4);
  Utils::MutateSNP(*g2, seq_len, pm_ratio/4);

  *n_dels = n_deletions;
}

// Follows:
// Weber, James L., and Carmen Wong.
// "Mutation of human short tandem repeats."
// Human molecular genetics 2.8 (1993): 1123-1128.
int Utils::MutateSTR(char **A ,
                     char **B,
                     char **C,
                     char **D,
                     size_t * len,
                     double mutation_ratio) {
  return MutateGENERIC_INSERT(A, B, C, D, len, mutation_ratio, 20, 60, true, true);
}

// Will not be used in the report.
int Utils::MutateREPEATS(char **A ,
                        char **B,
                        char **C,
                        char **D,
                        size_t * len,
                        double mutation_ratio) {
  return MutateGENERIC_INSERT(A, B, C, D, len, mutation_ratio, 5, 80, true, false);
}

// TODO(bio improvement): what is observed max size of "normal" repeat in the litetrature?
int Utils::MutateLONGINS(char **A ,
                         char **B,
                         char **C,
                         char **D,
                         size_t * len,
                         double mutation_ratio) {
  return MutateGENERIC_INSERT(A, B, C, D, len, mutation_ratio, min((size_t)50, (*len)/3), (*len)/2, false, false);
}

int Utils::MutateSHORTINS(char **A,
                          char **B,
                          char **C,
                          char **D,
                          size_t * len,
                          double mutation_ratio) {
  return MutateGENERIC_INSERT(A, B, C, D, len, mutation_ratio, 5, 49, false, false);
}

int Utils::MutateGENERIC_INSERT(char **A ,
                                char **B,
                                char **C,
                                char **D,
                                size_t * len,
                                double mutation_ratio,
                                size_t min_size,
                                size_t max_size,
                                bool do_repeat,
                                bool do_tandem) {
  if (do_tandem) assert(do_repeat);
  assert(max_size >= min_size);
  size_t min_indel = min_size;
  size_t max_indel = max_size;

  size_t to_insert = 4*(*len)*mutation_ratio;
  size_t new_len = (*len) + to_insert;
  vector<pair<size_t, size_t>> insert_lists;
  while (to_insert != 0) {
    max_indel = min(max_indel, to_insert);
    min_indel = min(min_indel, to_insert);
    assert(max_indel >= min_indel);

    size_t ins_len = min_indel + (size_t)rand()% (max_indel - min_indel + 1);

    assert(ins_len >= min_indel || ins_len == to_insert);
    assert(ins_len <= max_indel);
    size_t pos;
    if (do_repeat) {
      pos = ins_len + (size_t)rand()%((*len) - ins_len);
    } else {
      pos = (size_t)rand()%((*len) );
    }
    insert_lists.push_back({pos, ins_len});
    to_insert -= ins_len;
  }
  sort(insert_lists.begin(), insert_lists.end());
  if (do_repeat) {
    // TODO(improvement): if more than 4 inserts at the same positions,
    // some of them should be merged.
    // Not critical unless we use a very very high mutation ratio.
  }
  char * newA = new char[new_len];
  char * newB = new char[new_len];
  char * newC = new char[new_len];
  char * newD = new char[new_len];

  size_t t_pos = 0;
  size_t offset = 0;
  int prev_dice = 0;
  for (size_t i = 0; i < insert_lists.size(); i++) {
    for (; t_pos < insert_lists.at(i).first; t_pos++) {
      // copy others:
      newA[t_pos + offset] = (*A)[t_pos];
      newB[t_pos + offset] = (*B)[t_pos];
      newC[t_pos + offset] = (*C)[t_pos];
      newD[t_pos + offset] = (*D)[t_pos];
    }

    int dice = rand()%4;
    if (i > 0) {
      if (insert_lists.at(i).first == insert_lists.at(i-1).first) {
        dice = (prev_dice +1)%4;
      }
    }
    prev_dice = dice;
    char * target;
    char * source;
    if (dice == 0) {
      target = newA;
      source = *A;
    } else if (dice == 1) {
      target = newB;
      source = *B;
    } else if (dice == 2) {
      target = newC;
      source = *C;
    } else if (dice == 3) {
      target = newD;
      source = *D;
    } else {
      Debug::AbortPrint("This never happens!\n");
    }
    size_t insert_len = insert_lists.at(i).second;
    size_t tandem_length = 2 + (size_t)rand()%5;
    tandem_length = min(tandem_length, insert_len);
    assert(tandem_length <= insert_len);
    for (size_t j = 0; j < insert_len; j++) {
      newA[t_pos + offset + j] = '-';
      newB[t_pos + offset + j] = '-';
      newC[t_pos + offset + j] = '-';
      newD[t_pos + offset + j] = '-';
      if (do_repeat) {
        assert(t_pos >= insert_len);
        if (do_tandem) {
          target[t_pos + offset + j] = source[t_pos - tandem_length + (j % tandem_length)];
        } else {
          target[t_pos + offset + j] = source[t_pos - insert_len + j];
        }
      } else {
        target[t_pos + offset + j] = RandomChar();
      }
    }
    offset += insert_lists.at(i).second;
  }
  assert((*len) + offset == new_len);
  for (; t_pos < (*len); t_pos++) {
    newA[t_pos + offset] = (*A)[t_pos];
    newB[t_pos + offset] = (*B)[t_pos];
    newC[t_pos + offset] = (*C)[t_pos];
    newD[t_pos + offset] = (*D)[t_pos];
  }

  delete[] (*A);
  delete[] (*B);
  delete[] (*C);
  delete[] (*D);
  *A = newA;
  *B = newB;
  *C = newC;
  *D = newD;
  *len =new_len;
  return 0;
}


int Utils::MutateLONGDEL(char *A , size_t len, double mutation_ratio) {
  return MutateGENERIC_DEL(A, len, mutation_ratio, min((size_t)50, len/3), len/2);
}

int Utils::MutateSHORTDEL(char *A , size_t len, double mutation_ratio) {
  return MutateGENERIC_DEL(A, len, mutation_ratio, 5, 49);
}

int Utils::MutateGENERIC_DEL(char *A ,
                             size_t len,
                             double mutation_ratio,
                             size_t min_size,
                             size_t max_size) {
  size_t to_delete = len*mutation_ratio;
  size_t _to_delete = len*mutation_ratio;
  size_t real_deletes = 0;

  size_t min_indel = min(min_size, len/3);
  size_t max_indel = min(max_size, len/2);

  while (to_delete) {
    max_indel = min(max_indel, to_delete);
    min_indel = min(min_indel, to_delete);
    assert(max_indel >= min_indel);

    size_t del_size = min_indel + (size_t)rand()%(max_indel - min_indel + 1);

    size_t del_pos = (size_t)rand() % (len - del_size);
    to_delete -= del_size;
    while (del_size != 0 && del_pos != len) {
      if (A[del_pos] != '-') {
        A[del_pos]  = '-';
        del_size--;
        real_deletes++;
      }
      del_pos++;
    }
  }
  assert(real_deletes <= _to_delete);
  if (real_deletes < _to_delete) {
    printf("Aborted because we reach enf of seq.\n");  // TODO(debug): Check if synthetic trio 300 20 still trigger that. Define behavior.
  }
  return real_deletes;
}

int Utils::MutateSNP(char *A , size_t len, double mutation_ratio) {
  int err_count = 0;
  for (size_t i = 0; i < len; i++) {
    double r = rand() / (double)RAND_MAX;
    if (r < mutation_ratio && A[i] != '-') {
      char c = MutateChar(A[i]);
      A[i] = c;
      err_count++;
    }
  }
  return err_count;
}

int Utils::MutatePointDEL(char *A , size_t len, double mutation_ratio) {
  int err_count = 0;
  for (size_t i = 0; i < len; i++) {
    double r = rand() / (double)RAND_MAX;
    if (r < mutation_ratio && A[i] != '-') {
      A[i] = '-';
      err_count++;
    }
  }
  return err_count;
}

int Utils::MutateNoise(char *A , size_t len, double mutation_ratio) {
  int err_count = 0;
  for (size_t i = 0; i < len; i++) {
    double r = rand() / (double)RAND_MAX;
    if (r < mutation_ratio) {
      if (A[i] == '-') {
        A[i] = RandomChar();
      } else {
        A[i] = MutateChar(A[i]);
      }
      err_count++;
    }
  }
  return err_count;
}

char Utils::MutateChar(char a) {
  int pos = -1;
  for (int i = 0; i < 4; i++) {
    if (a == gen_alph[i]) pos = i;
  }
  if (pos == -1)
    Debug::AbortPrint("Invalid character\n");
  char tmp = gen_alph[3];
  gen_alph[3] = a;
  gen_alph[pos] = tmp;

  int new_pos = rand()%3;
  assert(gen_alph[new_pos] != a);
  return gen_alph[new_pos];
}

char Utils::RandomChar() {
  return gen_alph[rand()%4];
}

void Utils::Compact(char ** A,
                    char ** B,
                    size_t *len) {
  size_t new_len = *len;
  for (size_t i = 0; i < *len; i++) {
    if ((*A)[i] == '-' && (*B)[i] == '-') {
      new_len--;
    }
  }
  if (new_len == (*len))
    return;

  char * new_A = new char[new_len];
  char * new_B = new char[new_len];

  size_t pos = 0;
  for (size_t i = 0; i < *len; i++) {
    if ((*A)[i] != '-' || (*B)[i] != '-') {
      new_A[pos] = (*A)[i];
      new_B[pos] = (*B)[i];
      pos++;
    }
  }
  assert(pos == new_len);
  delete[] (*A);
  delete[] (*B);
  *A = new_A;
  *B = new_B;
  *len = new_len;
}

void Utils::CreateOffspring(char * M1, char * M2, size_t M_len,
                            char * F1, char * F2, size_t F_len,
                            double pm_ratio,
                            char **C1, char **C2, size_t *C_len,
                            int * non_inherited_dels, int * non_inherited_snps) {
  assert(M_len == F_len);
  *C_len = M_len;
  size_t n_cuts;
  if (M_len > 500)
    n_cuts = 1 + M_len/100;
  else
    n_cuts = 1 + M_len/40;

  // meiotic shuffle:
  *C1 = Utils::MeioticShuffle(M1, M2, M_len, n_cuts);
  *C2 = Utils::MeioticShuffle(F1, F2, F_len, n_cuts);
  int dels = 0;
  int snps = 0;
  dels += Utils::MutatePointDEL(*C1, *C_len, pm_ratio/2);
  snps += Utils::MutateSNP(*C1, *C_len, pm_ratio/2);
  dels += Utils::MutatePointDEL(*C2, *C_len, pm_ratio/2);
  snps += Utils::MutateSNP(*C2, *C_len, pm_ratio/2);
  *non_inherited_dels = dels;
  *non_inherited_snps = snps;
}

char * Utils::MeioticShuffle(char *A , char *B, size_t len, size_t n_cuts) {
  char* son = new char[len];
  size_t* cuts = Utils::GenerateCuts(len, n_cuts);
  bool father = true;
  char* ref = A;
  size_t prev_cut = 0;
  size_t next_cut;
  for (size_t i = 0; i < n_cuts; i++) {
    next_cut = cuts[i];
    ref = father? A : B;
    for (size_t pos = prev_cut; pos < next_cut; pos++) {
      son[pos] = ref[pos];
    }
    prev_cut = next_cut;
    father = !father;
  }
  for (size_t pos = prev_cut; pos < len; pos++) {
    son[pos] = ref[pos];
  }
  delete[]cuts;
  return son;
}

size_t * Utils::GenerateCuts(size_t len, size_t n_cuts) {
  assert(n_cuts <= len - 2);
  size_t * cuts = new size_t[n_cuts];
  size_t i = 0;
  while (i < n_cuts) {
    size_t pos = 1 + (size_t)rand() % (len-2);
    cuts[i] = pos;

    bool is_new = true;
    for (size_t j = 0; j < i; j++) {
      if (cuts[j] == cuts[i])
        is_new = false;
    }
    if (is_new)
      i++;
  }
  sort(cuts, cuts+n_cuts);
  return cuts;
}

void Utils::Scramble(char *A , char *B, size_t len, size_t chunk_len) {
  // bool reversed = false;
  char *F = A;
  char *M = A;
  for (size_t i = 0; i < len; i++) {
    if (i%chunk_len == 0) {
      int coin = rand() % 2;
      if (coin == 0) {
        F = A;
        M = B;
      } else {
        F = B;
        M = A;
      }
    }
    char a = F[i];
    char b = M[i];
    A[i] = a;
    B[i] = b;
  }
}

void Utils::Truncate(char ** seq, size_t new_len) {
  char * ans = new char[new_len];
  for (size_t i = 0; i < new_len; i++) {
    ans[i] = (*seq)[i];
  }
  // delete[] (*seq);
  free(*seq);
  *seq = ans;
}

char * Utils::CopySeq(char * seq, size_t len) {
  char * ans = new char[len];
  for (size_t i = 0; i < len; i++) {
    ans[i] = seq[i];
  }
  return ans;
}

void Utils::SaveChar(char * seq, size_t length, char * file_name) {
  FILE * file_seq;
  file_seq = fopen(file_name, "w");
  if (file_seq == NULL) {
    Debug::AbortPrint("SaveChar: could not open file for: %s \n", file_name);
  }
  if (fwrite(seq, sizeof(char), length, file_seq) != (length)) Debug::AbortPrint("Error in SaveChar, write\n)");
  fclose(file_seq);
}

void Utils::SaveFastaFile(char * out_file_name, char *seq, size_t  len) {
  FILE * file;
  file = fopen(out_file_name, "w");
  if (file == NULL) {
    Debug::AbortPrint("Could not open file for: %s \n", out_file_name);
  }

  const char * header = "> tmp_file plastic header";
  size_t header_len = 25;
  if (fwrite(header, sizeof(char), header_len, file) != (header_len))
    Debug::AbortPrint("Error in SaveChar, write\n)");
  if (fwrite("\n", sizeof(char), 1, file) != 1)
    Debug::AbortPrint("Error in SaveChar, write\n)");

  size_t pos = 0;
  while (pos + 70 < len) {
    if (fwrite(&(seq[pos]), sizeof(char), 70, file) != (70))
      Debug::AbortPrint("Error in SaveChar, write\n)");
    if (fwrite("\n", sizeof(char), 1, file) != 1)
      Debug::AbortPrint("Error in SaveChar, write\n)");
    pos += 70;
  }
  if (pos < len) {
    if (fwrite(&(seq[pos]), sizeof(char), len-pos, file) != (len-pos))
      Debug::AbortPrint("Error in SaveChar, write\n)");
    if (fwrite("\n", sizeof(char), 1, file) != 1)
      Debug::AbortPrint("Error in SaveChar, write\n)");
  }
  fclose(file);
}

void Utils::ReadFastaFile(char * file_name, char **ans, size_t * len) {
  FASTAFILE *ffp;
  char* seq;
  char* name;
  size_t L;
  ffp = OpenFASTA(file_name);
  while (ReadFASTA(ffp, &seq, &name, &L)) {
    printf(">len %u\n", (uint)L);
    printf(">%s\n", name);
    // printf("%s\n",  seq);
    free(name);
  }
  CloseFASTA(ffp);

  *len = L;
  *ans = seq;
}

int Utils::PhaseErrors(char * C1, char* C2, char* phase1, char* phase2, size_t len) {
  int errors = 0;
  for (size_t i = 0; i < len; i++) {
    assert(phase1[i] == '0' || phase1[i] == '1');
    assert(phase2[i] == '0' || phase2[i] == '1');
    if (phase1[i] != phase2[i] && C1[i] != C2[i])
      errors++;
  }
  return errors;
}

int Utils::CountNonGaps(char* seq, size_t len) {
  int ans = 0;
  for (size_t i = 0; i < len; i++) {
    if (seq[i] != '-')
      ans++;
  }
  return ans;
}



