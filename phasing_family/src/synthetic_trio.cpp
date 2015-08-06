/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */


#include <iomanip>
#include "./phaser.h"
#include "./debug.h"
#include "./utils.h"

using std::cout;
using std::left;
using std::setw;
using std::endl;
using std::setfill;

score_t SCORE_GAP = -1;
score_t SCORE_MISMATCH = -1;
score_t SCORE_MATCH = 1;

bool verbose = false;
void Sensibility(char * sequence, size_t seq_len, int n_repeats);
void Visual(char * sequence, size_t seq_len, int n_repeats);
void Table(char * sequence, size_t seq_len, int n_repeats);
void Experiment(char * sequence,
                size_t seq_len,
                int n_repeats,
                double error_ratio,
                double long_indels,
                double short_indels,
                double repeat_ratio,
                double str_ratio,
                double parents_pm_ratio,
                double child_pm_ratio,
                bool visual);

void PrintTableHead();
void PrintTableRow(
    size_t length,
    double avg_C_length,
    double time,
    double score,
    double phase_errors,
    double error_ratio,
    double parents_pm_ratio,
    double parents_long_indel_ratio,
    double parents_short_indel_ratio,
    double parents_repeat_ratio,
    double parents_str_ratio,
    double child_pm_ratio);

void PrintDataFile(size_t length,
                   double avg_phase_errors,
                   double error_ratio,
                   double parents_pm_ratio,
                   double parents_long_indel_ratio,
                   double parents_short_indel_ratio,
                   double parents_repeat_ratio,
                   double parents_str_ratio,
                   double child_pm_ratio);

void printUssage();
void printUssage() {
  fprintf(stderr, "Ussage:\n");
  fprintf(stderr, "./phase_synthetic_trio seed.fa length n_repeats\n");
}

int main(int argc, char ** argv) {
  // prefer gap over mismatch:
  assert(SCORE_MISMATCH <= SCORE_GAP);
  // prefer mismatch over two gaps:
  assert(2*SCORE_GAP <= SCORE_MISMATCH);
  if (argc != 4) {
    printUssage();
    return EXIT_FAILURE;
  }
  char* file_name = argv[1];
  size_t length = (size_t)atol(argv[2]);
  int n_repeats = atol(argv[3]);
  char * sequence;
  size_t seq_len;
  Utils::ReadFastaFile(file_name, &sequence, &seq_len);
  if (seq_len > length) {
    Utils::Truncate(&sequence, length);
    seq_len = length;
  } else if (seq_len  < length) {
    Debug::AbortPrint("Sequence is shorter than asked length.\n");
  }

  Sensibility(sequence, seq_len, n_repeats);
  // Visual(sequence, seq_len, n_repeats);
  // Table(sequence, seq_len, n_repeats);

  printf("Experiments Successful\n");
  delete[] sequence;
}


void Sensibility(char * sequence, size_t seq_len, int n_repeats) {
  PrintTableHead();
  bool visual = false;
  double def_value = 0;
  double default_ratios[8] = {def_value ,
    def_value,
    def_value,
    def_value,
    def_value,
    def_value,
    def_value};
  const char * names[8] = {
    "nothing    :\n",
    "gral noise :\n",
    "long indels:\n",
    "short indels:\n",
    "repeats    :\n",
    "strs    :\n",
    "parents pm :\n",
    "child pm   :\n"};
  for (double error_ratio : {0.0, 0.05, 0.1}) {
    for (size_t col = 2; col < 8; col++) {
      for (double curr_ratio : {0.01, 0.03, 0.06, 0.09, 0.12, 0.15}) {
        default_ratios[col] = curr_ratio;
        if (verbose) {
          printf("\n------\n%s", names[col]);
        }
        Experiment(sequence,
                   seq_len,
                   n_repeats,
                   error_ratio,
                   default_ratios[2],
                   default_ratios[3],
                   default_ratios[4],
                   default_ratios[5],
                   default_ratios[6],
                   default_ratios[7],
                   visual);
        default_ratios[col] = def_value;
      }
    }
    Experiment(sequence,
               seq_len,
               n_repeats,
               error_ratio,
               default_ratios[2],
               default_ratios[3],
               default_ratios[4],
               default_ratios[5],
               default_ratios[6],
               default_ratios[7],
               visual);
  }
}

// "all against all"
void Table(char * sequence, size_t seq_len, int n_repeats) {
  PrintTableHead();
  bool visual = false;
  for (double error_ratio : {0.0, 0.1}) {
    for (double long_indels : {0.0, 0.1}) {
      for (double short_indels : {0.0, 0.1}) {
        for (double repeats : {0.0, 0.1}) {
          for (double strs : {0.0, 0.1}) {
            for (double child_point_mutation : {0.0, 0.1 }) {
              for (double parents_point_mutation : {0.0, 0.1 }) {
                Experiment(sequence,
                           seq_len,
                           n_repeats,
                           error_ratio,
                           long_indels,
                           short_indels,
                           repeats,
                           strs,
                           parents_point_mutation,
                           child_point_mutation,
                           visual);
              }
              // printf("\\hline\n");
            }
          }
        }
      }
    }
  }
}

void Visual(char * sequence, size_t seq_len, int n_repeats) {
  PrintTableHead();
  bool visual = true;
  size_t ratios_len = 8;
  double ratios[8] = {0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const char * names[8] = {
    "nothing    :\n",
    "gral noise :\n",
    "long indels:\n",
    "short indels:\n",
    "repeats    :\n",
    "strs    :\n",
    "parents pm :\n",
    "child pm   :\n"};
  for (size_t i = 0; i < ratios_len; i++) {
    ratios[i] = 0.4;
    printf("\n------\n%s", names[i]);
    Experiment(sequence,
               seq_len,
               n_repeats,
               ratios[1],
               ratios[2],
               ratios[3],
               ratios[4],
               ratios[5],
               ratios[6],
               ratios[7],
               visual);
    ratios[i] = 0.0;
  }
}


void Experiment(char * sequence,
                size_t seq_len,
                int n_repeats,
                double error_ratio,
                double long_indels,
                double short_indels,
                double repeat_ratio,
                double str_ratio,
                double parents_pm_ratio,
                double child_pm_ratio,
                bool visual) {
  int min_error = INT_MAX;
  int max_error = INT_MIN;
  int tot_error = 0;
  score_t tot_score = 0;
  double tot_time = 0;
  size_t tot_C_len = 0;
  for (int r = 0; r < n_repeats; r++) {
    char *M1, *M2, *F1, *F2, *C1, *C2;
    size_t M_len, F_len, C_len;
    Utils::CreateTrio(sequence,
                      seq_len,
                      error_ratio,
                      long_indels,
                      short_indels,
                      repeat_ratio,
                      str_ratio,
                      parents_pm_ratio,
                      child_pm_ratio,
                      &M1, &M2, &M_len,
                      &F1, &F2, &F_len,
                      &C1, &C2, &C_len);

    char * phase_real = new char[C_len];
    for (size_t i = 0; i < C_len; i++) phase_real[i] = '0';
    // Scramble child? shouldnt change anything.

    Phaser * phaser =  new Phaser(M1, M2, M_len,
                                  F1, F2, F_len,
                                  C1, C2, C_len);
    phaser->SetScoreGap(SCORE_GAP);
    phaser->SetScoreMismatch(SCORE_MISMATCH);
    phaser->SetScoreMatch(SCORE_MATCH);
    Utils::StartClock();
    score_t score = phaser->similarity_and_phase();
    double time = Utils::StopClock();
    char * phase_algor = phaser->GetPhaseString();
    int phase_errors = Utils::PhaseErrors(C1, C2, phase_real, phase_algor, C_len);
    //
    min_error = (phase_errors < min_error) ? phase_errors: min_error;
    max_error = (phase_errors > max_error) ? phase_errors: max_error;
    tot_error += phase_errors;
    tot_score += score;
    tot_time += time;
    tot_C_len += C_len;
    //
    if (visual) {
      phaser->PrintSequences();
    }
    delete(phaser);
    delete[]M1;
    delete[]M2;
    delete[]F1;
    delete[]F2;
    delete[]C1;
    delete[]C2;
    delete[]phase_real;
  }

  // TODO(code quality): Incosistent way of coputing normalized phasing errors.
  // PrintTable computes that itself, to DataFile we provide it computed.

  PrintTableRow(seq_len,
                (double)tot_C_len/(double)n_repeats,
                tot_time/(double)n_repeats,
                (double)tot_score/(double)n_repeats,
                (double)tot_error/(double)n_repeats,
                error_ratio,
                parents_pm_ratio,
                long_indels,
                short_indels,
                repeat_ratio,
                str_ratio,
                child_pm_ratio);
  if (!visual) {
    PrintDataFile(seq_len,
                  (double)tot_error/(double)tot_C_len,
                  error_ratio,
                  parents_pm_ratio,
                  long_indels,
                  short_indels,
                  repeat_ratio,
                  str_ratio,
                  child_pm_ratio);
  }
}

void PrintTableHead() {
  const char sep    = ' ';
  const int width   = 10;
  const int l_width   = 15;

  cout << left << setw(width) << setfill(sep)  << "Length" << "|";
  cout << left << setw(width) << setfill(sep)  << "Len (avg)" << "|";
  cout << left << setw(width) << setfill(sep)  << "Time" << "|";
  cout << left << setw(width) << setfill(sep)  << "Score" << "|";
  cout << left << setw(width) << setfill(sep)  << "Ph. errors" << "|";
  cout << left << setw(l_width) << setfill(sep)  << "Ph. errors (/)" << "||";

  cout << left << setw(width) << setfill(sep)  << "PPMs (/)" << "|";
  cout << left << setw(width) << setfill(sep)  << "LIs (/)" << "|";
  cout << left << setw(width) << setfill(sep)  << "SIs (/)" << "|";
  cout << left << setw(width) << setfill(sep)  << "REAP (/)" << "|";
  cout << left << setw(width) << setfill(sep)  << "STRs (/)" << "|";
  cout << left << setw(width) << setfill(sep)  << "cPMs (/)" << "|";
  cout << left << setw(width) << setfill(sep)  << "NOISE" << endl;
}

void PrintTableRow(size_t length,
                   double avg_C_length,
                   double avg_time,
                   double avg_score,
                   double avg_phase_errors,
                   double error_ratio,
                   double parents_pm_ratio,
                   double parents_long_indel_ratio,
                   double parents_short_indel_ratio,
                   double parents_repeat_ratio,
                   double parents_str_ratio,
                   double child_pm_ratio) {
  const char sep    = ' ';
  const int width   = 10;
  const int l_width   = 15;

  cout << left << setw(width) << setfill(sep)  << length << "|";
  cout << left << setw(width) << setfill(sep)  << avg_C_length << "|";
  cout << left << setw(width) << setfill(sep)  << avg_time << "|";
  cout << left << setw(width) << setfill(sep)  << avg_score << "|";
  cout << left << setw(width) << setfill(sep)  << avg_phase_errors << "|";
  cout << left << setw(l_width) << setfill(sep)  << avg_phase_errors/avg_C_length << "||";
  cout << left << setw(width) << setfill(sep)  << parents_pm_ratio << "|";
  cout << left << setw(width) << setfill(sep)  << parents_long_indel_ratio << "|";
  cout << left << setw(width) << setfill(sep)  << parents_short_indel_ratio << "|";
  cout << left << setw(width) << setfill(sep)  << parents_repeat_ratio << "|";
  cout << left << setw(width) << setfill(sep)  << parents_str_ratio << "|";
  cout << left << setw(width) << setfill(sep)  << child_pm_ratio << "|";
  cout << left << setw(width) << setfill(sep)  << error_ratio << endl;
}

void PrintDataFile(size_t length,
                   double avg_phase_errors,
                   double error_ratio,
                   double parents_pm_ratio,
                   double parents_long_indel_ratio,
                   double parents_short_indel_ratio,
                   double parents_repeat_ratio,
                   double parents_str_ratio,
                   double child_pm_ratio) {
  char fname[1024];
  sprintf (fname,"LEN-%lu-NOISE-%.2f-PARENTS_LONG_INDEL-%.2f-PSI-%.2f-PR-%.2f-PSTR-%.2f-PARENTS_PM-%.2f-CHILD_PM-%.2f.data",  // NOLINT
           length,
           error_ratio,
           parents_long_indel_ratio,
           parents_short_indel_ratio,
           parents_repeat_ratio,
           parents_str_ratio,
           parents_pm_ratio,
           child_pm_ratio);
  FILE *fp;
  fp = fopen(fname, "w");
  fprintf(fp, "%lu %.8f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
           length,
           avg_phase_errors,
           error_ratio,
           parents_long_indel_ratio,
           parents_short_indel_ratio,
           parents_repeat_ratio,
           parents_str_ratio,
           parents_pm_ratio,
           child_pm_ratio);

  fclose(fp);
}

