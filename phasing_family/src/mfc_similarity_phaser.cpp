/* Copyright (C) 2013, Daniel Valenzuela, all rights reserved.
 * dvalenzu@cs.helsinki.fi
 */


#include <iomanip>
#include "./phaser.h"
#include "./debug.h"
#include "./utils.h"

score_t SCORE_GAP = -1;
score_t SCORE_MISMATCH = -1;
score_t SCORE_MATCH = 1;

bool verbose = false;
void printUssage();
void printUssage() {
  fprintf(stderr, "Ussage:\n");
  fprintf(stderr, "./mfc_similarity_phaser fatherA.fa fatherB.fa motherA.fa motherB.fa childA.fa childB.fa n_paths\n");  // NOLINT
}

void negate(char * phase_str, size_t len);
void negate(char * phase_str, size_t len) {
  for (size_t i = 0; i < len; i++) {
    if (phase_str[i] == '0')
      phase_str[i] ='1';
    else if (phase_str[i] == '1')
      phase_str[i] ='0';
    else 
      assert(0);
  }
}

char * MultiPassPhaser(char * motherA,
                       char * motherB,
                       size_t  mother_len,
                       char * fatherA,
                       char * fatherB,
                       size_t father_len,
                       char * childA,
                       char * childB,
                       size_t child_len,
                       int n_paths,
                       score_t * score_ans);

char * MultiPassPhaser(char * motherA,
                       char * motherB,
                       size_t  mother_len,
                       char * fatherA,
                       char * fatherB,
                       size_t father_len,
                       char * childA,
                       char * childB,
                       size_t child_len,
                       int n_paths,
                       score_t * score_ans) {
  char * phase_string_1;
  char * phase_string_2;
  char * phase_string_3;
  char * phase_string_4;
  score_t score_1;
  score_t score_2;
  score_t score_3;
  score_t score_4;
  {
    Phaser * phaser =  new Phaser(motherA, motherB, mother_len,
                                  fatherA, fatherB, father_len,
                                  childA, childB, child_len);
    phaser->SetScoreGap(SCORE_GAP);
    phaser->SetScoreMismatch(SCORE_MISMATCH);
    phaser->SetScoreMatch(SCORE_MATCH);
    score_1 = phaser->similarity_and_phase();
    phase_string_1 = phaser->GetPhaseString();
    if (n_paths == 1) 
      return phase_string_1;
  }

  {
    // changed order of M and F:
    Phaser * phaser =  new Phaser(fatherA, fatherB, father_len,
                                  motherA, motherB, mother_len,
                                  childA, childB, child_len);
    phaser->SetScoreGap(SCORE_GAP);
    phaser->SetScoreMismatch(SCORE_MISMATCH);
    phaser->SetScoreMatch(SCORE_MATCH);
    score_2 = phaser->similarity_and_phase();
    phase_string_2 = phaser->GetPhaseString();
    negate(phase_string_2, child_len);
  }
  if (score_1 != score_2) {
    fprintf(stderr,"WARNING: Different scores after changing the order of params, this should not occur\n");
  }

  if (n_paths == 4) {
    {
      Phaser * phaser =  new Phaser(motherA, motherB, mother_len,
                                    fatherA, fatherB, father_len,
                                    childB, childA, child_len);
      phaser->SetScoreGap(SCORE_GAP);
      phaser->SetScoreMismatch(SCORE_MISMATCH);
      phaser->SetScoreMatch(SCORE_MATCH);
      score_3 = phaser->similarity_and_phase();
      phase_string_3 = phaser->GetPhaseString();
      negate(phase_string_3, child_len);
    }
    if (score_1 != score_3) {
      fprintf(stderr,"WARNING: Different scores after changing the order of params, this should not occur\n");
    }

    {
      // changed order of M and F:
      Phaser * phaser =  new Phaser(fatherA, fatherB, father_len,
                                    motherA, motherB, mother_len,
                                    childB, childA, child_len);
      phaser->SetScoreGap(SCORE_GAP);
      phaser->SetScoreMismatch(SCORE_MISMATCH);
      phaser->SetScoreMatch(SCORE_MATCH);
      score_4 = phaser->similarity_and_phase();
      phase_string_4 = phaser->GetPhaseString();
    }
    if (score_1 != score_4) {
      fprintf(stderr,"WARNING: Different scores after changing the order of params, this should not occur\n");
    }
  } 
  
  int * sum = new int[child_len];
  
  char * consensus = new char[child_len];
  for (size_t i = 0; i < child_len; i++) {
    sum[i] = 0;
    if (phase_string_1[i] == '1') sum[i]++; 
    if (phase_string_2[i] == '1') sum[i]++; 
    if (n_paths == 4) {
      if (phase_string_3[i] == '1') sum[i]++; 
      if (phase_string_4[i] == '1') sum[i]++; 
    }
    //if (sum  == 3)
    //  std::cout << "sum:" << sum[i] << std::endl;
    
    int mid = 1;
    if (n_paths ==4 ) {
      mid = 2;
    }
    if (sum[i] < mid) {
      consensus[i] = '0';
    }
    else if (sum[i] > mid) {
      consensus[i] = '1';
    } else {
      assert(sum[i] == mid);
      consensus[i] = '?';
    }
  }
  *score_ans = score_1;
  return consensus;
}


int main(int argc, char ** argv) {
  if (argc != 8) {
    printUssage();
    return EXIT_FAILURE;
  }
  char * motherA;
  char * motherB;
  size_t mother_len;

  char * fatherA;
  char * fatherB;
  size_t father_len;

  char * childA;
  char * childB;
  size_t child_len;

  size_t seq_len;
  // mother:
  Utils::ReadFastaFile(argv[1], &motherA, &seq_len);
  mother_len = seq_len;

  Utils::ReadFastaFile(argv[2], &motherB, &seq_len);
  if (mother_len != seq_len)
    Debug::AbortPrint("motherA and motherB have different length. They must be an alignment.\n");

  // father:
  Utils::ReadFastaFile(argv[3], &fatherA, &seq_len);
  father_len = seq_len;

  Utils::ReadFastaFile(argv[4], &fatherB, &seq_len);
  if (father_len != seq_len)
    Debug::AbortPrint("fatherA and fatherB have different length. They must be an alignment.\n");

  // child:
  Utils::ReadFastaFile(argv[5], &childA, &seq_len);
  child_len = seq_len;

  Utils::ReadFastaFile(argv[6], &childB, &seq_len);
  if (child_len != seq_len)
    Debug::AbortPrint("childA and childB have different length. They must be an alignment.\n");

  if (!(argv[7][0] == '1' || argv[7][0] == '2' || argv[7][0] == '4')) {
    std::cout << " n_paths must be 1, 2 or 4" << std::endl;
    return 33;
  }
  
  int n_paths = atoi(argv[7]);
  std::cout << n_paths << " paths will be used" << std::endl;

  Utils::StartClock();
  score_t score;
  char * consensus = MultiPassPhaser(motherA, motherB, mother_len,
                                     fatherA, fatherB, father_len,
                                     childA, childB, child_len, 
                                     n_paths, &score);
  double time = Utils::StopClock();

  const char * output_filename = "phase_string.txt";
  Utils::SaveChar(consensus, child_len, (char *)output_filename);
  printf("Similarity score: %i\n", score);
  printf("Took in: %.2f seconds\n", time);

  delete[] motherA;
  delete[] motherB;
  delete[] fatherA;
  delete[] fatherB;
  delete[] childA;
  delete[] childB;
}

