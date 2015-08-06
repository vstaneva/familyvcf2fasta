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
  fprintf(stderr, "./mfc_similarity_phaser fatherA.fa fatherB.fa motherA.fa motherB.fa childA.fa childB.fa \n");  // NOLINT
}

int main(int argc, char ** argv) {
  if (argc != 7) {
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

  char * phase_string_1;
  char * phase_string_2;
  score_t score_1;
  score_t score_2;

  Utils::StartClock();
  {
    Phaser * phaser =  new Phaser(motherA, motherB, mother_len,
                                  fatherA, fatherB, father_len,
                                  childA, childB, child_len);
    phaser->SetScoreGap(SCORE_GAP);
    phaser->SetScoreMismatch(SCORE_MISMATCH);
    phaser->SetScoreMatch(SCORE_MATCH);
    score_1 = phaser->similarity_and_phase();
    phase_string_1 = phaser->GetPhaseString();
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
  }
  
  double time = Utils::StopClock();
  assert(score_1 == score_2);
  
  char * consensus = new char[child_len];
  for (size_t i = 0; i < child_len; i++) {
    if (phase_string_1[i] == '0' && phase_string_2[i] == '1') {
      consensus[i] = '0';
    }
    else if (phase_string_1[i] == '1' && phase_string_2[i] == '0') {
      consensus[i] = '1';
    } else {
      consensus[i] = '?';
    }
  }

  const char * output_filename = "phase_string.txt";
  Utils::SaveChar(consensus, child_len, (char *)output_filename);
  printf("Similarity score: %i\n", score_1);
  printf("Took in: %.2f seconds\n", time);

  delete[] motherA;
  delete[] motherB;
  delete[] fatherA;
  delete[] fatherB;
  delete[] childA;
  delete[] childB;
}

