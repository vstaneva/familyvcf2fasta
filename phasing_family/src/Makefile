CPP=g++
#CPP=g++-4.7
PARANOID=-pedantic -Wall -Wextra -Wcast-align -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Winline -Wno-error=unused-parameter -Wno-error=unused-variable

CPPFLAGS=-std=c++11 -ggdb -O0 -Wall -pedantic -Wunused-parameter $(PARANOID)
#CPPFLAGS=-std=c++11 -DNDEBUG -O3 -Wall -pedantic -Wunused-parameter $(PARANOID)

LIB_OBJECTS=phaser.o utils.o fasta.o
BIN_OBJECTS=test_phaser.o synthetic_trio.o mfc_similarity_phaser.o
OBJECTS=$(LIB_OBJECTS) $(BIN_OBJECTS)
BIN=test_phaser synthetic_trio mfc_similarity_phaser

LIB=$(LIB_OBJECTS)

%.o: %.cpp
	@echo " [$(CPP)] Compiling $<"
	@$(CPP) $(CPPFLAGS) $(INCS) -c $< -o $@

all: $(OBJECTS) $(BIN)
	@echo " [MSG] Done compiling"
	@echo " [FLG] $(CPPFLAGS)"

test_phaser: test_phaser.cpp $(OBJECTS)
	@echo " [LNK] Building test_phaser"
	@$(CPP) $(CPPFLAGS) -o test_phaser test_phaser.o $(LIB) 

synthetic_trio: synthetic_trio.cpp $(OBJECTS)
	@echo " [LNK] Building synthetic_trio"
	@$(CPP) $(CPPFLAGS) -o synthetic_trio synthetic_trio.o $(LIB) 

mfc_similarity_phaser: mfc_similarity_phaser.cpp $(OBJECTS)
	@echo " [LNK] Building mfc_similarity_phaser"
	@$(CPP) $(CPPFLAGS) -o mfc_similarity_phaser mfc_similarity_phaser.o $(LIB) 

clean:
	@echo " [CLN] Cleaning object, binary files."
	@rm -f $(OBJECTS) $(BIN); rm -f *.tmp; rm -f tmp.*;

types: types.vim
types.vim: *.[ch]
	ctags --c-kinds=gstu -o- *.[ch] |\
		awk 'BEGIN{printf("syntax keyword Type\t")}\
			{printf("%s ", $$1)}END{print ""}' > $@
