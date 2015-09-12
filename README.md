# vcf2fasta

[![Build Status](https://travis-ci.org/vstaneva/familyvcf2fasta.svg?branch=master)](https://travis-ci.org/vstaneva/familyvcf2fasta)

-- writeconfig.py creates a config file, named family.config

-- mfcVCFtoFASTA.py takes one argument, a configuration file, and produces FASTA files that can be analyzed with [phasing_family]


Known Issues
* We rely on vcf unique identifiers, if they are not present, the behavior is undefined.
