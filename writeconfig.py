import ConfigParser

config = ConfigParser.RawConfigParser()

# When adding sections or items, add them in the reverse order of
# how you want them to be displayed in the actual file.

#The common section includes a FASTA reference file and a .posinfo file
#The .posinfo file has three lines: the chromosome, and the start and end
#of the region that will be phased. Interesting regions can be obtained
#by running help_scripts/findInterestingRegions.py
#The Child, Mother and Father sections each contain a VCF file
#for that individual, as well as paths to the two aligned FASTA files that
#contain the non-overlapping variants,
#as described in [Makinen and Valenzuela, 2015]
config.add_section('Common')
config.set('Common', 'reference', 'build37-chr21.fa')
config.set('Common', 'window', 'NA12880-chr21.posinfo')
config.add_section('Child')
config.set('Child', 'vcf_file', 'NA12880-chr21-unphased-clean.vcf')
config.set('Child', 'fasta1', 'NA12880-chr21_1.fa')
config.set('Child', 'fasta2', 'NA12880-chr21_2.fa')
config.add_section('Mother')
config.set('Mother', 'vcf_file', 'NA12878-chr21-clean.vcf')
config.set('Mother', 'fasta1', 'NA12878-chr21_1.fa')
config.set('Mother', 'fasta2', 'NA12878-chr21_2.fa')
config.add_section('Father')
config.set('Father', 'vcf_file', 'NA12877-chr21-clean.vcf')
config.set('Father', 'fasta1', 'NA12877-chr21_1.fa')
config.set('Father', 'fasta2', 'NA12877-chr21_2.fa')

# Writing our configuration file to 'family.config'
with open('family.config', 'wb') as configfile:
    config.write(configfile)

