import ConfigParser

config = ConfigParser.RawConfigParser()

# When adding sections or items, add them in the reverse order of
# how you want them to be displayed in the actual file.
# In addition, please note that using RawConfigParser's and the raw
# mode of ConfigParser's respective set functions, you can assign
# non-string values to keys internally, but will receive an error
# when attempting to write to a file or when you get it in non-raw
# mode. SafeConfigParser does not allow such assignments to take place.
config.add_section('Common')
config.set('Common', 'reference', 'build37-chr21.fa')
config.set('Common', 'window', 'NA12880-chr21.posinfo')
config.add_section('Child')
config.set('Child', 'vcf_file', 'NA12880-chr21-unphased-clean.vcf')
config.set('Child', 'fasta1', 'NA12880-chr21_1.fa')
config.set('Child', 'fasta2', 'NA12880-chr21_2.fa')
config.add_section('Mother')
config.set('Mother', 'vcf_file', 'NA12878-chr21-unphased-clean.vcf')
config.set('Mother', 'fasta1', 'NA12878-chr21_1.fa')
config.set('Mother', 'fasta2', 'NA12878-chr21_2.fa')
config.add_section('Father')
config.set('Father', 'vcf_file', 'NA12877-chr21-unphased-clean.vcf')
config.set('Father', 'fasta1', 'NA12877-chr21_1.fa')
config.set('Father', 'fasta2', 'NA12877-chr21_2.fa')

# Writing our configuration file to 'example.cfg'
with open('family.config', 'wb') as configfile:
    config.write(configfile)

