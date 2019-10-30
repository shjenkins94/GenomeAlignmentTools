# dbSnpRs.sql was originally generated by the autoSql program, which also 
# generated dbSnpRs.c and dbSnpRs.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Information from dbSNP at the reference SNP level
CREATE TABLE dbSnpRsRn (
    rsId      varchar(255) not null, # dbSnp reference snp (rs) identifier
    avHet     float        not null, # the average heterozygosity from all observations
    avHetSE   float        not null, # the Standard Error for the average heterozygosity
    valid     set (
		'no-information',
		'by-2hit-2allele',
		'by-cluster',
		'by-frequency',
		'other-pop'
		) 	   not null default 'no-information', # the validation status of the SNP
    allele1   varchar(255) not null, # the sequence of the first allele
    allele2   varchar(255) not null, # the sequence of the second allele
    assembly  varchar(255) not null, # the sequence in the assembly
    alternate varchar(255) not null, # the sequence of the alternate allele
    func      set (
		'coding-synon',    # synonymous mutation
		'coding-nonsynon', # non-synonymous mutation
		'mrna-utr',        # untranslated region
		'intron',          # intronic region
		'splice-site',     # splice site
		'exception',       # alignment gap
		'coding',          # coding region
		'reference',       # 
		'locus-region',    # locus region
		'X'                # no functional classification
		) 	   not null default 'X', # the functional category of the SNP, if any
              #Indexes
    PRIMARY KEY(rsId)
);