/* stsInfoMouseNew.h was originally generated by the autoSql program, which also 
 * generated stsInfoMouseNew.c and stsInfoMouseNew.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2003 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef STSINFOMOUSENEW_H
#define STSINFOMOUSENEW_H

#define STSINFOMOUSENEW_NUM_COLS 25

struct stsInfoMouseNew
/* Constant STS marker information */
    {
    struct stsInfoMouseNew *next;  /* Next in singly linked list. */
    unsigned identNo;	/* UCSC identification number */
    char *name;	/* Official UCSC name */
    unsigned MGIId;	/* Marker's MGI Id */
    char *MGIName;	/* Marker's MGI name */
    unsigned UiStsId;	/* Marker's UiStsId */
    unsigned nameCount;	/* Number of alias */
    char *alias;	/* alias, or N/A */
    char *primer1;	/* primer1 sequence */
    char *primer2;	/* primer2 sequence */
    char *distance;	/* Length of STS sequence */
    unsigned sequence;	/* Whether the full sequence is available (1) or not (0) for STS */
    char *organis;	/* Organism for which STS discovered */
    char *wigName;	/* WI_Mouse_Genetic map */
    char *wigChr;	/* Chromosome in Genetic map */
    float wigGeneticPos;	/* Position in Genetic map */
    char *mgiName;	/* MGI map */
    char *mgiChr;	/* Chromosome in Genetic map */
    float mgiGeneticPos;	/* Position in Genetic map */
    char *rhName;	/* WhiteHead_RH map */
    char *rhChr;	/* Chromosome in Genetic map */
    float rhGeneticPos;	/* Position in Genetic map. */
    float RHLOD;	/* LOD score of RH map */
    char *GeneName;	/* Associated gene name */
    char *GeneID;	/* Associated gene Id */
    char *clone;	/* Clone sequence */
    };

void stsInfoMouseNewStaticLoad(char **row, struct stsInfoMouseNew *ret);
/* Load a row from stsInfoMouseNew table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct stsInfoMouseNew *stsInfoMouseNewLoad(char **row);
/* Load a stsInfoMouseNew from row fetched with select * from stsInfoMouseNew
 * from database.  Dispose of this with stsInfoMouseNewFree(). */

struct stsInfoMouseNew *stsInfoMouseNewLoadAll(char *fileName);
/* Load all stsInfoMouseNew from whitespace-separated file.
 * Dispose of this with stsInfoMouseNewFreeList(). */

struct stsInfoMouseNew *stsInfoMouseNewLoadAllByChar(char *fileName, char chopper);
/* Load all stsInfoMouseNew from chopper separated file.
 * Dispose of this with stsInfoMouseNewFreeList(). */

#define stsInfoMouseNewLoadAllByTab(a) stsInfoMouseNewLoadAllByChar(a, '\t');
/* Load all stsInfoMouseNew from tab separated file.
 * Dispose of this with stsInfoMouseNewFreeList(). */

struct stsInfoMouseNew *stsInfoMouseNewCommaIn(char **pS, struct stsInfoMouseNew *ret);
/* Create a stsInfoMouseNew out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new stsInfoMouseNew */

void stsInfoMouseNewFree(struct stsInfoMouseNew **pEl);
/* Free a single dynamically allocated stsInfoMouseNew such as created
 * with stsInfoMouseNewLoad(). */

void stsInfoMouseNewFreeList(struct stsInfoMouseNew **pList);
/* Free a list of dynamically allocated stsInfoMouseNew's */

void stsInfoMouseNewOutput(struct stsInfoMouseNew *el, FILE *f, char sep, char lastSep);
/* Print out stsInfoMouseNew.  Separate fields with sep. Follow last field with lastSep. */

#define stsInfoMouseNewTabOut(el,f) stsInfoMouseNewOutput(el,f,'\t','\n');
/* Print out stsInfoMouseNew as a line in a tab-separated file. */

#define stsInfoMouseNewCommaOut(el,f) stsInfoMouseNewOutput(el,f,',',',');
/* Print out stsInfoMouseNew as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* STSINFOMOUSENEW_H */
