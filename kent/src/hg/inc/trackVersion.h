/* trackVersion.h was originally generated by the autoSql program, which also 
 * generated trackVersion.c and trackVersion.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef TRACKVERSION_H
#define TRACKVERSION_H

#ifndef JKSQL_H
#include "jksql.h"
#endif

#define TRACKVERSION_NUM_COLS 9

struct trackVersion
/* version information for database tables to monitor data loading history */
    {
    struct trackVersion *next;  /* Next in singly linked list. */
    int ix;	/* auto-increment ID */
    char *db;	/* UCSC database name */
    char *name;	/* table name in database */
    char *who;	/* Unix userID that performed this update */
    char *version;	/* version string, whatever is meaningful for data source */
    char *updateTime;	/* YYYY-MM-DD HH:MM:SS most-recent-update time */
    char *comment;	/* other comments about version */
    char *source;	/* perhaps a URL for the data source */
    char *dateReference;	/* Ensembl date string for archive reference */
    };

void trackVersionStaticLoad(char **row, struct trackVersion *ret);
/* Load a row from trackVersion table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct trackVersion *trackVersionLoad(char **row);
/* Load a trackVersion from row fetched with select * from trackVersion
 * from database.  Dispose of this with trackVersionFree(). */

struct trackVersion *trackVersionLoadAll(char *fileName);
/* Load all trackVersion from whitespace-separated file.
 * Dispose of this with trackVersionFreeList(). */

struct trackVersion *trackVersionLoadAllByChar(char *fileName, char chopper);
/* Load all trackVersion from chopper separated file.
 * Dispose of this with trackVersionFreeList(). */

#define trackVersionLoadAllByTab(a) trackVersionLoadAllByChar(a, '\t');
/* Load all trackVersion from tab separated file.
 * Dispose of this with trackVersionFreeList(). */

struct trackVersion *trackVersionLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all trackVersion from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with trackVersionFreeList(). */

void trackVersionSaveToDb(struct sqlConnection *conn, struct trackVersion *el, char *tableName, int updateSize);
/* Save trackVersion as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Strings are automatically escaped to allow insertion into the database. */

struct trackVersion *trackVersionCommaIn(char **pS, struct trackVersion *ret);
/* Create a trackVersion out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new trackVersion */

void trackVersionFree(struct trackVersion **pEl);
/* Free a single dynamically allocated trackVersion such as created
 * with trackVersionLoad(). */

void trackVersionFreeList(struct trackVersion **pList);
/* Free a list of dynamically allocated trackVersion's */

void trackVersionOutput(struct trackVersion *el, FILE *f, char sep, char lastSep);
/* Print out trackVersion.  Separate fields with sep. Follow last field with lastSep. */

#define trackVersionTabOut(el,f) trackVersionOutput(el,f,'\t','\n');
/* Print out trackVersion as a line in a tab-separated file. */

#define trackVersionCommaOut(el,f) trackVersionOutput(el,f,',',',');
/* Print out trackVersion as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

struct trackVersion *getTrackVersion(char *database, char *track);
// Get most recent trackVersion for given track in given database

#endif /* TRACKVERSION_H */
