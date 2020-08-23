/*
 *  Test.c
 *
 *    Extracts the comments in CMISS and places them in a series of
 *    "database" files for later use.
 *
 *    Excerpts from the Tao of Programming: Book Four, Verse 1
 *
 *    A program should be light and agile, its subroutines connected like a
 *    strings of pearls. The spirit and intent of the program should be
 *    retained throughout. There should be neither too little nor too much,
 *    neither needless loops nor useless variables, neither lack of structure
 *    nor overwhelming rigidity.
 *    
 *    A program should follow the 'Law of Least Astonishment'. What is this
 *    law? It is simply that the program should always respond to the user
 *    in the way that astonishes him least.
 *    
 *    A program, no matter how complex, should act as a single unit. The
 *    program should be directed by the logic within rather than by outward
 *    appearances.
 *    
 *    If the program fails in these requirements, it will be in a state of
 *    disorder and confusion. The only way to correct this is to rewrite the
 *    program.
 *    
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: malloc(), free(), NULL */

#include <errno.h>

#include <stdlib.h>
    /* For: exit() */

#include <string.h>
    /* For: strcmp() */

#include <ctype.h>
    /* For: isprint(), tolower() */

#include "cmiss-lineio.h"
    /* For: FILESTATE */

#include "outputcommands.h"
    /* For  outputCommand(), outputCommand2() etc */

#include "wwwpaths.h"
    /* For CMISS_WWW_ROOT LOOKUP_URLPATH LOOKUP_PROGRAM_URLPATH */


#include "tags.h"
    /* For: TAG, getTag() */

#include "groups.h"
    /* For: GROUP */

#include "rules.h"
    /* For: RULES, newRules(), etc */
 
#include "rule-parse.h"
    /* For: parseRuleFile() etc */

#include "cgi-decode.h"
  /* For: CGI handling functions */

#include "htstring.h"
    /* For: htmlEscapedStrcpy(), htmlEscapedStrlen() */

#include "utilities.h"
    /* For: outputString() */

#include "output.h"
    /* For: OUTPUT */

#include "reducers.h"
    /* For: reduceVariable(), reduceModule(), reduceSubroutine(), etc */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), extract(), skip(), etc */

#include "regexp.h"
    /* For: regexp, regcomp(), regexec() */

/* -- Program Constants -------------------------------------------------------
*/
#define LINESIZE 200

char *DefaultBaseDirectory = CMISS_WWW_ROOT CM_URLPATH;

char *DefaultRulesFile = "cmiss-comment.rules";

char *DefaultModuleListFile  = "../MasterLists/module-list";
char *DefaultCModuleListFile = "../MasterLists/c-module-list";



/* -- CMISS Parser Method -----------------------------------------------------
*/

void parseFile( GROUPFILE *file, RULES *rules, OUTPUT *output )
  {
  GROUP *group;

  group = getGroup( file, rules );
  while( group != NULL )
    {
    addGroupToOutput( group, output );

    discardGroup( file );
    group = getGroup( file, rules );
    }
  }


void parseCMISS( OUTPUT *output, char *listName, char *ruleName )
  {
  RULES        *rules;
  FILESTATE    *file;
  FILESTATE    *list;
  GROUPFILE    *groupFile;
  char         *line;
  char          name[128];

  list = newFileState( listName );
  if( list == NULL )
    {
    fprintf( stderr, "newCmissFileState() returned NULL in parseCMISS()\n" );
    exit( 1 );
    }

  rules = newRules(); 
  if( rules == NULL )
    { 
    fprintf( stderr, "newRules() returned NULL in main()\n" );
    exit( 1 );
    }
  
  if( ! parseRuleFile( rules, ruleName ) )
    {
      fprintf( stderr, "failed to parse rule file %s\n", ruleName );
      exit( 1 );
    }

  line = getLine( list );
  while( line != NULL )
    {
    sscanf( line, "%s\n", name );

    file = newCmissFileState( name );
    if( file == NULL )
      {
      fprintf( stderr, "newCmissFileState() returned NULL in main()\n" );
      exit( 1 );
      }

    groupFile = newGroupFileFromFileState( file );
    if( groupFile == NULL )
      {
      fprintf( stderr, "newGroupFileFromFileState() returned NULL in "
        "main()\n" );
      exit( 1 );
      }

    parseFile( groupFile, rules, output );

    disposeFileState( file );

    discardLine( list );
    line = getLine( list );
    }

  disposeFileState( list );
  disposeRules( rules );
  }

void outputTag( FILE *file, TAG *tag )
  {
  char *s;

  fprintf( file, "%s:\n", tag->tag );
  
  s = tag->body;
  fprintf( file, "  " );
  while( *s != '\0' )
    {
    fputc( *s, file );
    if( *s == '\n' )
      fprintf( file, "  " );

    s++;
    }

  fputc( '\n', file );
  }


void outputGroupListFile( OUTPUT *output, char *groupName, 
  char *baseName )
  {
  FILE      *file;
  GROUPLIST *groupList;
  GROUP     *group;
  GROUPITEM *groupItem;
  TAG       *tag;
  TAGITEM   *tagItem;
  char      *fileName;


  fileName = malloc( strlen( baseName ) + 7 );
  if( fileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in outputGroupListFile()\n" );
    return;
    }

  groupList = lookupGroupList( output, groupName );
  if( groupList == NULL )
    {
    fprintf( stderr, "lookupGroupList() failed in outputGroupListFile()\n" );
    free( fileName );
    return;
    }
  
  sprintf( fileName, "%s.list", baseName );
  file = fopen( fileName, "wt" );
  if( file == NULL )
    {
    fprintf( stderr, "Failed to open file %s in outputGroupListFile(): %s\n",
	     fileName, strerror(errno) );
    free( fileName );
    return;
    }

  groupItem = groupList->firstGroupItem;
  while( groupItem != NULL )
    {
    group = groupItem->group;

    tag = group->headTag;
    outputTag( file, tag );

    tagItem = group->tagList->firstTagItem;
    while( tagItem != NULL )
      {
      tag = tagItem->tag;

      outputTag( file, tag );
      
      tagItem = tagItem->nextTagItem;
      }

    fprintf( file, "\n\n" );

    groupItem = groupItem->nextGroupItem;
    }

  fclose( file );
  free( fileName );
  }
  

void outputDatabaseFiles( char *baseDirectory, OUTPUT *output )
  {
  char *name;

  name = malloc( strlen( baseDirectory ) + 32 ); /* 32 is too much */
  if( name == NULL )
    {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
    }

  sprintf( name, "%s/variables", baseDirectory );
  outputGroupListFile( output, "VARIABLE",   name );

  sprintf( name, "%s/modules", baseDirectory );
  outputGroupListFile( output, "MODULE",     name );

  sprintf( name, "%s/subroutines", baseDirectory );
  outputGroupListFile( output, "SUBROUTINE", name );

  sprintf( name, "%s/functions", baseDirectory );
  outputGroupListFile( output, "FUNCTION",   name );

  sprintf( name, "%s/commands", baseDirectory );
  outputGroupListFile( output, "COMMAND",    name );

  sprintf( name, "%s/comments", baseDirectory );
  outputGroupListFile( output, "COMMENT",    name );

  sprintf( name, "%s/blockdata", baseDirectory );
  outputGroupListFile( output, "BLOCKDATA",    name );

  free( name );
  }


void outputCommandTree( char *baseDirectory, OUTPUT *output )
/*******************************************************************************
LAST MODIFIED : 1 July 1997 by Carey Stevens

DESCRIPTION :
outputCommandTree creates a tree of html files from the command group list
==============================================================================*/
{
  FILE      *file1,*file2,*file3,*file4,*file5;
  GROUPLIST *groupList;
  GROUP     *group;
  GROUPITEM *groupItem;
  TAG       *tag;
  char      *basename;
  char      *CurrentLevelOne,*CurrentLevelTwo,*CurrentLevelThree;
  char      *description_file,*fileName;
  char      *escapedParameters;
  char      *LevelOne,*LevelTwo,*LevelThree;
  char      *Parameters;
  char      *string;
  int       done;
  int       level5_counter;
  
  basename = malloc( strlen( baseDirectory ) + 32 ); 
  if( basename == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  CurrentLevelOne = (char *)malloc(LINESIZE);
  if( CurrentLevelOne == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  CurrentLevelTwo = (char *)malloc(LINESIZE);
  if( CurrentLevelTwo == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  CurrentLevelThree = (char *)malloc(LINESIZE);
  if( CurrentLevelThree == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  description_file = malloc( strlen( baseDirectory ) + 32 ); 
  if( description_file == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  escapedParameters = (char *)malloc(LINESIZE);
  if( escapedParameters == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  fileName = (char *)malloc(LINESIZE);
  if( fileName == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  LevelOne = (char *)malloc(LINESIZE);
  if( LevelOne == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  LevelTwo = (char *)malloc(LINESIZE);
  if( LevelTwo == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  LevelThree = (char *)malloc(LINESIZE);
  if( LevelThree == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  string = (char *)malloc(LINESIZE);
  if( string == NULL )
  {
    fprintf( stderr, "Memory allocation failure in outputDatabaseFiles()\n" );
    exit( 1 );
  }
  
  sprintf( basename, "%s/commands/", baseDirectory );
  
  groupList = lookupGroupList( output, "COMMAND" );
  if( groupList == NULL )
  {
    fprintf( stderr, "lookupGroupList() failed in outputGroupListFile()\n" );
    free( basename );
    free( CurrentLevelOne );
    free( CurrentLevelTwo );
    free( description_file );
    free( escapedParameters );
    free( fileName );
    free( LevelOne );
    free( LevelTwo );
    free( string );
    return;
  }
 
  sprintf( fileName, "%sindex.html", basename );
  file1 = fopen( fileName, "wt" );
  if( file1 == NULL )
  {
    fprintf( stderr, "Failed to open file %s in outputCommandTree(): %s\n",
	     fileName, strerror(errno) );
    free( basename );
    free( CurrentLevelOne );
    free( CurrentLevelTwo );
    free( description_file );
    free( escapedParameters );
    free( fileName );
    free( LevelOne );
    free( LevelTwo );
    free( string );
    return;
  }

  strcpy(LevelOne,"\0"); 
  strcpy(LevelTwo,"\0");
  strcpy(LevelThree,"\0");
  level5_counter=0;



/*   Loop over each command */
  groupItem = groupList->firstGroupItem;
  while( groupItem != NULL )
  {

    group = groupItem->group;    
    tag = group->headTag;
    done=0;
    
/*     Write command information into the final level file */
    level5_counter++;
    sprintf(string,"%sdescrpt_%d.html",basename,level5_counter);
    file5 = fopen( string, "wt" );
    if( file5 == NULL )
    {
      fprintf( stderr, "Failed to open file %s in outputCommandTree(): %s\n",
	       string, strerror(errno) );
      free( basename );
      free( CurrentLevelOne );
      free( CurrentLevelTwo );
      free( description_file );
      free( escapedParameters );
      free( fileName );
      free( LevelOne );
      free( LevelTwo );
      free( string );
      return;
    }
    outputCommand3( file5, group );  
    fclose( file5 );
    sprintf(description_file,"descrpt_%d.html",level5_counter);

    
/*     Parse the command  */
    strcpy(LevelOne,"\0");
    strcpy(LevelTwo,"\0");
    strcpy(LevelThree,"\0");
    sscanf(tag->body,"%[a-z,A-Z,-,_] %[a-z,A-Z,-,_] %[a-z,A-Z,-,_]",LevelOne,LevelTwo,LevelThree);

/*     Determine how many levels the command has, ie where do the parameters start, */
/*     and obtain a string of the parameters */
    if(!strcmp(LevelTwo,"\0"))
    {
      Parameters=tag->body+strlen(LevelOne); 
    }
    else if(!strcmp(LevelThree,"\0"))
    {
      Parameters=tag->body+strlen(LevelOne)+strlen(LevelTwo)+2; 
    }
    else
    {
      Parameters=tag->body+strlen(LevelOne)+strlen(LevelTwo)+strlen(LevelThree)+2; 
    }
    htmlEscapedStrcpy( escapedParameters, Parameters );


/*     Make a new Level 1 entry if it is required */
    if(strcmp(LevelOne,CurrentLevelOne)|| (!strcmp(LevelOne,"\0")))
    {
      sprintf(string,"<A HREF=\"%s.html\">%s</A><BR>",LevelOne,LevelOne);
      fprintf(file1,"%s\n",string);      
      sprintf(string,"%s%s.html",basename,LevelOne);
      if(!strcmp(CurrentLevelOne,"\0"))
      {
        fclose(file2);/* A file2 must already be open */
      }
      if((file2=fopen(string,"w"))==NULL)
      {
	fprintf( stderr, "Failed to open file %s in outputCommandTree(): %s\n",
		 string, strerror(errno) );
        free( basename );
        free( CurrentLevelOne );
        free( CurrentLevelTwo );
        free( description_file );
        free( escapedParameters );
        free( fileName );
        free( LevelOne );
        free( LevelTwo );
        free( string );
        return;
      }
      strcpy(CurrentLevelOne,LevelOne);
      strcpy(CurrentLevelTwo,"\0");
      strcpy(CurrentLevelThree,"\0");
    }
    
/*     Make a new Level 2 entry if it is required */
    if(strcmp(LevelTwo,CurrentLevelTwo) || (!strcmp(LevelTwo,"\0")))
    {
      if(!strcmp(LevelTwo,"\0"))
      {
        sprintf(string,"<A HREF=\"%s\">%s</A><BR>",description_file,CurrentLevelOne); 
        done=1;
      }
      else
      {
        sprintf(string,"%s <A HREF=\"%s%s.html\">%s</A><BR>",LevelOne,LevelOne,LevelTwo,
          LevelTwo);
      }      
      fprintf(file2,"%s\n",string);
      sprintf(string,"%s%s%s.html",basename,LevelOne,LevelTwo);      
      if(strcmp(CurrentLevelTwo,"\0"))
      {
        fclose(file3);
      }
      if((file3=fopen(string,"w"))==NULL)
      {
	fprintf( stderr, "Failed to open file %s in outputCommandTree(): %s\n",
		 string, strerror(errno) );
        free( basename );
        free( CurrentLevelOne );
        free( CurrentLevelTwo );
        free( description_file );
        free( escapedParameters );
        free( fileName );
        free( LevelOne );
        free( LevelTwo );
        free( string );
        return;
      }
      strcpy(CurrentLevelTwo,LevelTwo);
    }    
    
/*     Make a new Level 3 entry if it is required */
    if((strcmp(LevelThree,CurrentLevelThree) || (!strcmp(LevelThree,"\0"))) && (done==0)) 
    {
      if(!strcmp(LevelThree,"\0"))
      {
        sprintf(string,"%s <A HREF=\"%s\">%s</A><BR>",CurrentLevelOne,description_file,
          CurrentLevelTwo); 
      }
      else
      {
        sprintf(string,"%s %s <A HREF=\"%s%s%s.html\">%s</A><BR>",LevelOne,LevelTwo,LevelOne,LevelTwo,
          LevelThree,LevelThree);
      }      
      
      fprintf(file3,"%s\n",string);
      sprintf(string,"%s%s%s%s.html",basename,LevelOne,LevelTwo,LevelThree);

      if(strcmp(CurrentLevelThree,"\0"))
      {
        fclose(file4);
      }
      if((file4=fopen(string,"w"))==NULL)
      {
	fprintf( stderr, "Failed to open file %s in outputCommandTree(): %s\n",
		 string, strerror(errno) );
        free( basename );
        free( CurrentLevelOne );
        free( CurrentLevelTwo );
        free( description_file );
        free( escapedParameters );
        free( fileName );
        free( LevelOne );
        free( LevelTwo );
        free( string );
        return;
      }
      strcpy(CurrentLevelThree,LevelThree);
    }

    
    if(strcmp(LevelThree,"\0"))
    {
      sprintf(string,"%s %s <A HREF=\"%s\">%s</A>%s<BR>",CurrentLevelOne,CurrentLevelTwo,
        description_file,CurrentLevelThree,escapedParameters);    
      fprintf(file4,"%s",string);
      outputCommand2( file4, group );
    }



      groupItem = groupItem->nextGroupItem;
  }
  fclose( file1 );
  fclose( file2 );
  fclose( file3 );
  fclose( file4 );
  free( basename );
  free( CurrentLevelOne );
  free( CurrentLevelTwo );
  free( description_file );
  free( escapedParameters );
  free( fileName );
  free( LevelOne );
  free( LevelTwo );
  free( string );
}

void printHelp( char *programName )
  {
  printf( "Syntax: %s [-b <base directory>]\n", programName );
  }

void sortCommandGroup( OUTPUT *output )
/*******************************************************************************
LAST MODIFIED : 1 July 1997 by Carey Stevens

DESCRIPTION :
sortCommandGroup sorts the command group list into alphabetical order
==============================================================================*/
{
  GROUPLIST *groupList;
  GROUP     *group1,*group2;
  GROUPITEM *groupItem,*groupItem_tmp1,*groupItem_tmp2,*groupItem_tmp3;
  TAG       *tag1,*tag2;
  char      *word1,*word2,*word1_tmp,*word2_tmp;
  int       changed;
    
  groupList = lookupGroupList( output, "COMMAND" );
  if( groupList == NULL )
  {
    fprintf( stderr, "lookupGroupList() failed in outputGroupListFile()\n" );
    return;
  }
  
  word1 = (char *)malloc(LINESIZE);
  if( word1 == NULL )
  {
    fprintf( stderr, "Memory allocation failure in sortCommandGroup()\n" );
    exit( 1 );
  }
  word2 = (char *)malloc(LINESIZE);
  if( word2 == NULL )
  {
    fprintf( stderr, "Memory allocation failure in sortCommandGroup()\n" );
    exit( 1 );
  }
  word1_tmp = (char *)malloc(LINESIZE);
  if( word1_tmp == NULL )
  {
    fprintf( stderr, "Memory allocation failure in sortCommandGroup()\n" );
    exit( 1 );
  }
  word2_tmp = (char *)malloc(LINESIZE);
  if( word2_tmp == NULL )
  {
    fprintf( stderr, "Memory allocation failure in sortCommandGroup()\n" );
    exit( 1 );
  }
  
  changed=1;
  while(changed == 1)
  {
    changed=0;
    groupItem = groupList->firstGroupItem;
    while( groupItem->nextGroupItem->nextGroupItem != NULL )
    {
      if(groupItem == groupList->firstGroupItem)
      {
        groupItem_tmp1=groupItem;
        groupItem_tmp2=groupItem->nextGroupItem;
        groupItem_tmp3=groupItem->nextGroupItem->nextGroupItem;
        group1 = groupItem->group;
        group2 = groupItem->nextGroupItem->group;
        tag1 = group1->headTag;      
        tag2 = group2->headTag;
        if(strcmp(tag1->body,tag2->body) > 0)
        {
          groupItem->nextGroupItem->nextGroupItem = groupItem_tmp1;
          groupItem->nextGroupItem = groupItem_tmp3;
          groupList->firstGroupItem = groupItem_tmp2;
          changed=1;
        }
      }
      
      groupItem_tmp1=groupItem->nextGroupItem;
      groupItem_tmp2=groupItem->nextGroupItem->nextGroupItem;
      groupItem_tmp3=groupItem->nextGroupItem->nextGroupItem->nextGroupItem;
      group1 = groupItem->nextGroupItem->group;
      group2 = groupItem->nextGroupItem->nextGroupItem->group;
      tag1 = group1->headTag;      
      tag2 = group2->headTag;
      if(strcmp(tag1->body,tag2->body) > 0)
      {
        groupItem->nextGroupItem->nextGroupItem->nextGroupItem = groupItem_tmp1;
        groupItem->nextGroupItem->nextGroupItem->nextGroupItem->nextGroupItem = groupItem_tmp3;
        groupItem->nextGroupItem= groupItem_tmp2;
        changed=1;        
      }
      groupItem = groupItem->nextGroupItem;
    }
  }
}


/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  int     arg;
  char   *optPtr;

  char   *baseDirectory   = DefaultBaseDirectory;
  char   *rulesFile       = DefaultRulesFile;
  char   *moduleListFile  = DefaultModuleListFile;
  char   *cModuleListFile = DefaultCModuleListFile;

  OUTPUT *output;

  if( argc > 1 )
    {
    arg = 1;

    while( arg < argc )
      {
      if( argv[ arg ] != NULL && argv[ arg ][0] == '-' )
        {
        /* an option */
        optPtr = &(argv[ arg ][1]);

        while( *optPtr != NULL )
          {
          switch( *optPtr )
            {
            case 'b': 
              arg++;
              if( argv[ arg ] != NULL )
                {
                baseDirectory = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -b option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'r':
              arg++;
              if( argv[ arg ] != NULL )
                {
                rulesFile = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -r option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'l':
              arg++;
              if( argv[ arg ] != NULL )
                {
                moduleListFile = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -l option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'c':
              arg++;
              if( argv[ arg ] != NULL )
                {
                cModuleListFile = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -c option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'h':
              printHelp( argv[0] );
              return 0;

            default: 
              fprintf( stderr, "Illegal Option: %c\n", *optPtr );
              printHelp( argv[0] );
              exit( 1 );
            }

          optPtr++;
          }
        }
      else
        {
        if( argv[ arg ] != NULL )
          {
          fprintf( stderr, "Unknown Parameter: %s\n", argv[ arg ] );
          printHelp( argv[0] );
          exit( 1 );
          }
        }

      arg++;
      }
    }

  output = newOutput();
  if( output == NULL )
    {
    fprintf( stderr, "newOutput() return NULL in main()\n" );
    exit( 1 );
    }

  parseCMISS( output, moduleListFile, rulesFile );
  parseCMISS( output, cModuleListFile, rulesFile );
  outputDatabaseFiles( CMISS_WWW_ROOT LOOKUP_URLPATH, output );
  sortCommandGroup( output );
  outputCommandTree( baseDirectory,output );
  
  disposeOutput( output );

  return 0;
  }

