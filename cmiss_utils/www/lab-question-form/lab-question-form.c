/* -- Include Files -- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <cgi-decode.h> /* For: CGI handling functions */

/* -- Defines -- */

#define NUM_LABS 6
#define MAX_NUM_QUESTIONS 20
#define MAX_NUM_OPTIONS 6
#define STR_LEN 1000

#define INPUT_TYPE 1
#define TEXT_AREA_TYPE 2

#define CGI_BOUNDARY "cmisscgiboundary"

#define PRINT_TO_SCREEN 1

/* -- Structs -- */

typedef struct {
  int number_of_questions;
  int number_of_options[MAX_NUM_QUESTIONS];
  char question_string[MAX_NUM_QUESTIONS][STR_LEN];
  char option_strings[MAX_NUM_QUESTIONS][MAX_NUM_OPTIONS][STR_LEN];
  int answer[MAX_NUM_QUESTIONS];
  int correct_answer[MAX_NUM_QUESTIONS];
} QUESTIONAIRE;

/* -- Global Variables -- */

/* -- Main Program -- */

void write_cgi_start()
{
  printf("Content-type: multipart/mixed; boundary=\""CGI_BOUNDARY"\"\n\n");
  fflush(stdout);
}

void write_internal_boundary(char *contenttype)
{
  printf("\n--"CGI_BOUNDARY"\n");
  printf("Content-type: %s\n\n",contenttype);
  fflush(stdout);  
}

void write_internal_boundary_name(char *contenttype,
  char *name)
{
  printf("\n--"CGI_BOUNDARY"\n");
  printf("Content-type: %s\n",contenttype);
  printf("Content-Disposition: filename=\"%s\"\n\n",name);
  fflush(stdout);  
}

void write_cgi_end()
{
  printf("\n--"CGI_BOUNDARY"--\n");
  fflush(stdout);  
}

void write_html_header(char *titlestring,
  char *headingstring)
{
  write_internal_boundary("text/html");
  printf("<html><head>\n");
  printf("<title>%s</title>\n",titlestring);
  printf("</head><body>\n");
  printf("<h1>%s</h1>\n",headingstring);
  printf("<p></p><hr>\n");
}

void write_html_footer()
{
  printf("</body><hr>\n");
  printf("<a href=\""CMISS_WWW_URL"\">Return to CMISS Tutorials page</a>\n");
  printf("</html>\n");
}

void write_error(char *error_str)
{
  write_html_header("CMISS error","Error:");
  printf("%s\n",error_str);
  fflush(stdout);
/*   write_html_footer(); */
}

int get_laboratory(CGI *cgi,
  int *laboratory)
{
  char *laboratorystring;
  int error=1;

  if(laboratorystring = lookupString(cgi,"laboratory"))
  {
    *laboratory=atoi(laboratorystring);
    if(*laboratory < 1 || *laboratory > NUM_LABS)
    {
      write_error("Invalid laboratory");
      error=0;
    }
  }
  else
  {
    write_error("No laboratory selected");
    error=0;
  }

  return error;
}

void write_laboratory_form()
{
  int i;
  
  printf("<p>Welcome to the CMISS Tutorials Quiz page. Please select the Quiz you wish to take.</p><hr>\n");
  printf("<p><form method=\"GET\" action=\""CMISS_CGI_URL
    "cgi-bin/lab-question-form.cgi\">\n");
  printf("<input type=\"hidden\" name=\"stage\" value=\"1\">\n");
  for(i=1;i<=NUM_LABS;i++)
  {
    printf("<input type=\"radio\" name=\"laboratory\" value=\"%d\"> "
      "Tutorial %d Quiz<br>\n",i,i);
  }
  printf("<br>\n");
  printf("<input type=\"submit\" value=\"Select Quiz\">\n");
  printf("<input type=\"reset\" value=\"Reset\">\n");
  printf("</form></p>\n");
}

void write_question_form(int laboratory, QUESTIONAIRE *questionaire)
{
  int i,j;
  FILE *questionfile;
  char questionfilename[STR_LEN],line[STR_LEN];

  sprintf(questionfilename,"%s/questions/lab%d/questions.txt",QUESTION_PATH,laboratory);
  if(questionfile=fopen(questionfilename,"r"))
  {
    fscanf(questionfile,"Number of questions: %d\n",&(questionaire->number_of_questions));
    for(i=1;i<=questionaire->number_of_questions;i++)
    {
      fscanf(questionfile,"Question string: ");
      fgets(questionaire->question_string[i],STR_LEN,questionfile);
      fscanf(questionfile,"Number of options: %d\n",&(questionaire->number_of_options[i]));
      for(j=1;j<=questionaire->number_of_options[i];j++)
      {
        sprintf(line,"Option %d: ",j);
        fscanf(questionfile,line,line);
        fgets(questionaire->option_strings[i][j],STR_LEN,questionfile);
      }
      fscanf(questionfile,"Correct answer: %d\n",&(questionaire->correct_answer[i]));
    }
    fclose(questionfile);
    printf("<p><form method=\"GET\" action=\""CMISS_CGI_URL
      "cgi-bin/lab-question-form.cgi\">\n");
    printf("<input type=\"hidden\" name=\"stage\" value=\"2\">\n");
    printf("<input type=\"hidden\" name=\"laboratory\" value=\"%d\">\n",laboratory);
    printf("UPI: <input name=\"UPI\" size=\"7\" maxlength=\"7\"> <br>\n");
    printf("<ol>\n");
    for(i=1;i<=questionaire->number_of_questions;i++)
    {
      printf("<li> %s\n",questionaire->question_string[i]);
      printf("<ol>\n");
      for(j=1;j<=questionaire->number_of_options[i];j++)
      {
        printf("<li><input type=\"radio\" name=\"question%d\" value=\"%d\"> %s\n",i,j,questionaire->option_strings[i][j]);
      }
      printf("</ol>\n");
    }
    printf("</ol>\n");
//    printf("<p>By clicking on the \"Submit answers\" button I acknowledge that I understand the School of Engineering's policy on cheating and that all information supplied is correct and these answers are solely my work.</p>\n");
    printf("<input type=\"submit\" value=\"Submit answers\">\n");
    printf("<input type=\"reset\" value=\"Reset\">\n");
    printf("</form></p>\n");
  }
  else
  {
    write_error("Could not open questions file");
  }
}

void process_question_form(CGI *cgi, int laboratory, QUESTIONAIRE *questionaire)
{
  int i,j,num_correct,result;
  FILE *file;
  char filename[STR_LEN],line[STR_LEN],*string,*UPI;
  time_t timeval;
  mode_t mode;
  struct tm *time_structure;
  static char *month[12]={
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
  };

  UPI = lookupString(cgi,"UPI");
//  if(! strcmp(UPI,""))
//  {
//    write_error("No UPI specified");
//  }
  if( strlen(UPI)<6 || strlen(UPI)>7 ||
      strncasecmp(&UPI[0],"a",1)<0 || strncasecmp(&UPI[0],"z",1)>0 ||
      strncasecmp(&UPI[1],"a",1)<0 || strncasecmp(&UPI[1],"z",1)>0 ||
      strncasecmp(&UPI[2],"a",1)<0 || strncasecmp(&UPI[2],"z",1)>0 ||
      strncasecmp(&UPI[strlen(UPI)-3],"0",1)<0 || strncasecmp(&UPI[strlen(UPI)-3],"9",1)>0 ||
      strncasecmp(&UPI[strlen(UPI)-2],"0",1)<0 || strncasecmp(&UPI[strlen(UPI)-2],"9",1)>0 ||
      strncasecmp(&UPI[strlen(UPI)-1],"0",1)<0 || strncasecmp(&UPI[strlen(UPI)-1],"9",1)>0
      )
  {
    write_error("You must specify your correct UPI.<BR>Use the back button to retrieve your form answers.");
  }
  else
  {
    sprintf(filename,"%s/questions/lab%d/questions.txt",QUESTION_PATH,laboratory);
    if(file=fopen(filename,"r"))
    {
      fscanf(file,"Number of questions: %d\n",&(questionaire->number_of_questions));
      for(i=1;i<=questionaire->number_of_questions;i++)
      {
        fscanf(file,"Question string: ");
        fgets(questionaire->question_string[i],STR_LEN,file);
        fscanf(file,"Number of options: %d\n",&(questionaire->number_of_options[i]));
        for(j=1;j<=questionaire->number_of_options[i];j++)
        {
          sprintf(line,"Option %d: ",j);
          fscanf(file,line,line);
          fgets(questionaire->option_strings[i][j],STR_LEN,file);
        }
        fscanf(file,"Correct answer: %d\n",&(questionaire->correct_answer[i]));
        sprintf(line,"question%d",i);
        if(string = lookupString(cgi,line))
        {
          questionaire->answer[i]=atoi(string);
          free(string);
        }
        else
        {
          questionaire->answer[i]=0;
        }
      }
      fclose(file);
      mode=504;
      sprintf(filename,"%s/answers/%s",QUESTION_PATH,UPI);
      result=mkdir(filename,mode);
      if(result == 0 || errno == EEXIST)
      {
        sprintf(filename,"%s/answers/%s/lab%d",QUESTION_PATH,UPI,laboratory);
        result=mkdir(filename,mode);
        if(result == 0 || errno == EEXIST)
        {
          sprintf(filename,"%s/answers/%s/lab%d/answers.txt",QUESTION_PATH,UPI,laboratory);
          if(file=fopen(filename,"a+"))
          {
            fprintf(file,"CMISS Tutorial Quiz %d Solutions\n",laboratory);
            timeval=time(NULL);
            time_structure=localtime(&timeval);
            fprintf(file,"%02d%3s%02d,%02d:%02d:%02d\n\n", 
              time_structure->tm_mday, month[time_structure->tm_mon], 
              time_structure->tm_year+1900, time_structure->tm_hour+13, 
              time_structure->tm_min, time_structure->tm_sec);
            printf("<B>%02d%3s%02d,%02d:%02d:%02d</B><P>\n\n", 
              time_structure->tm_mday, month[time_structure->tm_mon], 
              time_structure->tm_year+1900, time_structure->tm_hour+13, 
              time_structure->tm_min, time_structure->tm_sec);
            fprintf(file,"UPI: %s\n",UPI);
            printf("<B>UPI:</B> %s<P><PRE>\n",UPI);
            num_correct=0;
            for(i=1;i<=questionaire->number_of_questions;i++)
            {
              fprintf(file,"%d. %s",i,questionaire->question_string[i]);
              printf("%d. %s",i,questionaire->question_string[i]);
              for(j=1;j<=questionaire->number_of_options[i];j++)
              {
                fprintf(file,"   %d. %s",j,questionaire->option_strings[i][j]);
                printf("   %d. %s",j,questionaire->option_strings[i][j]);
              }
              fprintf(file,"   ANSWER: %d\n",questionaire->answer[i]);
              printf("   YOUR ANSWER: %d\n",questionaire->answer[i]);
              if(questionaire->answer[i] == questionaire->correct_answer[i])
              {
                fprintf(file,"   CORRECT\n");
#ifdef	PRINT_TO_SCREEN		
//print answers to screen
		printf("   CORRECT\n");
#endif
                num_correct++;
              }
              else
              {
                fprintf(file,"   INCORRECT\n");
#ifdef  PRINT_TO_SCREEN
//print answers to screen
                printf("   INCORRECT\n");
#endif
                fprintf(file,"   Correct answer: %d\n",questionaire->correct_answer[i]);
#ifdef  PRINT_TO_SCREEN
//print answers to screen
                printf("   Correct answer: %d\n",questionaire->correct_answer[i]);
#endif
              }
              fprintf(file,"   --------------\n");
              printf("   --------------\n");
            }
            fprintf(file,"\nNumber of correct answers: %d/%d\n------------------------------\n",num_correct,questionaire->number_of_questions);
#ifdef  PRINT_TO_SCREEN
//print answers to screen
            printf("\nNumber of correct answers: %d/%d\n",num_correct,questionaire->number_of_questions);
#endif
            printf("</PRE>\n");
            fclose(file);
          }
          else
          {
            write_error("Could not create answers laboratory directory");
          }
        }
        else
        {
          write_error("Could not create answers UPI directory");
        }
      }
      else
      {
        write_error("Could not open questions file");
      }
    }
    else
    {
      write_error("No email address specified");
    }
    free(UPI);
  }
}

/* -- The Program Entry Point -- */
  
int main(int argc,
  char *argv[])
{
  CGI *cgi;
  char *stagestring,headerstring[STR_LEN];
  int laboratory,stage;
  QUESTIONAIRE questionaire;
  
  write_cgi_start();
  
  cgi = getCGIEnvironment(argc,argv);
  if(cgi == NULL)
  {
    write_error("getCGIEnvironment() failed in main()");
  }
  else
  {
    if(!(stagestring = lookupString(cgi,"stage")))
    {
      stage=0;
    }
    else
    {
      stage=atoi(stagestring);
    }    
    switch (stage)
    {
      case 0:
      {
        write_html_header("CMISS Tutorial Quiz Selection","CMISS Tutorial Quiz Selection");
        write_laboratory_form();
        write_html_footer();
        break;
      }
      case 1:      
      {
        if(get_laboratory(cgi,&laboratory))
        {
          sprintf(headerstring,"CMISS Tutorial %d Quiz",laboratory);
          write_html_header(headerstring,headerstring);
          write_question_form(laboratory,&questionaire);
          write_html_footer();
        }
        break;
      }
      case 2:
      {
        if(get_laboratory(cgi,&laboratory))
        {
          sprintf(headerstring,"CMISS Tutorial %d Quiz Answers",laboratory);
          write_html_header(headerstring,headerstring);
          process_question_form(cgi,laboratory,&questionaire);
          write_html_footer();
        }
        break;
      }
      default:
      {
        write_error("Invalid stage number");
        break;
      }
    }
  } 

  write_cgi_end();
  disposeCGI(cgi);
  fflush(stdout);

  return 0;
}
