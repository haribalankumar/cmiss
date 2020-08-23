
#ifndef WEBCMISS_H
#define WEBCMISS_H

/* -- Module Datatypes --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif


/* -- Module Constants --------------------------------------------------------
*/
#undef TESTING
#undef SETGID

#define WEBMASTER "webmaster@www.esc.auckland.ac.nz"
#define HTTP_URL  "http://www.esc.auckland.ac.nz/"
#define LOGIN_URL "http://www.bioeng.auckland.ac.nz/cmiss/hpc_access.php"

#ifdef TESTING

#define CMISS_WWW_URL           HTTP_URL "People/Staff/Norris/"
#define EXAMPLE_FILES_URL       HTTP_URL "Groups/Bioengineering/CMISS/help/examples/"

#define ESU_ROOT                "/www/People/Staff/Norris/"
#define ESU_WEB                 ESU_ROOT "webcmiss/"
#define ESU_NOWEB               "/usr/people/norris/webcmiss/"
            
#else

#define CMISS_WWW_URL           HTTP_URL "Groups/Bioengineering/CMISS/"
#define EXAMPLE_FILES_URL       CMISS_WWW_URL "help/examples/"

#define ESU_ROOT                "/product/cmiss/"
#define ESU_WEB                 ESU_ROOT "www/webcmiss/"
#define ESU_NOWEB               "/product/cmiss/webcmiss/"

#endif

#define WEBCMISS_URL            CMISS_WWW_URL "webcmiss/"

#define ESU_ETC                 ESU_NOWEB "etc/"
#define ESU_IMAGES              ESU_NOWEB "images/"
#define ESU_HTML                ESU_NOWEB "html/"
#define ESU_EXAMPLES_PATH       ESU_ROOT "examples/"
#define ESU_RULE_FILE_NAME      ESU_ETC "job.rules"
#define ESU_MASTER_USER_LIST    ESU_ETC "password"
#define ESU_MASTER_HOST_LIST    ESU_ETC "host"


/* Script/Server names */

#define ACCOUNT_DETAILS         "account-details.cgi"
#define CHANGE_EXAMPLE          "change-example.cgi"
#define CHANGE_PASSWORD         "change-password.cgi"
#define CMGUI_PROGRAM           "cmgui-deliver.cgi"
#define DELETE_DIR              "delete-dir.cgi"
#define DELETE_FILE             "delete-file.cgi"
#define EDIT_FILE               "edit-file.cgi"
#define JOB_CONTROL             "job-control.cgi"
#define LIST_PROBLEM            "list-problem.cgi"
#define REVERT_FILE             "revert-file.cgi"
#define REVERT_FILECONF         "revert-confirmed.cgi"
#define SAVE_FILE               "save-file.cgi"
#define VIEW_FILE               "view-file.cgi"

#define SUBMIT_SERVER           "submission-server"
#define ABORT_SERVER            "abort-server"

#define BORDER_HEAD_FILE        ESU_HTML "border_head.html"
#define BORDER_HEADER_FILE      ESU_HTML "border_header.html"
#define BORDER_FOOTER_FILE      ESU_HTML "border_footer.html"

#define JOB_FILE_NAME           "job.info"
#define IMAGE_DIR               "../webcmiss/images/"

/* Scripts for the simplified parameter forms */

#define PARAM_MENU              "parameter-menu.cgi"
#define PARAM_EXAMPLE           "parameter-example.cgi"
#define PARAM_SET_EXAMPLE       "parameter-set-example.cgi"
#define PARAM_JOB               "parameter-job.cgi"
#define PARAM_VIEW              "parameter-view.cgi"

#define PARAM_DIR                ESU_WEB "param_examples/"
#define PARAM_MENU_FILE          PARAM_DIR "parameter_menu.html"

#endif
