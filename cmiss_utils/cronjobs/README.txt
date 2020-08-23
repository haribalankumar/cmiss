This directory contains files used to run the daily builds of cm, unemap and cmgui plus the example testing.


In February 2009 a lot of the jobs were updated. The master cronjob now runs 
on bioeng85 not esu8. Stuff that ran on bioeng22 now runs on bioeng1031.


There is now only one crontab file used to control the overnight build and 
testing. A backup of this crontab file is listed here as crontab_bioeng85.file
and is under CVS control. All scripts are also under CVS control.

crontab_bioeng85.file contains the commands used by cron to call the scripts 
that perform the required overnight tasks including:
- an svn update
- starting the main build jobs
- updating the cm documentation
- the final job (which updates the website and emails results)

If you modify the crontab file on bioeng85 you should update this backup file 
too and commit your changes. This also applies to any of the scripts mentioned
below.

The key script file is cmiss_cronjobs_daily_master.sh which starts builds on 
several machines and also sets scripts running on other machines.


A brief summary of what the daily scripts do is as follows:

cmiss_cronjobs_daily_master.sh  - Updates master log files
                                Starts bioeng85 cron job (see below)
                                Starts bioeng1031 cron job (see below)
                                Starts esu8 cron job (see below)
                                Checks cmgui, unemap are up to date

cmiss_cronjobs_daily_bioeng85.sh 
                                - Updates bioeng85 log files
                                Runs x86 32 bit build jobs
                                Makes the following: perl interpreter,
                                linear solvers, cm, cmgui and zinc  
                                Tests cm/cmgui x86 32 bit examples 

cmiss_cronjobs_daily_bioeng1031.sh    
                                - Updates bioeng1031 log files
                                Runs x86 64 bit build jobs
                                Makes the following: perl interpreter,
                                linear solvers, cm, cmgui and zinc  
                                Starts testing on bioeng22

cmiss_cronjobs_daily_esu8.sh    - Updates esu8 log files
                                Runs 64 bit IRIX build jobs
                                Makes the following: linear solves and
                                perl interpreter
                                Tests cmgui examples

cmiss_cronjobs_daily_bioeng22.sh  
                                - Updates bioeng22 log files
                                Tests cm/cmgui x86 64 bit examples

cmiss_cronjobs_daily_final.sh   - Updates final log files
                                Updates examples website
                                Emails results of daily build


cmiss_cronjobs_daily_cm-docs.sh - Create CM documentation

For more information see http://www.cmiss.org/cmgui/wiki/CmguiNightlyBuild/view

