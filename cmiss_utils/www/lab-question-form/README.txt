CMISS Tutorials
http://www.cmiss.org/cm/wiki/CMISSTutorials/
(Revamped by Martyn Nash Feb 2006)


Descriptive Tutorials:
----------------------

These CMISS (cm) tutorials were created some time ago and have gone through
several iterations of refinement.  The latest has involved transferring
them to the CMISS content management system in Feb 2006. Latest files at:
http://www.cmiss.org/cm/wiki/CMISSTutorials/


On-line Tutorial Quiz's:
------------------------

There are a set of multi-choice questions associated with each
tutorial (except #1).  Source code for the question form CGI script is
stored in: 
/hpc/cmiss/cmiss_utils/www/lab-question-form/

The source files (lab-question-form.c and Makefile) are controlled
using CVS, so you need to check it out to your own directory to make
changes (eg. "cvs co cmiss_utils/www/lab-question-form"). Once edits
have beem committed, the CGI script needs to be 're-made':
> log in to bioeng1031 as cmiss
> cd cmiss_utils/www/lab-question-form
> make

The Makefile compiles lab-question-form.cgi and copies it to:
/hpc/cmissweb/cgi-bin/
from where it is accessible to the live web server:
http://cmiss.bioeng.auckland.ac.nz/

The CGI script reads/writes data from/to:
/hpc/cmissweb/lab_questions/

Question data is read from the 'questions' subdirectory (surprise!),
which contains subdirectories for each tutorial. If additions or typos 
are needed, these files can be edited directly.

Even more surprisingly, form data filled out on the web (ie. answers)
are stored under the 'answers' subdirectory, which contains a
subdirectory named using the UPI entry in the web form.  Within each
UPI directory, there is a subdirectory for each tutorial.  Corrected answer
sheets ('answers.txt') are written in each subdirectory, along with a
total mark.  Neither answers nor marks are written to the browser
screen.

If a student has submitted one set of answers, then attempts to
resubmit another set for the same tutorial, then the new answers will be
appended to the previous answers.txt file for that tutorial.
