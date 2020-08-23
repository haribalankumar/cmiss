#
#  Makefile for common-browser
#

SHELL =/bin/sh

CGI_DIR=../CGIDecode/
UTE_DIR=../Utilities/
 
# Normal installation 
INSTALL_DIR=/www/Groups/Bioengineering/CMISS/scripts/
BIN_DIR=/product/cmiss/webcmiss/bin/
HPC_BIN=webcmiss@hpc2:bin/

# When Testing
#INSTALL_DIR=/www/People/Staff/Norris/cgi-bin/
#BIN_DIR=$(HOME)/webcmiss/bin/
#HPC_BIN=norris@hpc2:webcmiss/bin/

# -- Processor Directives -----------------------------------------------------
#     

CC       =cc
CFLAGS   =-n32 -O0 -pedantic -fullwarn -woff 841
CINCLUDES=-I$(CGI_DIR) -I$(UTE_DIR)

LIBFLAGS =-L. -lform -L$(CGI_DIR) -lcgi -L$(UTE_DIR) -lUte

AR       =ar
AFLAGS   =rvf


.c.o: 
	$(CC) -c $(CFLAGS) $(CINCLUDES) -o $@ $<

.c: 
	$(CC) -o $@ $< $(CFLAGS) $(CINCLUDES)

LIB = libform.a
OBJ = directory.o print-error.o user-managment.o get-user.o job.o transfer.o\
      password.o host.o border.o resource.o

CGI = list-problem.cgi job-control.cgi account-details.cgi\
      delete-dir.cgi change-password.cgi check-password.cgi\
      change-example.cgi revert-file.cgi\
      edit-file.cgi save-file.cgi view-file.cgi delete-file.cgi cmgui-deliver.cgi\
      parameter-menu.cgi parameter-example.cgi parameter-set-example.cgi\
      parameter-job.cgi

SERVER = submission-server abort-server
PASSWD = webcmiss_passwd
TESTER = server-test


all: lib cgip utils cgi server passwd

lib: $(LIB)
cgi: $(CGI)
server: $(SERVER)
passwd: $(PASSWD)

cgip:
	cd $(CGI_DIR) ; make

utils:
	cd $(UTE_DIR) ; make


# -- Specific Targets ---------------------------------------------------------
#

$(LIB): $(OBJ)
	$(AR) $(AFLAGS) $@ $(OBJ)

server-test: server-test.o $(LIB)
	$(CC) -o $@ server-test.o $(CFLAGS) $(LIBFLAGS)

submission-server: submission-server.o $(LIB)
	$(CC) -o $@ submission-server.o $(CFLAGS) $(LIBFLAGS)

abort-server: abort-server.o $(LIB)
	$(CC) -o $@ abort-server.o $(CFLAGS) $(LIBFLAGS)

webcmiss_passwd: webcmiss_passwd.o $(LIB)
	$(CC) -o $@ webcmiss_passwd.o $(CFLAGS) $(LIBFLAGS)



list-problem.cgi: list-problem.o $(LIB)
	$(CC) -o $@ list-problem.o $(CFLAGS) $(LIBFLAGS)

delete-dir.cgi: delete-dir.o $(LIB)
	$(CC) -o $@  delete-dir.o $(CFLAGS) $(LIBFLAGS)

change-example.cgi: change-example.o $(LIB)
	$(CC) -o $@  change-example.o $(CFLAGS) $(LIBFLAGS)

account-details.cgi: account-details.o $(LIB)
	$(CC) -o $@  account-details.o $(CFLAGS) $(LIBFLAGS)

job-control.cgi: job-control.o $(LIB)
	$(CC) -o $@  job-control.o $(CFLAGS) $(LIBFLAGS)

change-password.cgi: change-password.o $(LIB)
	$(CC) -o $@  change-password.o $(CFLAGS) $(LIBFLAGS)

check-password.cgi: check-password.o $(LIB)
	$(CC) -o $@  check-password.o $(CFLAGS) $(LIBFLAGS)

edit-file.cgi: edit-file.o $(LIB)
	$(CC) -o $@ edit-file.o $(CFLAGS) $(LIBFLAGS)

save-file.cgi: save-file.o $(LIB)
	$(CC) -o $@ save-file.o $(CFLAGS) $(LIBFLAGS)

view-file.cgi: view-file.o $(LIB)
	$(CC) -o $@  view-file.o $(CFLAGS) $(LIBFLAGS)

delete-file.cgi: delete-file.o $(LIB)
	$(CC) -o $@  delete-file.o $(CFLAGS) $(LIBFLAGS)

revert-file.cgi: revert-file.o $(LIB)
	$(CC) -o $@ revert-file.o $(CFLAGS) $(LIBFLAGS)

cmgui-deliver.cgi: cmgui-deliver.o $(LIB)
	$(CC) -o $@  cmgui-deliver.o $(CFLAGS) $(LIBFLAGS)

parameter-menu.cgi: parameter-menu.o $(LIB)
	$(CC) -o $@  parameter-menu.o $(CFLAGS) $(LIBFLAGS)

parameter-example.cgi: parameter-example.o $(LIB)
	$(CC) -o $@  parameter-example.o $(CFLAGS) $(LIBFLAGS)

parameter-set-example.cgi: parameter-set-example.o $(LIB)
	$(CC) -o $@  parameter-set-example.o $(CFLAGS) $(LIBFLAGS)

parameter-job.cgi: parameter-job.o $(LIB)
	$(CC) -o $@  parameter-job.o $(CFLAGS) $(LIBFLAGS)


install: $(INSTALL_DIR) $(CGI) $(SERVER)
	cp $(CGI) $(INSTALL_DIR)
	chmod uog+rx $(INSTALL_DIR)/*.cgi
	chmod ug+s $(INSTALL_DIR)/change-password.cgi
	cp $(PASSWD) hpc_cmgui.sh $(BIN_DIR)
	chmod u+rx $(BIN_DIR)/hpc_cmgui.sh
# Only when scripts are setgid to cmiss
#	chgrp cmiss $(INSTALL_DIR)/change-example.cgi $(INSTALL_DIR)/revert-confirmed.cgi
#	chmod g+s $(INSTALL_DIR)/change-example.cgi $(INSTALL_DIR)/revert-confirmed.cgi

rinstall: $(SERVER)
	rcp $(SERVER) mail-user.sh $(HPC_BIN)

clean:
	rm -f $(SERVER) $(PASSWD) *.cgi *.o *.a *~

