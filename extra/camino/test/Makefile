# Make file for the camino tests.
# 
# $Id: Makefile,v 1.5 2006/06/02 16:31:27 ucacpco Exp $
UNAME=$(shell uname)
CYGWIN=$(findstring CYGWIN,$(UNAME))
SCRIPTDIR=$(shell pwd)

ifeq ($(CYGWIN),)
    export CLASSPATH=$(SCRIPTDIR)/..:$(SCRIPTDIR)/test:$(SCRIPTDIR)/junit.jar
else 
    CAMINO=$(shell cygpath -w $(SCRIPTDIR)/..)
    TEST=$(shell cygpath -w $(SCRIPTDIR)/test)
    JUNIT=$(shell cygpath -w $(SCRIPTDIR)/junit.jar)
    export CLASSPATH="$(CAMINO);$(TEST);$(JUNIT)"    	
endif


all:				AllTests tractography/Images


AllTests:			AllTests.class

AllTests.class::		
				javac -classpath $(CLASSPATH) AllTests.java

tractography/Images:		tractography/Images.class
tractography/Images.class:	test/tractography/Images.java
				javac -classpath $(CLASSPATH) test/tractography/Images.java

clean:
				rm -f `find . -name "*~"` `find . -name "*\.class"` `find . -name "tractography/*.Bdouble"`
