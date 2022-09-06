
# we don't need to build the python version of ssq
TARGET_APP=ssq-cxx 

all::
	chmod +x ssq-java ssq-python

ssq-cxx: codingssqs.o

