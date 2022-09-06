Starting Code for lga-coding-ssqs Programming Questions
=======================================================

For using C++, run make(1) or gmake(1) at the command line, make changes to 
codingssqs.cxx.  The resulting simulation is ssq-cxx

For Java, run make(1) or gmake(1) at the command line, make changes to 
codingssqs.java.  The resulting simulation is ssq-java

For Python3, simply edit ssq-python; which is itself the simulation.


All the programs simulate a single-server FIFO service node using
arrival times and service times read from a text file (provided as the
first command line argument).  The server is assumed to be idle when the
first job arrives.  All jobs are processed completely so that the server
is again idle at the end of the simulation.   The output statistics are
the average interarrival time, average service time, the average delay
in the queue, and the average wait (sojourn time) in the service node. 

Additionally, all the programs are faithful representations of Algorithm
1.2.1 of the text, complete with the very tersely named variables.  The
only substantial change is the reporting of summary statistics instead
of the n delay times. Consult the text for details.

