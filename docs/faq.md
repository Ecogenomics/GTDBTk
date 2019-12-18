# FAQ

**Q:** Why is there a discrepancy in the naming system between GTDB-Tk and NCBI or Silva taxonomic names?  
**A:** GTDB-Tk uses the GTDB taxonomy ([http://gtdb.ecogenomic.org/](http://gtdb.ecogenomic.org/)). 
This taxonomy is similar, but not identical to NCBI and Silva. 
In many cases the GTDB taxonomy more strictly follows the nomenclatural rules for rank suffixes which is why there is Nitrospirota instead of Nitrospirae.

**Q:** GTDB-Tk crashes after reaching the job memory limit on my Job Scheduler / why does pplacer fail?    
**A:** We believe that this is an issue with how Linux reports memory usage when `pplacer` is run with multiple CPUs. 
It appears that Linux believes that each thread is consuming the total amount of memory allocated, not sharing it. We have typically oberseved this issue affecting queuing systems/HPCs.

A solution is to run the classify step (and thus pplacer) with a single CPU until pplacer can be run using an alternate thread count (see [#195](https://github.com/Ecogenomics/GTDBTk/issues/195) to track the status of that feature).
