# FAQ

**Q:** Why is there a discrepancy in the naming system between GTDB-Tk and NCBI or Silva taxonomic names?  
**A:** GTDB-Tk uses the GTDB taxonomy ([http://gtdb.ecogenomic.org/](http://gtdb.ecogenomic.org/)). 
This taxonomy is similar, but not identical to NCBI and Silva. 
In many cases the GTDB taxonomy more strictly follows the nomenclatural rules for rank suffixes which is why there is Nitrospirota instead of Nitrospirae.

**Q:** GTDB-Tk crashes after reaching the job memory limit on my Job Scheduler?    
**A:** We believe this is an issue with how Linux reports memory usage when `pplacer` is run with multiple CPUs. 
It appears that Linux believes the amount of memory being requested is factored by the number of CPUs. As such, 
this can cause issues with queuing systems. It will naturally run slower, but a solution is to use a single CPU.

**Q:** How is the translation table automatically determined? <br>
**A:** See [this page](prodigal.md#translation-table-selection)
