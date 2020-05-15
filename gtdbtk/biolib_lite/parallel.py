###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import logging
import multiprocessing as mp
import sys
import traceback

from .seq_io import read_seq


class Parallel(object):
    """Processes data in parallel.

    This class implements a producer/consumer model for
    processing data in parallel. In order to keep data
    handling efficient the main process acts as the
    consumer. This keeps all consumed data in the
    address space of the main process.

    Producer and consumer functions are passed to
    the run() method. The return value of the producer
    function is provided to the consumer as an argument.
    The producer function must not return None. If a generic
    return value is required simply return True. A progress
    function can also be specified which allows a progress
    string to be defined.

    Example:

        class Test(object):
            def _producer(self, data):
                return data * data

            def _consumer(self, produced_data, consumer_data):
                if consumer_data == None:
                    # setup structure for consumed data
                    consumer_data = 0

                consumer_data += produced_data

                return consumer_data

            def _progress(self, processed_items, total_items):
                return 'Processed %d of %d items.' % (processed_items, total_items)

            def run(self):
                parallel = Parallel(cpus = 2)
                consumer_data = parallel.run(self._producer,
                                                self._consumer,
                                                data_items = [1, 2, 3, 4, 5],
                                                self._progress)

                print consumer_data
    """

    def __init__(self, cpus=1):
        """Initialization.

        cpus : int
            Number of processes to create.

        """

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def __producer(self, producer_callback, producer_queue, consumer_queue):
        """Process data items with producer processes.

        Parameters
        ----------
        producer_callback : function
            Function to process data items.
        producer_queue : mp.Queue
            Data items to process.
        consumer_queue : mp.Queue
            Queue for holding processed data items to be consumed.
        """
        while True:
            dataItem = producer_queue.get(block=True, timeout=None)
            if dataItem is None:
                break

            rtn = producer_callback(dataItem)

            consumer_queue.put(rtn)

    def __process_manager(self, producer_callback, producer_queue, consumer_queue):
        """Setup producer processes and manage completion of processes.

        Parameters
        ----------
        producer_callback : function
            Function to process data items.
        producer_queue : queue
            Data items to process.
        consumer_queue : mp.Queue
            Queue for holding processed data items to be consumed.
        """

        try:
            producer_proc = [mp.Process(target=self.__producer, args=(
                producer_callback, producer_queue, consumer_queue)) for _ in range(self.cpus)]

            for p in producer_proc:
                p.start()

            for p in producer_proc:
                p.join()

            consumer_queue.put(None)
        except:
            for p in producer_proc:
                p.terminate()

    def run(self, producer, consumer, data_items, progress=None):
        """Setup producers and consumer processes.

        The producer function must be specified and must
        not return None. Consumer and progress can be set to None.

        Parameters
        ----------
        producer : function
            Function to process data items.
        consumer : function
            Function to consumed processed data items.
        data_items : list
            Data items to process.
        progress : function
            Function to report progress string.

        Returns
        -------
        <user specified>
            Set by caller in the consumer function.
        """

        # populate producer queue with data to process
        producer_queue = mp.Queue()
        for d in data_items:
            producer_queue.put(d)

        for _ in range(self.cpus):
            producer_queue.put(None)  # signal processes to terminate

        try:
            consumer_queue = mp.Queue()
            manager_proc = mp.Process(target=self.__process_manager, args=(
                producer, producer_queue, consumer_queue))

            manager_proc.start()

            # process items produced by workers
            items_processed = 0
            consumer_data = None
            while True:
                if progress:
                    status = progress(items_processed, len(data_items))
                    sys.stdout.write('\r%s' % status)
                    sys.stdout.flush()

                produced_data = consumer_queue.get(block=True, timeout=None)
                if produced_data is None:
                    break

                if consumer:
                    consumer_data = consumer(produced_data, consumer_data)

                items_processed += 1

            if progress:
                sys.stdout.write('\n')

            manager_proc.join()

            return consumer_data
        except Exception as _err:
            print(sys.exc_info()[0])
            print(traceback.format_exc())
            self.logger.warning('Exception encountered while processing data.')
            manager_proc.terminate()

    def run_seqs_file(self, producer, consumer, seq_file, progress=None):
        """Process sequences in parallel.

        The producer function must be specified and must
        not return None. Consumer and progress can be set to None.

        Parameters
        ----------
        producer : function
            Function to process data items.
        consumer : queue
            Function to consumed processed data items.
        seq_file : str
            Name of fasta/q file to read.
        progress : function
            Function to report progress string.

        Returns
        -------
        <user specified>
            Set by caller in the consumer function.
        """

        # populate producer queue with data to process
        seq_iter = read_seq(seq_file)
        producer_queue = mp.Queue()
        read_all_seqs = False
        for _ in range(self.cpus):
            try:
                seq_data = next(seq_iter)
                producer_queue.put(seq_data)
            except StopIteration:
                read_all_seqs = True
                for _ in range(self.cpus):
                    producer_queue.put(None)  # signal processes to terminate
                break

        data_items = sum(1 for _ in read_seq(seq_file))
        try:
            consumer_queue = mp.Queue()
            manager_proc = mp.Process(target=self.__process_manager, args=(
                producer, producer_queue, consumer_queue))

            manager_proc.start()

            # process items produced by workers
            items_processed = 0
            consumer_data = None
            while True:
                if progress:
                    status = progress(items_processed, data_items)
                    sys.stdout.write('\r%s' % status)
                    sys.stdout.flush()

                produced_data = consumer_queue.get(block=True, timeout=None)
                if produced_data is None:
                    break

                if not read_all_seqs:
                    try:
                        seq_data = next(seq_iter)
                        producer_queue.put(seq_data)
                    except StopIteration:
                        read_all_seqs = True
                        for _ in range(self.cpus):
                            # signal processes to terminate
                            producer_queue.put(None)

                if consumer:
                    consumer_data = consumer(produced_data, consumer_data)

                items_processed += 1

            if progress:
                sys.stdout.write('\n')

            manager_proc.join()

            return consumer_data
        except Exception as _err:
            print(sys.exc_info()[0])
            print(traceback.format_exc())
            self.logger.warning('Exception encountered while processing data.')
            manager_proc.terminate()
