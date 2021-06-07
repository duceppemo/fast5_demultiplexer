import os
import subprocess
from argparse import ArgumentParser
from multiprocessing import cpu_count
from pathlib import Path
import shutil
import gzip
# from concurrent import futures

# TODO: start with the multi fast5 folder from MinKNOW. Add function to convert multi to single.
# TODO: list reads names form files in parallel
# TODO: convert demultiplexed fast5 from single to multi in parallel

__author__ = 'duceppemo'
__version__ = '0.1'


class Demul(object):
    def __init__(self, args):
        self.args = args
        self.basecalled_folder = args.basecalled
        self.single_fast5_folder = args.singles
        self.demultimplexed_folder = args.demultiplexed
        self.threads = args.threads

        # Run
        self.run()

    def run(self):
        print('Listing single fast5 files...')
        fast5_dict = Demul.get_fast5_tree(self.single_fast5_folder)
        print('Listing demultiplexed fastq files...')
        fastq_dict = Demul.get_fastq_tree(self.basecalled_folder)
        print('Reorganizing fast5 according to barcodes...')
        Demul.move_fast5(fast5_dict, fastq_dict, self.single_fast5_folder, self.demultimplexed_folder)
        print('Converting single fast5 to multi fast5...')
        Demul.single_to_multi_fast5_loop(fastq_dict, self.demultimplexed_folder,
                                         self.demultimplexed_folder + '_multi', self.threads)
        print('Removing temporary folders')

    @staticmethod
    def get_fast5_tree(single_fasta5_folder):
        fast5_dict = dict()
        for root, dirs, files in os.walk(single_fasta5_folder):
            for name in files:
                if name.endswith('.fast5'):
                    read_name = '.'.join(name.split('.')[:-1])
                    fast5_dict[read_name] = os.path.join(root, name)

        return fast5_dict

    @staticmethod
    def get_fastq_tree(basecalled_folder):
        fastq_dict = dict()

        for root, dirs, files in os.walk(basecalled_folder):
            for name in files:
                if name.endswith(('.fastq', '.fastq.gz')):
                    absolute_path = os.path.join(root, name)
                    path_list = absolute_path.split('/')
                    status = path_list[-3]  # pass or fail
                    barcode = path_list[-2]  # barcode01, ..., unclassified

                    read_list = Demul.extract_read_names(absolute_path)

                    if status not in fastq_dict:
                        fastq_dict[status] = dict()
                        fastq_dict[status][barcode] = read_list
                    elif barcode not in fastq_dict[status]:
                        fastq_dict[status][barcode] = read_list

        return fastq_dict

    @staticmethod
    def extract_read_names(input_fastq):
        name_list = list()
        line_counter = 0
        with gzip.open(input_fastq, 'rb') if input_fastq.endswith('.gz') else open(input_fastq, 'r') as f:
            for line in f:
                line_counter += 1
                if input_fastq.endswith('.gz'):
                    line = line.decode()
                if line_counter == 1:
                    name_list.append(line.split()[0].replace('@', ''))
                elif line_counter == 4:
                    line_counter = 0
        return name_list

    @staticmethod
    def make_folders(folder):
        Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def move_fast5(fast5_dict, fastq_dict, single_fasta5_folder, demultiplexed_fast5_folder):
        # Maker folders
        for status, info_dict in fastq_dict.items():
            for barcode, name_list in info_dict.items():
                Demul.make_folders(os.path.join(demultiplexed_fast5_folder, status, barcode))
                for name in name_list:
                    source = fast5_dict[name]
                    destination = os.path.join(demultiplexed_fast5_folder, status, barcode, name + '.fast5')
                    shutil.move(source, destination)

    @staticmethod
    def multi_to_single_fast5(input_fast5_folder, output_fast5_folder, n_cpu):
        cmd = ['multi_to_single_fast5',
               '--recursive',
               '-i', input_fast5_folder,
               '-s', output_fast5_folder,
               '-t', str(n_cpu)]
        subprocess.run(cmd)

    @staticmethod
    def single_to_multi_fast5(input_fast5_folder, output_fast5_folder, n_cpu):
        cmd = ['single_to_multi_fast5',
               '-i', input_fast5_folder,
               '-s', output_fast5_folder,
               # '-f', os.path.dirname(input_fast5_folder).split('/')[-1],
               '-t', str(n_cpu)]
        subprocess.run(cmd)

    @staticmethod
    def single_to_multi_fast5_loop(fastq_dict, single_fast5_demultuplexed_folder,
                                   multi_fast5_demultuplexed_folder, n_cpu):
        for status, info_dict in fastq_dict.items():
            for barcode, name_list in info_dict.items():
                Demul.make_folders(os.path.join(multi_fast5_demultuplexed_folder, status, barcode))
                Demul.single_to_multi_fast5(os.path.join(single_fast5_demultuplexed_folder, status, barcode),
                                            os.path.join(multi_fast5_demultuplexed_folder, status, barcode),
                                            n_cpu)

    # @staticmethod
    # def single_to_multi_fast5_parallel(input_fast5_folder, out_fast5_folder, n_cpu):
    #     out_folder = os.path.dirname(input_fast5_folder).split('/')[-1]
    #     Demul.make_folders(out_folder)
    #     pass
    #
    #     with futures.ThreadPoolExecutor(max_workers=n_cpu/4) as executor:
    #         args = ((a, b) for a, b in XXX.items())
    #         for results in executor.map(lambda p: Demul.single_to_multi_fast5(*p), args):
    #             pass


if __name__ == '__main__':
    cpu = cpu_count()

    parser = ArgumentParser(description='Reorganize fast5 single files according to Guppy demultiplexing results.'
                                        ' WARNING: about times the disk footprint of the fast5 folder is required!')
    parser.add_argument('-b', '--basecalled', metavar='/basecalled/',
                        required=True,
                        help='"Save_path" folder used for Guppy.')
    parser.add_argument('-s', '--singles', metavar='/single_fast5/',
                        required=True,
                        help='parent folder where all the single fast5 files are located.')
    parser.add_argument('-d', '--demultiplexed', metavar='/demultiplexed_fast5',
                        required=True,
                        help='Folder where to move sinfle fast5 files according to Guppy demultiplexing.')
    parser.add_argument('-t', '--threads', metavar='{}'.format(cpu),
                        required=False, default=cpu,
                        type=int,
                        help='Number of threads to use. Default max available.')

    arguments = parser.parse_args()  # Get the arguments into an object
    Demul(arguments)
