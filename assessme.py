#!/usr/bin/env python

# Assess Statistical Estimates from Shotgun Proteomics Scoring Methods: ASSESSMe.py

#   Copyright 2013 Viktor Granholm, Stockholm University
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


import sys
import random
import matplotlib.pyplot as plt

class Documentation(object):
    '''General documentation of the program'''
    def __init__(self, argv):
        self.argv = argv
        self.program_name = argv[0].split('/')[-1]  # Only UNIX systems
        self.entrapment_prefix = 'entrapment_'  # Entrapment proteins always starts with this prefix

    def get_greeting(self):
        return 'Assess Statistical Estimates from Shotgun Proteomics Scoring Methods: %s' % (self.program_name)

    def get_mode_description(self):
        return 'Running in mode: %s' % (self.mode.title())

    def print_mode_help(self):
        '''Print documentation about this mode'''
        help = []
        help.append('')
        help.append('Usage: %s %s [options]' % (self.argv[0], self.mode))
        help.append('')
        help.append('Options:')
        for option in self.mode_options:
            help.append(option.get_documentation_line())
        print '\n'.join(help)


class FastaEntry(object):
    '''Holds fasta entry information'''
    def __init__(self, header, sequence):
        '''Stores a fasta header and a sequence
        Arguments:
        header - A string with everything except initial '>' of a fasta header
        sequence - A string with a complete fasta entry protein sequence
        '''
        self.header = header
        self.sequence = sequence

    def make_fasta_string(self):
        '''Return a fasta entry string'''
        line_length = 80
        lines = [self.header]
        lines.extend([self.sequence[i:i+line_length] for i in range(0, len(self.sequence), line_length)])
        return '\n'.join(lines) + '\n'


class Option(list):
    '''Extended list-object that holds some information about the command line options'''
    def __init__(self, short_flag, long_flag, description, argument_description='<filename>'):
        '''Initialize an option object, with flags, and description
        Arguments:
        short_flag - String of a short version flag, e.g. '-i'
        long_flag - String of the long version flag, e.g. '--input'
        description - String to describe the option
        Keyword arguments:
        argument_description - String to describe whether flag requires some input (e.g. <filename>)
        '''
        list.__init__(self, [short_flag, long_flag])
        self.description = description
        self.argument_description = argument_description

    def get_documentation_line(self):
        '''Return a small documentation of the option'''
        return '%s\n%s %s\t\t%s' % (self[0], self[1], self.argument_description, self.description)


class DatabaseMode(Documentation):
    '''Create entrapment sequence from a sample database, output bipartite database and a reversed decoy'''
    def __init__(self, argv):
        '''Set up Database mode class
        Argument:
        argv - The list of arguments that invoked the program
        '''
        Documentation.__init__(self, argv)
        self.mode = 'database'
        print self.get_greeting()
        print self.get_mode_description()
        # Default parameters
        self.sample_fasta_path = 'known_proteins.fasta'
        self.target_fasta_path = 'target.bipartite.fasta'
        self.decoy_fasta_path = 'decoy.bipartite.fasta'
        self.entrapment_size = 25  # Number of times to shuffle sample partition to get entrapment
        # Define options
        input_options = Option('-i', '--input', 'Fasta file with sequences of proteins present in the sample')
        target_options = Option('-t', '--target', 'Path to output bipartite target fasta file')
        decoy_options = Option('-d', '--decoy', 'Path to output decoy fasta file')
        entrapment_size_options = Option('-e', '--entrap_size', 'Number of times to shuffle sample database to generate entrapments', '<integer>')
        self.mode_options = [input_options, target_options, decoy_options, entrapment_size_options]
        # Parse options
        if len(argv) < 3:
            self.print_mode_help()
            sys.exit()
        index = 2
        while index < len(argv):
            if argv[index] in input_options:
                index += 1
                self.sample_fasta_path = argv[index]
            elif argv[index] in target_options:
                index += 1
                self.target_fasta_path = argv[index]
            elif argv[index] in decoy_options:
                index += 1
                self.decoy_fasta_path = argv[index]
            elif argv[index] in entrapment_size_options:
                index += 1
                self.entrapment_size = int(argv[index])
            else:
                print 'Error: Unknown parameter %s' % (argv[index])
                self.print_mode_help()
                sys.exit()
            index += 1

    def run(self):
        '''Generate a bipartite target database (and decoy) from a smaller database of known proteins'''
        # Generate and write bipartite target
        sample_sequences = self.import_fasta(self.sample_fasta_path)
        entrapment_sequences = self.shuffle_sequences(sample_sequences, 'shuffle', self.entrapment_size, self.entrapment_prefix)
        bipartite_sequences = sample_sequences + entrapment_sequences
        self.write_fasta(bipartite_sequences, self.target_fasta_path)
        print '\nWrote target bipartite database to %s' % self.target_fasta_path
        # Generate and write reversed decoy
        decoy_sequences = self.shuffle_sequences(bipartite_sequences, 'reverse', 1, 'reverse_')
        self.write_fasta(decoy_sequences, self.decoy_fasta_path)
        print 'Wrote reversed decoy database to %s' % self.decoy_fasta_path
        
    def import_fasta(self, path):
        '''Open a fasta-file, return entries in a list of FastaEntry objects
        Arguments:
        path - string with path to fasta-file to read
        '''
        sequences = []
        sequence = ''
        for line in open(path):
            if line.startswith('>'):
                if not sequence == '':
                    sequences.append(FastaEntry(header, sequence))
                sequence = ''
                header = line.strip()[1:]
            else:
                sequence = sequence + line.strip()
        # Don't forget last line
        sequences.append(FastaEntry(header, sequence))
        return sequences

    def shuffle_sequences(self, sequences, method, repeat_number, prefix):
        '''Take a list of FastaEntry objects, repeatedly shuffle sequences and output new, longer list
        Arguments:
        sequence - list of FastaEntry objects
        method - string, either 'shuffle' or 'reverse'
        repeat_number - integer of the number of multiples of shuffled sequence to produce
        prefix - the prefix to give the shuffled sequences
        '''
        shuffled_sequences = []
        sequence_count = 1
        for repeat in range(repeat_number):
            for entry in sequences:
                if method == 'shuffle':
                    shuffled_sequence = self.shuffle(entry.sequence)
                elif method == 'reverse':
                    shuffled_sequence = self.reverse(entry.sequence)
                else:
                    raise Exception('shuffle_sequence does not recognize method: %s' % (method))
                header = prefix + str(sequence_count).zfill(6)
                shuffled_sequences.append(FastaEntry(header, shuffled_sequence))
                sequence_count += 1
        return shuffled_sequences

    def reverse(self, sequence):
        '''Take a string, and reverse it'''
        sequence = list(sequence)
        sequence.reverse()
        return ''.join(sequence)

    def shuffle(self, sequence):
        '''Take a string, and shuffle it. This method is taken from PepHype.py'''
        l = list(sequence)
        # i ranges from the last index to zero, decrementing by 1
        for i in range(len(l)-1,0,-1) :
            # j (swap index) goes from 0 to i
            j = int(random.random() * (i+1))
            # Swap indices i,j
            (l[i],l[j]) = (l[j],l[i])
        return ''.join(l)

    def write_fasta(self, sequences, path):
        '''Take a list of FastaEntry objects, write them to a fasta'''
        outfile = open(path, 'w')
        for entry in sequences:
            outfile.write(entry.make_fasta_string())
        outfile.close()


class CalibrationMode(Documentation):
    '''Checks the uniformity of null p-values, to evaluate their calibration'''
    def __init__(self, argv):
        '''Set up class to perform calibration
        Argument:
        argv - The list of arguments that invoked the program
        '''
        Documentation.__init__(self, argv)
        self.mode = 'calibration'
        print self.get_greeting()
        print self.get_mode_description()
        # Default parameters
        self.identification_path = 'known_proteins.fasta'
        self.figure_path = 'qqplot.png'
        # Define options
        input_options = Option('-i', '--input', 'File with a pvalue and a protein ID on each line')
        plot_options = Option('-p', '--plot_path', 'Path to output Q-Q plot')
        self.mode_options = [input_options, plot_options]
        # Parse options
        if len(argv) < 3:
            self.print_mode_help()
            sys.exit()
        index = 2
        while index < len(argv):
            if argv[index] in input_options:
                index += 1
                self.identification_path = argv[index]
            elif argv[index] in plot_options:
                index += 1
                self.figure_path = argv[index]
            else:
                print 'Error: Unknown parameter %s' % (argv[index])
                self.print_mode_help()
                sys.exit()
            index += 1

    def run(self):
        '''Assess the calibration of null (entrapment) pvalues'''
        pvalues = self.import_pvalues(self.identification_path, self.entrapment_prefix)
        ideal_pvalues = self.get_uniform_distribution(len(pvalues))
        self.make_qqplot(pvalues, ideal_pvalues, self.figure_path)
        dvalue = self.run_kstest(pvalues, ideal_pvalues)
        self.make_calibration_statement(dvalue, self.figure_path)
        
    def import_pvalues(self, filepath, protein_prefix):
        '''Take a file with pvalues and protein IDs, output list of entrapment pvalues
        Arguments:
        filepath - string path to file with a p-value and a protein ID on each line
        protein_prefix - string, only store p-values with protein ID's starting with protein_prefix
        '''
        pvalues = []
        for line in open(filepath):
            words = line.split()
            pvalue = float(words[0])
            protein_id = words[1]
            if protein_id.startswith(protein_prefix):
                pvalues.append(pvalue)
        return pvalues

    def make_qqplot(self, reported_pvalues, ideal_pvalues, figure_path):
        '''Make a log-scale quantile-quantile plot of pvalues, save to filepath
        Arguments:
        reported_pvalues - list of pvalues reported from statistical estimation method
        ideal_pvalues - list of ideal (uniform) pvalues
        figure_path - string to path to print figure
        '''
        reported_pvalues.sort()
        lower_limit = min(reported_pvalues+ideal_pvalues)*0.1
        # Plot
        plt.scatter(ideal_pvalues, reported_pvalues, s=30, c='red', edgecolor='red', marker='o')
        line_a = line_b = [lower_limit, 10.0]
        over_line = [i*2 for i in line_a]
        lower_line = [i/2 for i in line_a]
        plt.plot(line_a, line_b, c='black')
        plt.plot(line_a, over_line, c='grey')
        plt.plot(line_a, lower_line, c='grey')
        plt.xlim([lower_limit, 1])
        plt.ylim([lower_limit, 1])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Ideal $p$ values', fontsize='x-large')
        plt.ylabel('Reported $p$ values', fontsize='x-large')
        plt.savefig(figure_path)

    def run_kstest(self, a, b):
        '''Kolmogorov-Smirnov tests between two distributions
        Arguments:
        a - list of values in first distribution
        b - list of values in second distribution
        '''
        a = [(value, 1) for value in a]
        b = [(value, 0) for value in b]
        all_values = sorted(a + b)
        sorted_indicators = [t[1] for t in all_values]
        max_dvalue = 0
        total_a_count = len(a)
        total_b_count = len(b)
        a_count = b_count = 0.0
        for indicator in sorted_indicators:
            a_count += indicator
            b_count += (1-indicator)
            a_ratio = a_count/total_a_count
            b_ratio = b_count/total_b_count
            max_dvalue = max(max_dvalue, abs(a_ratio-b_ratio))
        return max_dvalue

    def get_uniform_distribution(self, length):
        '''Make an ideal uniform distribution between 0 and 1 of given length (no value = 0)
        Arguments:
        length - the number of values in the distribution
        '''
        interval = 1.0/length
        value = interval
        values = []
        while len(values) < length:
            values.append(value)
            value += interval
        return values

    def make_calibration_statement(self, dvalue, figure_path):
        '''Print a statement about the calibration
        Arguments:
        dvalue - float of KS-test D-value from calibration
        figure_path - string with path to outputted Q-Q plot
        '''
        statement = ['']
        statement.append('A Q-Q plot of the calibration has been printed to: %s' % (figure_path))
        statement.append('The more closely the red point lie to the diagonal, the better calibrated are the p-values')
        statement.append('')
        statement.append('A Kolmogorov-Smirnov test estimated the distance between estimated and ideal p-values to be:')
        statement.append('D-value: %s' % (dvalue))
        statement.append('This number should be as low as possible, preferably less than 0.1')
        print '\n'.join(statement)

class NoMode(Documentation):
    def __init__(self, argv):
        '''Set up class when no mode is choosen
        Argument:
        argv - The list of arguments that invoked the program
        '''
        Documentation.__init__(self, argv)
        self.mode = 'no mode'

    def run(self):
        '''Print the most basic documentation about the program and exits'''
        help = []
        help.append(self.get_greeting())
        help.append('')
        help.append('Usage: %s mode [options]' % (self.argv[0]))
        help.append('')
        help.append('Select a mode for more information:')
        help.append('database (or d):\tGenerate bipartite database')
        help.append('calibration (or c):\tTest the statistical calibration of p-values')
        print '\n'.join(help)
        sys.exit()


def main():
    '''
    ASSESSMe: Assess Statistical Estimates from Shotgun Proteomics Scoring Methods
    '''
    # Read through first arguments and decide on in which mode to run
    if len(sys.argv) < 2:
        no_mode = NoMode(sys.argv)
        no_mode.run()
    # Run in database mode
    elif sys.argv[1] in ['d', 'database']:
        database = DatabaseMode(sys.argv)
        database.run()
    # Run in calibration mode
    elif sys.argv[1] in ['c', 'calibration']:
        calibration = CalibrationMode(sys.argv)
        calibration.run()
    # If unrecognized parameter, run in no-mode mode
    else:
        no_mode = NoMode(sys.argv)
        no_mode.run()

if __name__ == "__main__":
    main()
