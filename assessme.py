#!/usr/bin/env python

import sys
import random
import matplotlib.pyplot as plt

class Documentation(object):
    '''General documentation of the program'''
    def __init__(self, argv):
        self.argv = argv
        self.program_name = argv[0].split('/')[-1]
        self.greeting = 'Assess Statistical Estimates from Shotgun Proteomics Scoring Methods: %s' % (self.program_name)
        self.entrapment_prefix = 'entrapment_'

    def print_introductory_help(self):
        '''Print the most basic documentation about the program'''
        help = []
        help.append(self.greeting)
        help.append('')
        help.append('Usage: %s mode [options]' % (self.argv[0]))
        help.append('')
        help.append('Select a mode for more information:')
        help.append('database (or d):\tGenerate bipartite database')
        help.append('calibration (or c):\tTest the statistical calibration of p-values')
        print '\n'.join(help)

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
        self.header = header
        self.sequence = sequence

    def make_fasta_string(self):
        '''Return a fasta entry'''
        line_length = 80
        lines = [self.header]
        lines.extend([self.sequence[i:i+line_length] for i in range(0, len(self.sequence), line_length)])
        return '\n'.join(lines) + '\n'


class Option(list):
    def __init__(self, short_flag, long_flag, description, argument_description='<filename>'):
        list.__init__(self, [short_flag, long_flag])
        self.description = description
        self.argument_description = argument_description

    def get_documentation_line(self):
        return '%s\n%s %s\t\t%s' % (self[0], self[1], self.argument_description, self.description)


class Database(Documentation):
    '''Produces bipartite databases'''
    def __init__(self, argv):
        '''Set up Database generation class
        Argument:
        argv - The list of arguments that invoked the program
        '''
        Documentation.__init__(self, argv)
        self.mode = 'database'
        print self.greeting
        print 'Running in mode: Database'
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
        '''Generate a bipartite database from a smaller database of known proteins'''
        # Generate and write bipartite target
        sample_sequences = self.import_fasta(self.sample_fasta_path)
        entrapment_sequences = self.shuffle_sequences(sample_sequences, self.entrapment_size, self.entrapment_prefix)
        bipartite_sequences = sample_sequences + entrapment_sequences
        self.write_fasta(bipartite_sequences, self.target_fasta_path)
        # Generate and write reversed decoy
        decoy_sequences = self.reverse_sequences(bipartite_sequences, 'reverse_')
        self.write_fasta(decoy_sequences, self.decoy_fasta_path)
        
    def import_fasta(self, path):
        '''Open a fasta-file, store entries in a list of FastaEntry objects'''
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

    def shuffle_sequences(self, sequences, repeat_number, prefix):
        '''Take a list of FastaEntry objects and shuffle it into a new list'''
        shuffled_sequences = []
        sequence_count = 1
        for repeat in range(repeat_number):
            for entry in sequences:
                shuffled_sequence = self.shuffle(entry.sequence)
                header = prefix + str(sequence_count).zfill(6)
                shuffled_sequences.append(FastaEntry(header, shuffled_sequence))
                sequence_count += 1
        return shuffled_sequences

    def reverse_sequences(self, sequences, prefix):
        '''Take a list of FastaEntry objects and returned their reversed versions'''
        reversed_sequences = []
        sequence_count = 1
        for entry in sequences:
            sequence = list(entry.sequence)
            sequence.reverse()
            sequence = ''.join(sequence)
            header = prefix + str(sequence_count).zfill(6)
            reversed_sequences.append(FastaEntry(header, sequence))
        return reversed_sequences

    def write_fasta(self, sequences, path):
        '''Take a list of FastaEntry objects, write them to a fasta'''
        outfile = open(path, 'w')
        for entry in sequences:
            outfile.write(entry.make_fasta_string())
        outfile.close()

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


class Calibration(Documentation):
    '''Checks the calibration of pvalues'''
    def __init__(self, argv):
        '''Set up class to perform calibration
        Argument:
        argv - The list of arguments that invoked the program
        '''
        Documentation.__init__(self, argv)
        self.mode = 'calibration'
        print self.get_greeting()
        print 'Running in mode: Calibration'
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
        '''Assess the calibration of pvalues'''
        pvalues = self.import_pvalues(self.identification_path, self.entrapment_prefix)
        ideal_pvalues = self.get_uniform_distribution(len(pvalues))
        self.make_qqplot(pvalues, ideal_pvalues, self.figure_path)
        dvalue = self.run_kstest(pvalues, ideal_pvalues)
        print 'Kolomogorov-Smirnov test D-value estimated to: %s' % (dvalue)
        
    def import_pvalues(self, filepath, protein_prefix):
        '''Take a file with pvalues and protein IDs, output list of entrapment pvalues'''
        pvalues = []
        for line in open(filepath):
            words = line.split()
            pvalue = float(words[0])
            protein_id = words[1]
            if protein_id.startswith(protein_prefix):
                pvalues.append(pvalue)
        pvalues.sort()
        return pvalues

    def make_qqplot(self, reported_pvalues, ideal_pvalues, figure_path):
        '''Make a qq-plot of pvalues, save to filepath'''
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
        print 'Plotted Q-Q plot to %s' % (figure_path)

    def run_kstest(self, distribution_a, distribution_b):
        '''Take a list of null pvalues, test similarity to uniform distribution using KS-test'''
        distribution_a = [(value, 1) for value in distribution_a]
        distribution_b = [(value, 0) for value in distribution_b]
        all_values = sorted(distribution_a + distribution_b)
        sorted_indicators = [t[1] for t in all_values]
        max_dvalue = 0
        total_a = len(distribution_a)
        total_b = len(distribution_b)
        lower_a = 0.0
        lower_b = 0.0
        for indicator in sorted_indicators:
            lower_a += indicator
            lower_b += (1-indicator)
            a_ratio = lower_a/total_a
            b_ratio = lower_b/total_b
            dvalue = abs(a_ratio-b_ratio)
            max_dvalue = max(max_dvalue, dvalue)
        return max_dvalue

    def get_uniform_distribution(self, length):
        '''Make an ideal uniform distribution of given length'''
        # Get an ideal, uniform, null pvalue distribution of same length as query distribution
        interval = 1.0/length
        value = interval
        values = []
        while len(values) < length:
            values.append(value)
            value += interval
        return values


class NoMode(Documentation):
    def __init__(self, argv):
        '''Set up class when no mode is choosen
        Argument:
        argv - The list of arguments that invoked the program
        '''
        Documentation.__init__(self, argv)
        self.mode = 'no mode'

    def run(self):
        self.print_introductory_help()
        sys.exit()


def main():
    '''
    ASSESMe: Assess Statistical Estimates from Shotgun Proteomics Methods
    '''
    if len(sys.argv) < 2:
        no_mode = NoMode(sys.argv)
        no_mode.run()
    # Run in database mode
    elif sys.argv[1] in ['d', 'database']:
        database = Database(sys.argv)
        database.run()
    # Run in calibration mode
    elif sys.argv[1] in ['c', 'calibration']:
        calibration = Calibration(sys.argv)
        calibration.run()
    # If unrecognized parameter, run in no-mode mode
    else:
        no_mode = NoMode(sys.argv)
        no_mode.run()

if __name__ == "__main__":
    main()
