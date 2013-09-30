#!/usr/bin/env python

# nullvalidator.py

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

import os
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
        return '%s\nValidate the calibration of null statistics for shotgun proteomics results' % (self.program_name)

    def get_mode_description(self):
        return 'Running in mode: %s\n' % (self.mode.title())

    def print_mode_help(self):
        '''Print documentation about this mode'''
        help = []
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


class PercolatorXML(object):
    '''Class that holds Percolator XML objects'''
    def __init__(self, filepath):
        from lxml import etree  # Importing here makes the program less dependent
        self.filepath = filepath
        parser = etree.XMLParser(ns_clean=False, huge_tree=False)        
        self.tree = etree.parse(self.filepath, parser)
        self.ns = self.get_namespace()

    def get_namespace(self):
        root_element = self.tree.getroot()
        schema_location = root_element.attrib['{http://www.w3.org/2001/XMLSchema-instance}schemaLocation']
        for format_name in ['percolator_in', 'percolator_out']:
            for ns_number in [11, 12, 13, 14]:
                namespace = 'http://per-colator.com/%s/%s' % (format_name, ns_number)
                if namespace in schema_location:
                    return namespace

    def get_type(self):
        if 'percolator_in' in self.ns:
            return 'percolator_in'
        else:
            return 'percolator_out'

    def get_feature_names(self):
        '''Output a list of all the feature names in the pin file'''
        feature_names = []
        for element in self.tree.findall('{%s}featureDescriptions/{%s}featureDescription' % (self.ns, self.ns)):
            name = element.attrib['name']
            feature_names.append(name)
        return feature_names

    def get_feature_values(self, feature_index, isDecoy=False, type_func=float):
        '''Extract all features values for a given feature, and output a list (first feature has index 0)'''
        values = []
        for element in self.tree.findall('{%s}fragSpectrumScan/{%s}peptideSpectrumMatch' % (self.ns, self.ns)):
            element_isDecoy = element.attrib['isDecoy'] == 'true'  # Convert element string boolean to boolean
            if isDecoy == element_isDecoy:
                feature_element = element.findall('{%s}features/{%s}feature' % (self.ns, self.ns))[feature_index]
                feature_value = feature_element.text
                values.append(type_func(feature_value))  # type_func: float(), int(), etc...
        return values


class PValueList(list):
    '''Extended list object that holds pvalues from a particular set of scores'''
    def __init__(self, title, pvalues=[], target_scores=[], decoy_scores=[]):
        '''Take a title, and either pre-estimates pvalues, or scores to estimate pvalues
        Arguments:
        title - A string title to print in the output
        pvalues - A list of float pvalues, could be empty
        target_scores - A list of float target scores, not needed if pvalues is populated
        decoy_scores - A list of float decoy scores, not needed if pvalues is populated
        '''
        list.__init__(self, pvalues)
        self.title = title
        self.target_scores = target_scores
        self.decoy_scores = decoy_scores
        if len(self) == 0:  # No pvalues were supplied
            self._run_target_decoy_pvalues()

    def _run_target_decoy_pvalues(self):
        '''Run a target-decoy analysis on target and decoy scores'''
        if len(self.target_scores) == 0 or len(self.target_scores) == 0:
            raise Exception("Can't calculate p values without target or decoy scores")
        targets = [(value, True) for value in self.target_scores]
        decoys = [(value, False) for value in self.decoy_scores]
        all_values = sorted(targets + decoys)
        sorted_indicators = [i[1] for i in all_values]  # Lowest value first
        total_decoys = len(decoys)
        higher_decoys = total_decoys
        for is_target in sorted_indicators:
            if is_target:
                pvalue = (higher_decoys+1.0)/(total_decoys+1.0)  # ((r+1)/(n+1))
                self.append(pvalue)
            else:
                higher_decoys -= 1

    def get_qqplot_filename(self):
        '''Produce a string from self.title that could be used as filename for plot'''
        words = [word.lower() for word in self.title.split()]
        return '_'.join(words) + '.png'


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
        self.identification_path = 'identifications.txt'
        self.figure_directory = '.'
        self.html_path = 'calibration.html'
        # Define options
        input_options = Option('-i', '--input', 'File with a pvalue and a protein ID on each line, or pin or pout XML')
        html_options = Option('-o', '--output', 'Path to output HTML file')
        plot_options = Option('-p', '--plot_dir', 'Directory to output Q-Q plot')
        self.mode_options = [input_options, html_options, plot_options]  # Make list of all options, for documentation
        # Parse options
        if len(argv) < 3:
            self.print_mode_help()
            sys.exit()
        index = 2
        while index < len(argv):
            if argv[index] in input_options:
                index += 1
                self.identification_path = argv[index]
            elif argv[index] in html_options:
                index += 1
                self.html_path = argv[index]
            elif argv[index] in plot_options:
                index += 1
                self.figure_directory = argv[index]
            else:
                print 'Error: Unknown parameter %s' % (argv[index])
                self.print_mode_help()
                sys.exit()
            index += 1

    def run(self):
        '''Assess the calibration of null (entrapment) pvalues'''
        pvalue_lists = self.import_pvalues()
        dvalues = []
        figure_paths = []
        for pvalues in pvalue_lists:
            # Get figures and D values
            ideal_pvalues = self.get_uniform_distribution(len(pvalues))
            figure_paths.append(self.make_qqplot(pvalues, ideal_pvalues))
            dvalues.append(self.run_kstest(pvalues, ideal_pvalues))
        # Make output
        titles = [pvalues.title for pvalues in pvalue_lists]
        self.make_html_output(titles, dvalues, figure_paths)
        
    def import_pvalues(self):
        '''Import pvalues or scores from self.identification_path, output list of entrapment PValueList objects'''
        infile = open(self.identification_path)
        first_line = infile.readline()
        infile.close()
        if first_line.strip().startswith('<?xml version'):  # Check first line
            # XML file
            pvalue_lists = self.read_input_xml(self.identification_path, self.entrapment_prefix)
        else:
            # Tab-delimited text file
            pvalue_lists = self.read_tab_file(self.identification_path, self.entrapment_prefix)
        return pvalue_lists

    def read_tab_file(self, filepath, protein_prefix):
        '''Take a path to a file with pvalues and proteins IDs, output list with a PValueList object
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
        pvalue_list = PValueList('Reported pvalues', pvalues)
        return [pvalue_list]  # Return in list to work with subsequent for loops

    def read_input_xml(self, filepath, protein_prefix):
        '''Take a path to a pin or pout XML file, output lists with PValueList object
        Arguments:
        filepath - string path to a pin- or pout-xml file
        protein_prefix - string, only store p-values with protein ID's starting with protein_prefix
        '''
        pvalue_lists = []
        xml = PercolatorXML(filepath)
        if xml.get_type() == 'percolator_in':
            features = xml.get_feature_names()
            for feature in features:
                targets = xml.get_feature_values(feature, is_decoy=False, with_prefix=protein_prefix)
                decoys = xml.get_features_values(feature, is_decoy=True, with_prefix='')
                pvalues = PValueList(feature_name, target_scores=targets, decoy_scores=decoys)
                pvalue_lists.append(pvalues)
        elif xml.get_type() == 'percolator_out':
            pvalues = xml.get_pvalues(is_decoy=False, with_prefix=protein_prefix)
            pvalues = PValue('Percolator pvalues', pvalues)
            pvalue_lists.append(pvalues)
        return pvalue_lists

    def make_qqplot(self, reported_pvalues, ideal_pvalues):
        '''Make a log-scale quantile-quantile plot of pvalues, save to filepath
        Arguments:
        reported_pvalues - a PValueList object of pvalues reported from a statistical estimation method
        ideal_pvalues - a PValueList object of ideal (uniform) pvalues
        '''
        figure_path = os.path.abspath('%s/%s' % (self.figure_directory, reported_pvalues.get_qqplot_filename()))
        reported_pvalues.sort()
        lower_limit = min(reported_pvalues+ideal_pvalues)*0.1
        # Plot
        plt.scatter(ideal_pvalues, reported_pvalues, s=30, c='red', edgecolor='red', marker='o')
        line_a = line_b = [lower_limit, 10.0]
        over_line = [i*2 for i in line_a]
        lower_line = [i/2 for i in line_a]
        plt.plot(line_a, line_b, c='black')
        plt.plot(line_a, over_line, c='grey', linestyle='--')
        plt.plot(line_a, lower_line, c='grey', linestyle='--')
        plt.xlim([lower_limit, 1])
        plt.ylim([lower_limit, 1])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Ideal $p$ values', fontsize='x-large')
        plt.ylabel('Reported $p$ values', fontsize='x-large')
        plt.title(reported_pvalues.title)
        plt.savefig(figure_path)
        return figure_path

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

    def make_html_output(self, titles, dvalues, figure_paths):
        '''Print a statement about the calibration
        Arguments:
        dvalue - float of KS-test D-value from calibration
        figure_path - string with path to outputted Q-Q plot
        '''
        html = open(self.html_path, 'w')
        html.write('<HTML>\n')
        html.write('<HEAD>\n<TITLE>Statistical calibration</TITLE>\n</HEAD>\n')
        html.write('<BODY>\n')
        html.write('<H2>nullvalidator.py automatic output</H2>\n')
        html.write('''
This page shows the results from a validation of the calibration of null statistics for
shotgun proteomics results. For a detailed description of the method, please refer to the publication:<br>
<a href="http://pubs.acs.org/doi/abs/10.1021/pr1012619">"On Using Samples of Known Protein Content to Assess
the Statistical Calibration of Scores Assigned to Peptide-Spectrum Matches in Shotgun Proteomics"</a> by <b>
Granholm <i>et al.</i></b>, <i>Journal of Proteome Research</i>, 2011.<br>''')
        html.write('<H4>Guidelines for interpretation</H4>\n')
        html.write('''
The program nullvalidator.py selects only <i>p</i> values of PSMs (or peptides) that map \
to the entrapment partition of a bipartite database, as these "incorrect" identifications make up a good null
model. If the statistical estimation procedure is accurate, these null <i>p</i> values should follow a uniform
distribution. The Kolmogorov-Smirnov test <i>D</i> value is a measure of how distantly the null <i>p</i> values
are distributed relative the uniform distribution. A low <i>D</i> value indicates high similarity to the uniform
distribution. The quantile-quantile shows the reported null <i>p</i> values plotted against an ideal null
distribution, the closer the points lie to the <i>x=y</i> diagonal, the more uniform they are. The dashed grey
lines represent the lines <i>x=2y</i> and <i>x=y/2</i>. In short:\n''')
        html.write('<UL>\n')
        html.write('<LI>A <i>D</i> value greater than 0.1 is generally not good (poor calibration).\n')
        html.write('<LI>If many points consistently lie outside the dashed grey lines, it\'s a bad sign (poor calibration).\n')
        html.write('</UL>\n')
        html.write('<H3>Summary of input</H3>\n')
        html.write('Path to input file with <i>p</i> values/scores and protein names: <b>%s</b><br>\n' % (self.identification_path))
        html.write('Protein prefix to recognize entrapment identifications: <b>%s</b><br>\n' % (self.entrapment_prefix))
        html.write('<H3>Summary of output</H3>\n')
        for index, title in enumerate(titles):
            dvalue = dvalues[index]
            figure_path = figure_paths[index]
            html.write('<H4>%s</H4>\n' % (title))
            html.write('Kolomogorov-Smirnov test <I>D</I> value: <B>%s</B><BR>\n' % (dvalue))
            html.write('Q-Q plot: <B><A HREF="%s">%s</A></B><BR>\n' % (figure_path, figure_path.split('/')[-1]))  # On UNIX systems
            html.write('Path to Q-Q plot: <B>%s</B><BR>\n' % (figure_path))
        html.write('</BODY>\n')
        html.write('</HTML>\n')
        html.close()
        print 'Wrote results summary to %s, go to this address from a browser' % (self.html_path)
        print 'file://%s' % (os.path.abspath(self.html_path))

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
    nullvalidator, validate the calibration of null statistics for shotgun proteomics results
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
