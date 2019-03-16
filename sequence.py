import re


class Sequence:
    """A class that represents DNA and RNA sequences and allows for comparisons and calculations on them.
    There is a method for reading files of single lined sequences with a single lined header.

     Attributes:
         __sequence: A string representing the actual sequence of DNA with bases (ATCG), or RNA with bases (AUCG)
         __type: String on format 'DNA' or 'RNA', specifying the input sequence type (casing is ignored).
         __is_valid: A boolean automatically set to True if the sequence is valid, that is if:
             (i)__type = 'DNA' and __sequence contains only the characters A, T, C, G
             (ii)__type = 'RNA' and __sequence contains only the characters A, U, C, G
             Otherwise __is_valid = False
        __T_or_U: String automatically set to 'T' if __type = 'DNA', or 'U' if __type = 'RNA'

     """

    def __init__(self, __sequence, __type='DNA'):
        """ Returns a Sequence object with sequence *__sequence* (may contain only characters ATCG for type 'DNA'
        OR only characters AUCG for type 'RNA') and type *__type*, default is 'DNA'. """
        # upper() makes sure all characters are uppercase to make the class insensitive to casing in the input
        self.__sequence = __sequence.upper()
        self.__type = __type.upper()
        if self.__type not in ('DNA', 'RNA'):
            raise Exception('Invalid type specified, acceptable types are: DNA(default) or RNA.')
        if self.__type == 'RNA':
            self.__T_or_U = 'U'
        else:
            self.__T_or_U = 'T'
        self.__is_valid = self.__is_valid()

    def __is_valid(self):
        """Returns Boolean to field *__is_valid* which is False if the object has incorrect format
        (contains incorrect characters for the specified type), or True otherwise"""
        if self.__type == 'DNA':
            for i in range(self.n_bases()):
                if self.__sequence[i] not in ('A', 'T', 'C', 'G'):
                    return False
        elif self.__type == 'RNA':
            for i in range(self.n_bases()):
                if self.__sequence[i] not in ('A', 'U', 'C', 'G'):
                    return False
        return True

    def get_type(self):
        """Returns value of private field *___type*(string)"""
        return self.__type

    def get_is_valid(self):
        """Returns value of private field *___is_valid*(Boolean)"""
        return self.__is_valid

    def get_sequence(self):
        """Returns value of private field *___sequence*(string)"""
        return self.__sequence

    def n_bases(self):
        """Computes and returns the number of bases(characters in *__sequence*) as an integer"""
        return len(self.__sequence)

    def complement(self, output_type='DNA'):
        """Returns new Sequence object that is complementary to *self*.
        Arg *output_type* specifies the sequence type of the output Sequence object, default is 'DNA'"""
        if not self.__is_valid:
            raise Exception('cannot find complement to invalid sequence')
        t_or_u = 'T'
        if output_type.upper() == 'RNA':
            t_or_u = 'U'
        complement = ''
        for i in range(self.n_bases()):
            if self.__sequence[i] == 'A':
                complement += t_or_u
            if self.__sequence[i] == self.__T_or_U:
                complement += 'A'
            if self.__sequence[i] == 'C':
                complement += 'G'
            if self.__sequence[i] == 'G':
                complement += 'C'
        return Sequence(complement, output_type.upper())

    def first_unmatched_basepair(self, other):
        """Locates and returns the index of the first unmatched pair of bases between two Sequence objects
        of equal type and length, returns -1 if the sequences are identical.
        Arg *other* other Sequence object to be compared with. """
        if self.__type != other.__type:
            raise Exception('cannot compare sequences of different type')
        if self.n_bases() != other.n_bases():
            raise Exception('cannot compare sequences of different lengths')
        for i in range(self.n_bases()):
            if self.__sequence[i] != other.__sequence[i]:
                return i
        return -1

    def gene_separation(self):
        """Splits a Sequence object into a list of Sequence objects based on occurrences of 10 A followed by 10 T"""
        genes = re.split('A{10}T{10}', self.__sequence, flags=re.IGNORECASE)
        # In the cases of a sequence beginning, ending or has two repeated sequences of 10 A and 10 B, re.split() will
        # return empty strings at those places. Filter removes all empty strings from the list but returns
        # a filter object, and needs to be converted back to a list object.
        genes = list(filter(None, genes))
        gene_sequences = []
        for i in range(len(genes)):
            gene_sequences.append(Sequence(genes[i]))
        return gene_sequences


    def sequence_from_file(input_file):
        """Reads a two line file and returns a Sequence object with the sequence
        read in the second line.
        Arg *input_file* file to be read."""
        with open(input_file, encoding='ASCII') as file:
            # skip header
            next(file)
            output_sequence = str(file.read())
            return Sequence(output_sequence)

    def swap_mutation(self, other):
        """Counts and returns the number of swap mutations (differences) between two sequence objects of the same type
        Arg *other* other Sequence to be compared with"""
        if self.__type != other.__type:
            raise Exception('cannot compare sequences of different type')
        if self.n_bases() != other.n_bases():
            raise Exception('cannot compare sequences of different lengths')
        n_swaps = 0
        for i in range(self.n_bases()):
            if self.__sequence[i] != other.__sequence[i]:
                n_swaps = n_swaps + 1
        return n_swaps

    def __eq__(self, other):
        """ Overrides default operator '==' and '!=' for comparison of Sequence objects.
        Arg *other* other Sequence object to be compared with. """
        if isinstance(other, self.__class__) and (self.__sequence == other.__sequence):
            return True
        return False

    def __str__(self):
        """Returns a string representation of the object"""
        return '<' + self.__sequence + '>'
