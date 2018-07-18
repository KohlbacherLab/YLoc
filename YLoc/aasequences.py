import types;

class AASequences(object):
	
	"Class contains a list of amino acid sequences, each sequence consists of a name and aa sequence."

	def __init__(self, sequences = None):
		self.__sequence_list = sequences;
		if (sequences == None):
			self.__sequence_list = [];

	def __add__(self, right):
		return AASequences(self.__sequence_list + right.__sequence_list);

	def get(self,nr):
		return self.__sequence_list[nr];

	def set(self,sequence_tuple, nr):
		if (types.TupleType == type(sequence_tuple)):
			self.__sequence_list[nr] = sequence_tuple;
		else:
			raise "AASequences: set: argument is not a tuple";
	
	def size(self):
		return len(self.__sequence_list);

	def clear(self):
		self.__sequence_list = [];

	def append(self, sequence_tuple):
		if (types.TupleType == type(sequence_tuple)):
			self.__sequence_list.append(sequence_tuple);
		else:
			raise "AASequences: append: argument is not a tuple";

	def write_fasta_file(self,ff_name):
		ff = open(ff_name,'w');
		for elem in self.__sequence_list:
			(name, seq) = elem;
			ff.write('>'+str(name)+"\n");
			ff.write(seq+"\n");
		ff.close();
			
	def read_fasta_file(self,ff_name):
		ff = open(ff_name);
		entry = self.__getFastaEntry(ff);
		while entry:
			self.append(entry);
			entry = self.__getFastaEntry(ff);

	def __getFastaEntry(self,fh):
		header = fh.readline()

		# eof detection
		if not header:
			return header

		# no fasta format
		if header[0] != '>':
			return None

		header = header.strip();
		seq = ""
		line = fh.readline()
		while line:
			if line[0] == '>':
				# go back to the start of the header line
				fh.seek(-len(line), 1)
				break
			seq += line[:-1]
			line = fh.readline()
        
		return header[1:], seq.upper()
	
		
